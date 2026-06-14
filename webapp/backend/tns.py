from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, field
from typing import Any

import requests

from webapp.backend.models import TargetResolveResponse


TNS_BASE_URL = "https://www.wis-tns.org/api/get"


class TNSLookupError(RuntimeError):
    """Raised when a TNS target cannot be resolved."""


class TNSCredentialsError(TNSLookupError):
    """Raised when required server-side TNS credentials are missing."""


@dataclass(frozen=True)
class ResolvedTarget:
    name: str
    ra_deg: float
    dec_deg: float
    tns_name: str = ""
    prefix: str = ""
    objid: str = ""
    transient_type: str = ""
    redshift: str = ""
    host_name: str = ""
    aliases: list[str] = field(default_factory=list)

    def to_response(self) -> TargetResolveResponse:
        return TargetResolveResponse(
            name=self.name,
            ra_deg=self.ra_deg,
            dec_deg=self.dec_deg,
            tns_name=self.tns_name,
            prefix=self.prefix,
            objid=self.objid,
            transient_type=self.transient_type,
            redshift=self.redshift,
            host_name=self.host_name,
            aliases=self.aliases,
        )


def normalize_query_name(name: str) -> str:
    text = name.strip()
    for prefix in ("SN", "AT"):
        if text.upper().startswith(prefix + " "):
            return text.split(None, 1)[1].strip()
    return text


class TNSClient:
    def __init__(self, timeout: float = 45.0, sleep: float = 1.0) -> None:
        self.timeout = timeout
        self.sleep = sleep
        self.api_key = os.environ.get("TNS_API_KEY", "")
        self.bot_id = os.environ.get("TNS_BOT_ID", "")
        self.bot_name = os.environ.get("TNS_BOT_NAME", "")
        self.tns_type = os.environ.get("TNS_TYPE", "bot")
        self.session = requests.Session()

    @property
    def has_credentials(self) -> bool:
        return bool(self.api_key and self.bot_id and self.bot_name)

    def require_credentials(self) -> None:
        if not self.has_credentials:
            raise TNSCredentialsError("TNS_API_KEY, TNS_BOT_ID, and TNS_BOT_NAME must be set on the server.")

    def lookup(self, name: str) -> ResolvedTarget:
        self.require_credentials()
        clean_name = normalize_query_name(name)
        target = self._lookup_api(clean_name)
        if target is None:
            raise TNSLookupError(f"No TNS target found for {name!r}.")
        return target

    def _marker_header(self) -> dict[str, str]:
        marker = f'tns_marker{{"tns_id": "{self.bot_id}", "type": "{self.tns_type}", "name": "{self.bot_name}"}}'
        return {"User-Agent": marker}

    def _post_api(self, endpoint: str, payload: dict[str, Any]) -> Any:
        response = self.session.post(
            f"{TNS_BASE_URL}/{endpoint}",
            headers=self._marker_header(),
            data={"api_key": self.api_key, "data": json.dumps(payload)},
            timeout=self.timeout,
        )
        if response.status_code >= 400:
            raise TNSLookupError(response.text.replace("\n", " ")[:500])
        time.sleep(self.sleep)
        return response.json()

    @staticmethod
    def _reply(payload: Any) -> Any:
        if isinstance(payload, list):
            payload = payload[0] if payload else {}
        if not isinstance(payload, dict):
            return {}
        data = payload.get("data", {})
        return data.get("reply", data) if isinstance(data, dict) else {}

    def _lookup_api(self, clean_name: str) -> ResolvedTarget | None:
        if clean_name.lower().startswith("ztf"):
            search_payload = {"internal_name": clean_name, "num_page": 10}
        else:
            search_payload = {"objname": clean_name, "num_page": 10}

        reply = self._reply(self._post_api("search", search_payload))
        rows = reply if isinstance(reply, list) else reply.get("objects", []) if isinstance(reply, dict) else []
        candidates = [clean_name]
        for row in rows:
            row_name = _scalar(_first(row, ("objname", "name"), ""))
            if row_name and row_name not in candidates:
                candidates.append(row_name)

        for candidate in candidates:
            try:
                obj_reply = self._reply(
                    self._post_api("object", {"objname": candidate, "photometry": "0", "spectra": "0"})
                )
            except TNSLookupError:
                continue
            if isinstance(obj_reply, dict):
                target = target_from_row(obj_reply, candidate)
                if target is not None:
                    return target

        for row in rows:
            target = target_from_row(row, clean_name)
            if target is not None:
                return target
        return None


def target_from_row(row: dict[str, Any], fallback_name: str) -> ResolvedTarget | None:
    raw_name = _scalar(_first(row, ("objname", "name", "Name", "TNS Name"), fallback_name))
    prefix = _scalar(_first(row, ("prefix", "Prefix", "name_prefix"), ""))
    if raw_name.upper().startswith("SN "):
        prefix, raw_name = "SN", raw_name.split(None, 1)[1]
    elif raw_name.upper().startswith("AT "):
        prefix, raw_name = "AT", raw_name.split(None, 1)[1]

    ra = _parse_ra(_first(row, ("radeg", "ra", "RA", "ra_deg", "objra")))
    dec = _parse_dec(_first(row, ("decdeg", "dec", "DEC", "declination", "objdec")))
    if ra is None or dec is None:
        return None

    aliases = _scalar(_first(row, ("internal_names", "internal_name", "aliases"), ""))
    display = f"{prefix} {raw_name}".strip()
    return ResolvedTarget(
        name=display,
        tns_name=raw_name,
        prefix=prefix,
        objid=_scalar(_first(row, ("objid", "id", "ID"), "")),
        ra_deg=ra,
        dec_deg=dec,
        transient_type=_scalar(_first(row, ("object_type.name", "objtype.name", "type", "Type"), "")),
        redshift=_scalar(_first(row, ("redshift", "Redshift", "z"), "")),
        host_name=_scalar(_first(row, ("hostname", "Host Name", "host.name", "host"), "")),
        aliases=[part.strip() for part in aliases.replace(",", ";").split(";") if part.strip()],
    )


def _first(row: dict[str, Any], keys: tuple[str, ...], default: Any = "") -> Any:
    for key in keys:
        if key in row and row[key] not in (None, "", [], {}):
            return row[key]
        current: Any = row
        for part in key.split("."):
            if not isinstance(current, dict) or part not in current:
                current = None
                break
            current = current[part]
        if current not in (None, "", [], {}):
            return current
    return default


def _scalar(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, dict):
        return str(_first(value, ("name", "value", "id"), ""))
    if isinstance(value, list):
        return "; ".join(part for item in value if (part := _scalar(item)))
    return str(value).strip()


def _parse_ra(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        pass
    parts = text.replace("h", ":").replace("m", ":").replace("s", "").split(":")
    if len(parts) < 3:
        parts = text.split()
    if len(parts) < 3:
        return None
    try:
        return (float(parts[0]) + float(parts[1]) / 60.0 + float(parts[2]) / 3600.0) * 15.0
    except ValueError:
        return None


def _parse_dec(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        pass
    sign = -1.0 if text.startswith("-") else 1.0
    parts = text.lstrip("+-").replace("d", ":").replace("m", ":").replace("s", "").split(":")
    if len(parts) < 3:
        parts = text.lstrip("+-").split()
    if len(parts) < 3:
        return None
    try:
        return sign * (float(parts[0]) + float(parts[1]) / 60.0 + float(parts[2]) / 3600.0)
    except ValueError:
        return None
