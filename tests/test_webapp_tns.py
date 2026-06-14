import pytest

from webapp.backend.tns import TNSClient, TNSCredentialsError, target_from_row


def test_tns_client_requires_server_credentials(monkeypatch):
    monkeypatch.delenv("TNS_API_KEY", raising=False)
    monkeypatch.delenv("TNS_BOT_ID", raising=False)
    monkeypatch.delenv("TNS_BOT_NAME", raising=False)

    with pytest.raises(TNSCredentialsError):
        TNSClient().lookup("SN 2023ixf")


def test_target_from_row_parses_tns_payload():
    target = target_from_row(
        {
            "objname": "2023ixf",
            "prefix": "SN",
            "ra": "14:03:38.562",
            "dec": "+54:18:41.94",
            "object_type": {"name": "SN II"},
            "redshift": "0.0008",
            "hostname": "M101",
            "internal_names": "ZTF23abc; ATLAS23abc",
        },
        "2023ixf",
    )

    assert target is not None
    assert target.name == "SN 2023ixf"
    assert target.ra_deg == pytest.approx(210.910675)
    assert target.dec_deg == pytest.approx(54.31165, rel=1e-5)
    assert target.transient_type == "SN II"
    assert target.aliases == ["ZTF23abc", "ATLAS23abc"]


def test_tns_lookup_uses_authenticated_api(monkeypatch):
    monkeypatch.setenv("TNS_API_KEY", "key")
    monkeypatch.setenv("TNS_BOT_ID", "42")
    monkeypatch.setenv("TNS_BOT_NAME", "NightSkyPlan")

    def fake_post_api(self, endpoint, payload):
        if endpoint == "search":
            return {"data": {"reply": [{"objname": "2023ixf"}]}}
        return {
            "data": {
                "reply": {
                    "objname": "2023ixf",
                    "prefix": "SN",
                    "radeg": 210.910675,
                    "decdeg": 54.311651,
                }
            }
        }

    monkeypatch.setattr(TNSClient, "_post_api", fake_post_api)

    target = TNSClient(sleep=0).lookup("SN 2023ixf")

    assert target.name == "SN 2023ixf"
    assert target.ra_deg == 210.910675
    assert target.dec_deg == 54.311651
