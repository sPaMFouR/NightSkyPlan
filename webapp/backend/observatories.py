from __future__ import annotations

from dataclasses import asdict

from nightskyplan.models import Observatory


OBSERVATORIES: dict[str, Observatory] = {
    "HCT": Observatory(
        name="HCT / Hanle",
        latitude_deg=32.7794,
        longitude_deg=78.9642,
        elevation_m=4486.0,
        timezone="Asia/Kolkata",
        horizon_deg=25.0,
        zenith_deg=85.0,
    ),
    "DOT": Observatory(
        name="DOT / Devasthal",
        latitude_deg=29.3608,
        longitude_deg=79.6844,
        elevation_m=2450.0,
        timezone="Asia/Kolkata",
        horizon_deg=15.0,
        zenith_deg=87.5,
    ),
    "KANATA": Observatory(
        name="Kanata / Higashi-Hiroshima",
        latitude_deg=34.3764,
        longitude_deg=132.7767,
        elevation_m=503.0,
        timezone="Asia/Tokyo",
        horizon_deg=10.0,
        zenith_deg=88.0,
    ),
    "LAPALMA": Observatory(
        name="Roque de los Muchachos / La Palma",
        latitude_deg=28.7606,
        longitude_deg=-17.8792,
        elevation_m=2396.0,
        timezone="Atlantic/Canary",
        horizon_deg=12.0,
        zenith_deg=88.0,
    ),
    "MAUNAKEA": Observatory(
        name="Mauna Kea",
        latitude_deg=19.8206,
        longitude_deg=-155.4681,
        elevation_m=4205.0,
        timezone="Pacific/Honolulu",
        horizon_deg=15.0,
        zenith_deg=88.0,
    ),
    "PARANAL": Observatory(
        name="Paranal",
        latitude_deg=-24.6272,
        longitude_deg=-70.4042,
        elevation_m=2635.0,
        timezone="America/Santiago",
        horizon_deg=12.0,
        zenith_deg=88.0,
    ),
    "LASILLA": Observatory(
        name="La Silla",
        latitude_deg=-29.2567,
        longitude_deg=-70.7346,
        elevation_m=2400.0,
        timezone="America/Santiago",
        horizon_deg=12.0,
        zenith_deg=88.0,
    ),
    "PALOMAR": Observatory(
        name="Palomar Observatory",
        latitude_deg=33.3563,
        longitude_deg=-116.8650,
        elevation_m=1712.0,
        timezone="America/Los_Angeles",
        horizon_deg=15.0,
        zenith_deg=88.0,
    ),
}


def get_observatory(observatory_id: str) -> Observatory:
    key = observatory_id.upper()
    if key not in OBSERVATORIES:
        raise KeyError(f"Unknown observatory id: {observatory_id}")
    return OBSERVATORIES[key]


def list_observatories() -> list[dict[str, object]]:
    return [{"id": key, **asdict(value)} for key, value in OBSERVATORIES.items()]
