from __future__ import annotations

from nightskyplan.models import Observatory


DEFAULT_OBSERVATORY = Observatory(
    name="HCT / Hanle",
    latitude_deg=32.7794,
    longitude_deg=78.9642,
    elevation_m=4486.0,
    timezone="Asia/Kolkata",
    horizon_deg=25.0,
    zenith_deg=85.0,
)
