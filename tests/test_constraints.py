import pandas as pd

from nightskyplan.constraints import visibility_mask


def test_visibility_mask_rejects_failed_constraints():
    tracks = pd.DataFrame(
        {
            "alt_deg": [40.0, 20.0, 50.0, 60.0],
            "airmass": [1.4, 1.2, 3.1, 1.1],
            "moon_sep_deg": [60.0, 60.0, 60.0, 20.0],
            "ha_hour": [1.0, 1.0, 1.0, 1.0],
        }
    )

    mask = visibility_mask(
        tracks,
        horizon_deg=25.0,
        zenith_deg=85.0,
        max_airmass=2.5,
        min_moon_sep_deg=30.0,
        ha_limit_hour=6.0,
    )

    assert mask.tolist() == [True, False, False, False]
