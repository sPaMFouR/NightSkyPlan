from datetime import datetime

import pandas as pd
from astropy.time import Time

from nightskyplan.observatory import DEFAULT_OBSERVATORY
from nightskyplan.plotting import visibility_figure


def test_visibility_figure_accepts_datetime_twilight_annotations():
    tracks = pd.DataFrame(
        {
            "name": ["A", "A"],
            "time": [datetime(2026, 1, 1, 18), datetime(2026, 1, 1, 19)],
            "alt_deg": [35.0, 45.0],
            "airmass": [1.7, 1.4],
            "ha_hour": [-1.0, 0.0],
            "moon_sep_deg": [80.0, 85.0],
            "valid": [True, True],
        }
    )
    context = {
        "dusk": Time(datetime(2026, 1, 1, 18)),
        "dawn": Time(datetime(2026, 1, 2, 6)),
    }

    fig = visibility_figure(tracks, context, ["A"], DEFAULT_OBSERVATORY)

    assert len(fig.layout.annotations) >= 2
