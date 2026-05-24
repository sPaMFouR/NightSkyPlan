from datetime import date

from nightskyplan.ephemeris import make_time_grid, solar_lunar_context
from nightskyplan.observatory import DEFAULT_OBSERVATORY


def test_solar_lunar_context_allows_dates_outside_iers_table():
    times, _ = make_time_grid(date(2026, 5, 25), DEFAULT_OBSERVATORY.timezone, cadence_min=120)

    context = solar_lunar_context(times, DEFAULT_OBSERVATORY, twilight_alt_deg=-18)

    assert len(context["sun_alt"]) == len(times)
    assert len(context["moon_alt"]) == len(times)
