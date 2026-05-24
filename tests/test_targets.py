from astropy import units as u
import pandas as pd
import pytest

from nightskyplan.targets import parse_angle, parse_targets


def test_parse_sexagesimal_ra_and_dec():
    ra = parse_angle("12:00:00", is_ra=True)
    dec = parse_angle("-30:00:00", is_ra=False)

    assert ra.to_value(u.deg) == pytest.approx(180.0)
    assert dec.to_value(u.deg) == pytest.approx(-30.0)


def test_parse_targets_defaults_optional_columns():
    targets = parse_targets(pd.DataFrame({"name": ["A"], "ra": ["180"], "dec": ["-30"]}))

    assert targets.loc[0, "name"] == "A"
    assert targets.loc[0, "exposure_min"] == 30.0
    assert targets.loc[0, "priority"] == 1.0
