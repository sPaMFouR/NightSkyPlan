from __future__ import annotations

from pathlib import Path
from typing import IO

import pandas as pd

from nightskyplan.targets import parse_targets


def read_target_csv(source: str | Path | IO[str]) -> pd.DataFrame:
    return parse_targets(pd.read_csv(source))
