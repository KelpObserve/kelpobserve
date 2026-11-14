import json
import math

import pandas as pd
import pytest

from kelpobserve.monthly.evaluation import INPUTS
from kelpobserve.monthly.field import aggregate_events


def rows(species, covers):
    return pd.DataFrame(
        {
            "source_row": range(2, len(species) + 2),
            "Species": species,
            "Date": ["2026-07-06"] * len(species),
            "Depth": [7] * len(species),
            "Transect": [1] * len(species),
            "Quadrat": [1] * len(species),
            "% cover": covers,
            "Juvenile (< 30cm)": [1] * len(species),
            "Adult (>= 30cm)": [0] * len(species),
        }
    )


def test_plant_rows_do_not_multiply_cover_and_absence_is_explicit():
    schema = json.loads((INPUTS / "field_schema.json").read_text())
    plant = "saccorhiza polyschides"
    result = aggregate_events(rows([plant] * 3, [60.0, None, None]), schema)
    assert len(result) == 1 and result.cover_pct.iloc[0] == 60
    assert result.plant_rows.iloc[0] == 3
    assert aggregate_events(rows(["EMPTY"], [None]), schema).cover_pct.iloc[0] == 0
    assert math.isnan(aggregate_events(rows([plant], [None]), schema).cover_pct.iloc[0])
    with pytest.raises(ValueError, match="both EMPTY"):
        aggregate_events(rows([plant, "EMPTY"], [60.0, None]), schema)
    with pytest.raises(ValueError, match="Repeated cover"):
        aggregate_events(rows([plant, plant], [60.0, 60.0]), schema)
