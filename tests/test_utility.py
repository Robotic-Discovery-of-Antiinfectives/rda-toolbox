# TODO: Write tests for every function in rda_toolbox/utility.py

# There are zero automated checks that mapapply_96_to_384,
# background_normalize_zfactor, Altair plot builders,
# or MIC result formatting behave as expected.
# Adding lightweight tests
# (e.g., fixtures that feed synthetic plate data through background_normalize_zfactor
#  and verify normalization/Z‑factor math, or verifying plateheatmaps leaves input untouched)
# will catch the stability issues above before they regress again.

import numpy as np
import pandas as pd

from rda_toolbox.utility import (
    mapapply_96_to_384
)


def test_mapapply_96_to_384_maps_quadrants():
    """
    Test that mapapply_96_to_384 correctly maps 96-well plate positions in a Z-Scheme
    """
    rows_96 = list("ABCDEFGH")
    cols_96 = list(range(1, 13))
    df = pd.DataFrame(
        {
            "Row_96": [row for _ in range(4) for row in rows_96 for _ in cols_96],
            "Column_96": cols_96 * (len(rows_96) * 4),
            "Quadrant": [
                quadrant
                for quadrant in range(1, 5)
                for _ in range(len(rows_96) * len(cols_96))
            ],
        }
    )

    result = mapapply_96_to_384(df)

    expected_rows_384 = [
        row
        for quadrant in range(1, 5)
        for row in ("ACEGIKMO" if quadrant in (1, 2) else "BDFHJLNP")
        for _ in cols_96
    ]
    expected_cols_384 = sum([], [
        col
        for quadrant in range(1, 5)
        for col in (list(range(1, 25, 2)) if quadrant in (1, 3) else list(range(2, 25, 2))) * len(rows_96)
    ])
    assert result["Row_384"].tolist() == expected_rows_384
    assert result["Col_384"].tolist() == expected_cols_384

