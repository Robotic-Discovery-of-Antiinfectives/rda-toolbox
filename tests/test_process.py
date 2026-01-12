import numpy as np
import pandas as pd

from rda_toolbox.process import (
    background_normalize_zfactor
)


def test_background_normalize_zfactor_computes_expected_values():
    measurement = "Optical Density"
    negative_controls = "Negative Control"
    blanks = "Blank"
    substance_id = "ID"
    norm_by_barcode = "Barcode"

    df = pd.DataFrame(
        {
            substance_id: [
                blanks,
                blanks,
                negative_controls,
                negative_controls,
                negative_controls,
                "Sample A",
                "Sample B",
            ],
            f"Raw {measurement}": [0.1, 0.2, 0.8, 0.9, 1.0, 0.5, 1.2],
            norm_by_barcode: ["Plate-1"] * 7,
        }
    )

    result = background_normalize_zfactor(
        df.copy(),
        substance_id,
        measurement,
        negative_controls,
        blanks,
        norm_by_barcode,
    )

    blanks_mean = np.mean([0.1, 0.2])
    denoised = df[f"Raw {measurement}"] - blanks_mean
    neg_mean = denoised[df[substance_id] == negative_controls].mean()
    expected_relative = (denoised / neg_mean) * 100
    expected_zfactor = 1 - (
        3 * (np.std([0.8, 0.9, 1.0]) + np.std([0.1, 0.2])) / abs(np.mean([0.8, 0.9, 1.0]) - blanks_mean)
    )
    expected_robust_zfactor = 1 - (
        3
        * (
            np.median(abs(np.array([0.8, 0.9, 1.0]) - np.median([0.8, 0.9, 1.0])))
            + np.median(abs(np.array([0.1, 0.2]) - np.median([0.1, 0.2])))
        )
        / abs(np.median([0.8, 0.9, 1.0]) - np.median([0.1, 0.2]))
    )

    pd.testing.assert_series_equal(
        result[f"Denoised {measurement}"].reset_index(drop=True),
        denoised.reset_index(drop=True),
        check_names=False,
    )
    pd.testing.assert_series_equal(
        result[f"Relative {measurement}"].reset_index(drop=True),
        expected_relative.reset_index(drop=True),
        check_names=False,
    )
    assert np.allclose(result["Z-Factor"].iloc[0], expected_zfactor)
    assert np.allclose(result["Robust Z-Factor"].iloc[0], expected_robust_zfactor)
