#!/usr/bin/env python3

import numpy as np
import pandas as pd


def zfactor(positive_controls, negative_controls):
    return 1 - (
        3
        * (np.std(positive_controls) + np.std(negative_controls))
        / abs(np.mean(positive_controls - np.mean(negative_controls)))
    )


def minmax_normalization(x, minimum, maximum):
    return ((x - minimum) / (maximum - minimum)) * 100


def max_normalization(x, maximum):
    return (x / maximum) * 100


def background_normalize_zfactor(
    grp: pd.DataFrame,
    substance_id,
    measurement,
    negative_controls,
    blanks,
    norm_by_barcode,
) -> pd.DataFrame:
    """
    This function is supposed to be applied to a grouped DataFrame.
    It does the following operations:
    - Background subtraction by subtracting the mean of the blanks per plate
    - Normalization by applying max-normalization using the 'Negative Controls'
    - Z-Factor calculation using negative controls and blanks

    *`negative_controls` are controls with organism (e.g. bacteria) and medium*
    *and are labeled in the input DataFrame as 'Negative Controls'.*
    *`blanks` are controls with only medium and are labeled*
    *in the input DataFrame as 'Medium'.*
    """
    plate_blanks_mean = grp[grp[substance_id] == blanks][f"Raw {measurement}"].mean()
    # Subtract background noise:
    grp[f"Denoised {measurement}"] = grp[f"Raw {measurement}"] - plate_blanks_mean
    plate_denoised_negative_mean = grp[grp[substance_id] == negative_controls][
        f"Denoised {measurement}"
    ].mean()
    plate_denoised_blank_mean = grp[grp[substance_id] == blanks][
        f"Denoised {measurement}"
    ].mean()
    # Normalize:
    grp[f"Relative {measurement}"] = grp[f"Denoised {measurement}"].apply(
        lambda x: max_normalization(x, plate_denoised_negative_mean)
    )
    # Z-Factor:
    plate_neg_controls = grp[grp[substance_id] == negative_controls][
        f"Raw {measurement}"
    ]
    plate_blank_controls = grp[grp[substance_id] == blanks][f"Raw {measurement}"]
    grp["Z-Factor"] = zfactor(plate_neg_controls, plate_blank_controls)

    return grp


def preprocess(
    raw_df: pd.DataFrame,
    input_df: pd.DataFrame,
    substance_id: str = "ID",
    measurement: str = "Optical Density",
    negative_controls: str = "Negative Control",
    blanks: str = "Blank",
    norm_by_barcode="Barcode",
) -> pd.DataFrame:
    """
    - raw_df: raw reader data obtained with `rda.readerfiles_rawdf()`
    - input_df: input specifications table with required columns:
        - Dataset (with specified references as their own dataset 'Reference')
        - ID (substance_id) (with specified blanks and negative_controls)
        - Assay Transfer Barcode
        - Row_384 (or Row_96)
        - Col_384 (or Col_96)
        - Concentration
        - Replicate (specifying replicate number)
        - Organism (scientific organism name i.e. with strain)
    ---
    Processing function which merges raw reader data (raw_df)
    with input specifications table (input_df) and then
    normalizes, calculates Z-Factor per plate (norm_by_barcode)
    and rounds to sensible decimal places.
    """
    # merging reader data and input specifications table
    df = pd.merge(raw_df, input_df, how="outer")
    df = (
        df.groupby(norm_by_barcode)[df.columns]
        .apply(
            lambda grp: background_normalize_zfactor(
                grp,
                substance_id,
                measurement,
                negative_controls,
                blanks,
                norm_by_barcode,
            )
        )
        .reset_index(drop=True)
    )
    return df.round(
        {
            "Denoised Optical Density": 2,
            "Relative Optical Density": 2,
            "Z-Factor": 2,
            "Concentration": 2,
        }
    )


def get_thresholded_subset(
    df: pd.DataFrame,
    negative_controls: str = "Negative Control",
    blanks: str = "Medium",
    blankplate_organism: str = "Blank",
    threshold=None,
) -> pd.DataFrame:
    """
    Expects a DataFrame with a mic_cutoff column
    """
    # TODO: hardcode less columns

    # Use only substance entries, no controls, no blanks etc.:
    substance_df = df[
        (df["ID"] != blanks)
        & (df["ID"] != negative_controls)
        & (df["Organism"] != blankplate_organism)
    ]
    # Apply threshold:
    if threshold:
        substance_df["Cutoff"] = threshold
    else:
        if "mic_cutoff" not in substance_df:
            raise KeyError("Noo 'mic_cutoff' column in Input.xlsx")
    selection = substance_df[
        substance_df["Relative Optical Density"] < substance_df["Cutoff"]
    ]
    # Apply mean and std in case of replicates:
    result = selection.groupby(["ID", "Organism"], as_index=False).agg(
        {
            "Relative Optical Density": ["mean", "std"],
            "ID": ["first", "count"],
            "Organism": "first",
            "Cutoff": "first",
        }
    )
    result.columns = [
        "Relative Optical Density mean",
        "Relative Optical Density std",
        "ID",
        "Replicates",
        "Organism",
        "Cutoff",
    ]
    return result
