#!/usr/bin/env python3

import numpy as np
import pandas as pd

from .parser import read_platemapping
from .utility import mic_assaytransfer_mapping


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


def mic_process_inputs(
    substances_file: str,
    ast_mapping_file: str,
    acd_mapping_file: str,
):
    substances = pd.read_excel(substances_file, sheet_name="Substances")
    organisms = pd.read_excel(substances_file, sheet_name="Organisms")
    dilutions = pd.read_excel(substances_file, sheet_name="Dilutions")
    controls = pd.read_excel(substances_file, sheet_name="Controls")

    # Split control position:
    controls["Row_384"] = controls["Position"].apply(lambda x: x[0])
    controls["Col_384"] = controls["Position"].apply(lambda x: x[1:])

    organisms = list(organisms["Organism"])

    # input_df = pd.read_excel(substances_file)
    ast_platemapping, _ = read_platemapping(
        ast_mapping_file, substances["MP Barcode 96"].unique()
    )

    # Do some sanity checks:
    necessary_columns = [
        "Dataset",
        "Internal ID",
        "MP Barcode 96",
        "MP Position 96",
    ]
    # Check if all necessary column are present in the input table:
    if not all(column in substances.columns for column in necessary_columns):
        raise ValueError(
            f"Not all necessary columns are present in the input table.\n(Necessary columns: {necessary_columns})"
        )
    # Check if all of the necessary column are complete:
    if substances[necessary_columns].isnull().values.any():
        raise ValueError(
            "Input table incomplete, contains NA (missing) values."
        )
    # Check if there are duplicates in the internal IDs (apart from references)
    if any(substances[substances["Dataset"] != "Reference"]["Internal ID"].duplicated()):
        raise ValueError("Duplicate Internal IDs.")

    # Map AssayTransfer barcodes to the motherplate barcodes:
    substances["Row_384"], substances["Col_384"], substances["AsT Barcode 384"] = (
        zip(
            *substances.apply(
                lambda row: mic_assaytransfer_mapping(
                    row["MP Position 96"],
                    row["MP Barcode 96"],
                    ast_platemapping,
                ),
                axis=1,
            )
        )
    )
    acd_platemapping, replicates_dict = read_platemapping(
        acd_mapping_file, substances["AsT Barcode 384"].unique()
    )

    num_replicates = list(set(replicates_dict.values()))[0]
    print(f"""
Rows expected without concentrations:\n
{len(substances["Internal ID"].unique())} (unique substances) * {len(organisms)} (organisms) * {num_replicates} (replicates) = {len(substances["Internal ID"].unique()) * 5 * 3}
    """)
    print(f"""
Rows expected with concentrations:\n
{len(substances["Internal ID"].unique())} (unique substances) * {len(organisms)} (organisms) * {num_replicates} (replicates) * (11 (concentrations) + 1 (Medium/Blank or Negative Control)) = {len(substances["Internal ID"].unique()) * len(organisms) * num_replicates * (11 + 1) }
    """)
    single_subst_concentrations = []

    for substance, subst_row in substances.groupby("Internal ID"):
        # Collect the concentrations each as rows for a single substance:
        single_subst_conc_rows = []
        init_pos = int(subst_row["Col_384"].iloc[0]) - 1
        col_positions_384 = [list(range(1, 23, 2)), list(range(2, 23, 2))]
        for col_i, conc in enumerate(list(dilutions["Concentration"])):
            # Add concentration:
            subst_row["Concentration"] = conc
            # Add corresponding column:
            subst_row["Col_384"] = str(col_positions_384[init_pos][col_i])
            single_subst_conc_rows.append(subst_row.copy())

        # Concatenate all concentrations rows for a substance in a dataframe
        single_subst_concentrations.append(pd.concat(single_subst_conc_rows))
    # Concatenate all substances dataframes to one whole
    input_w_concentrations = pd.concat(single_subst_concentrations)

    acd_dfs_list = []
    for ast_barcode, ast_plate in input_w_concentrations.groupby("AsT Barcode 384"):
        controls["AsT Barcode 384"] = list(ast_plate["AsT Barcode 384"].unique())[0]
        ast_plate = pd.concat([ast_plate, controls])
        for org_i, organism in enumerate(organisms):
            for replicate in range(num_replicates):
                # Add the AcD barcode
                ast_plate["AcD Barcode 384"] = acd_platemapping[ast_barcode][
                    replicate
                ][org_i]

                ast_plate["Replicate"] = replicate + 1
                # Add the scientific Organism name
                ast_plate["Organism"] = organism
                acd_dfs_list.append(ast_plate.copy())
                # Add concentrations:
    acd_single_concentrations_df = pd.concat(acd_dfs_list)
    return acd_single_concentrations_df
