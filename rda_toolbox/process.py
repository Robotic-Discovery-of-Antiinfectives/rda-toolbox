#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import pathlib

from scipy.stats import median_abs_deviation

from .parser import (
    read_platemapping,
    parse_mappingfile,
    parse_readerfiles,
    read_inputfile,
)
from .utility import mic_assaytransfer_mapping, mapapply_96_to_384


def zfactor(positive_controls, negative_controls):
    return 1 - (
        3
        * (np.std(positive_controls) + np.std(negative_controls))
        / abs(np.mean(positive_controls - np.mean(negative_controls)))
    )


def median_absolute_deviation(x, scale=1.0):
    med = x.median()
    x = abs(x - med)
    mad = x.median()
    return mad / scale


def zfactor_median(positive_controls, negative_controls):
    return 1 - (
        3
        * (
            median_absolute_deviation(positive_controls)
            + median_absolute_deviation(negative_controls)
        )
        / abs(np.median(positive_controls - np.median(negative_controls)))
    )


def median_polish_df(
    plate_df: pd.DataFrame,
    max_iter: int = 100,
    tol: float = 1e-4,
    measurement_header: str = "Raw Optical Density",
    row_header: str = "Row_384",
    col_header: str = "Col_384",
) -> pd.DataFrame:
    plate_df = plate_df.copy()

    plate_df["row_effect"] = 0
    plate_df["col_effect"] = 0

    overall_effect = plate_df[measurement_header].median()
    plate_df[measurement_header] -= overall_effect

    for iteration in range(max_iter):

        # ------------------------------------------------------------
        plate_df["row_median"] = plate_df.groupby(row_header)[
            measurement_header
        ].transform("median")
        plate_df[measurement_header] -= plate_df["row_median"]
        plate_df["row_effect"] += plate_df["row_median"]
        overall_effect += plate_df["row_effect"].median()
        plate_df["row_effect"] -= plate_df["row_effect"].median()
        # ------------------------------------------------------------
        plate_df["col_median"] = plate_df.groupby(col_header)[
            measurement_header
        ].transform("median")
        plate_df[measurement_header] -= plate_df["col_median"]
        plate_df["col_effect"] += plate_df["col_median"]
        overall_effect += plate_df["col_effect"].median()
        plate_df["col_effect"] -= plate_df["col_effect"].median()
        # ------------------------------------------------------------

        if (
            plate_df["row_median"].abs().max() < tol
            and plate_df["col_median"].abs().max() < tol
        ):
            # print(f"tol break at iteration {iteration}")
            break
    return plate_df


def add_b_score(
    plate_df: pd.DataFrame,
    measurement_header: str = "Raw Optical Density",
    row_header: str = "Row_384",
    col_header: str = "Col_384",
) -> pd.DataFrame:
    """
    Expects a Dataframe comprising a **whole** plate (without controls!).
    """
    # We could also collect iterations of the median polish function and plot the results to show progress of normalization
    plate_df = median_polish_df(plate_df)
    mad_value = median_absolute_deviation(plate_df[measurement_header], scale=1.4826)
    if not np.isfinite(mad_value) or np.isclose(mad_value, 0.0):
        raise ValueError("Median absolute deviation is zero; cannot compute B-scores.")
    plate_df["b_scores"] = plate_df[measurement_header] / mad_value
    return plate_df.drop(
        columns=["row_effect", "col_effect", "row_median", "col_median", measurement_header]
        ).round({"b_scores": 2})


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

    # Work on a copy and ensure the measurement column is numeric so subtraction is supported
    grp = grp.copy()
    raw_col = f"Raw {measurement}"
    if raw_col not in grp:
        raise KeyError(f"Column '{raw_col}' not found in input DataFrame.")
    grp[raw_col] = pd.to_numeric(grp[raw_col], errors="coerce")

    plate_neg_controls = grp[grp[substance_id] == negative_controls][raw_col]
    plate_blank_controls = grp[grp[substance_id] == blanks][raw_col]

    # Check inputs :)
    if len(plate_neg_controls) == 0:
        raise KeyError("Please check if keyword 'negative_controls' is matching with input table.")
    if len(plate_blank_controls) == 0:
        raise KeyError("Please check if keyword 'blanks' is matching with input table.")
    if grp[raw_col].isna().all():
        raise ValueError("Raw measurement column contains no numeric values.")

    plate_blanks_mean = plate_blank_controls.mean()
    if not np.isfinite(plate_blanks_mean):
        raise ValueError("Blank controls contain non-finite values.")
    # Subtract background noise:
    grp[f"Denoised {measurement}"] = grp[f"Raw {measurement}"] - plate_blanks_mean
    plate_denoised_negative_mean = grp[grp[substance_id] == negative_controls][
        f"Denoised {measurement}"
    ].mean()
    plate_denoised_blank_mean = grp[grp[substance_id] == blanks][
        f"Denoised {measurement}"
    ].mean()
    if not np.isfinite(plate_denoised_negative_mean) or np.isclose(
        plate_denoised_negative_mean, 0.0
    ):
        plate_label = grp[norm_by_barcode].iat[0] if norm_by_barcode in grp else "Unknown"
        raise ValueError(
            f"Plate {plate_label} cannot be normalized: negative controls after background subtraction have near-zero mean."
        )
    # Normalize:
    grp[f"Relative {measurement}"] = (
        grp[f"Denoised {measurement}"] / plate_denoised_negative_mean
    ) * 100
    # Z-Factor:
    grp["Z-Factor"] = zfactor(plate_neg_controls, plate_blank_controls)

    # Robust Z-Factor using median instead of mean:
    grp["Robust Z-Factor"] = zfactor_median(plate_neg_controls, plate_blank_controls)

    return grp


def preprocess(
    df: pd.DataFrame,  # mapped inputs
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
    df[substance_id] = df[substance_id].astype(str)
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

    # df[substance_id] = df[substance_id].astype(str)

    # detect and report NA values (defined in input, not in raw data)
    orgs_w_missing_data = df[df[f"Raw {measurement}"].isna()]["Organism formatted"].unique()
    if orgs_w_missing_data.size > 0:
        print(
            f"""Processed data:
      Organisms with missing data, excluded from processed data: {orgs_w_missing_data}.
      If this is not intended, please check the Input.xlsx or if raw data files are complete.
              """
        )
        df = df.dropna(subset=[f"Raw {measurement}"])
    # Report missing
    # Remove missing from "processed" dataframe
    return df.round(
        {
            "Denoised Optical Density": 2,
            "Relative Optical Density": 2,
            "Z-Factor": 2,
            "Robust Z-Factor": 2,
            # "Concentration": 2,
        }
    )


def get_thresholded_subset(
    df: pd.DataFrame,
    id_column="ID",
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
    substance_df = df.loc[
        (df[id_column] != blanks)
        & (df[id_column] != negative_controls)
        & (df["Organism"] != blankplate_organism),
        :,
    ].copy()
    # Apply threshold:
    if threshold:
        substance_df["Cutoff"] = threshold
    else:
        if "mic_cutoff" not in substance_df:
            raise KeyError("No 'mic_cutoff' column in Input.xlsx")
    selection = substance_df[
        substance_df["Relative Optical Density"] < substance_df["Cutoff"]
    ]
    # Apply mean and std in case of replicates:
    result = selection.groupby([id_column, "Organism", "Dataset"], as_index=False).agg(
        {
            "Relative Optical Density": ["mean", "std"],
            id_column: ["first", "count"],
            "Organism": "first",
            "Cutoff": "first",
            "Dataset": "first",
        }
    )
    result.columns = [
        "Relative Optical Density mean",
        "Relative Optical Density std",
        id_column,
        "Replicates",
        "Organism",
        "Cutoff",
        "Dataset",
    ]
    return result


def mic_results(df, filepath, thresholds=[20, 50]):
    """
    Expects the results from rda.preprocess() function.
    Means measurements between replicates and obtains the MIC values per substance and organism.
    Saves excel files per dataset and sheets per organism with Minimum Inhibitory Concentrations (MICs)
    at the given thresholds.
    """

    df = df[(df["Dataset"] != "Negative Control") & (df["Dataset"] != "Blank")].dropna(
        subset=["Concentration"]
    )
    # the above should remove entries where Concentration == NAN

    # Pivot table to get the aggregated values:
    pivot_df = pd.pivot_table(
        df,
        values=["Relative Optical Density", "Replicate", "Z-Factor"],
        index=[
            "Internal ID",
            "External ID",
            "Organism",
            "Concentration",
            "Dataset",
        ],
        aggfunc={
            "Relative Optical Density": ["mean"],
            "Replicate": ["count"],
            "Z-Factor": ["mean", "std"],  # does this make sense? with std its usable.
            # "Z-Factor": ["std"],
        },
    ).reset_index()

    # merge pandas hirarchical column index (wtf is this pandas!?)
    pivot_df.columns = [" ".join(x).strip() for x in pivot_df.columns.ravel()]

    mic_records = []
    for group_names, grp in pivot_df.groupby(
        ["Internal ID", "External ID", "Organism", "Dataset"]
    ):
        internal_id, external_id, organism, dataset = group_names
        # Sort by concentration just to be sure:
        grp = grp[
            [
                "Concentration",
                "Relative Optical Density mean",
                "Z-Factor mean",
                "Z-Factor std",
            ]
        ].sort_values(by=["Concentration"])  # This sorting is very important for the mic determination!
        # Get rows where the OD is below the given threshold:
        record = {
            "Internal ID": internal_id,
            "External ID": external_id,
            "Organism": organism,
            "Dataset": dataset,
            "Z-Factor mean": list(grp["Z-Factor mean"])[0],
            "Z-Factor std": list(grp["Z-Factor std"])[0],
        }

        for threshold in thresholds:
            values_below_threshold = grp[
                grp["Relative Optical Density mean"] < threshold
            ]
            # thx to jonathan - check if the OD at maximum concentration is below threshold (instead of any concentration)
            max_conc_below_threshold = list(
                grp[grp["Concentration"] == max(grp["Concentration"])][
                    "Relative Optical Density mean"
                ]
                < threshold
            )[0]
            if max_conc_below_threshold:
                mic = values_below_threshold.iloc[0]["Concentration"] # the maximum concentration which is below the threshold
            else:
                mic = None
            record[f"MIC{threshold} in µM"] = mic
        mic_records.append(record)
    # Drop entries where no MIC could be determined
    mic_df = pd.DataFrame.from_records(mic_records)
    # mic_df.dropna(
    #     subset=[f"MIC{threshold} in µM" for threshold in thresholds],
    #     how="all",
    #     inplace=True,
    # )
    mic_df.round(2).to_excel(
        os.path.join(filepath, "MIC_Results_AllDatasets_longformat.xlsx"), index=False
    )
    for dataset, dataset_grp in mic_df.groupby(["Dataset"]):
        pivot_multiindex_df = pd.pivot_table(
            dataset_grp,
            values=[f"MIC{threshold} in µM" for threshold in thresholds]
            + ["Z-Factor mean", "Z-Factor std"],
            index=["Internal ID", "External ID", "Dataset"],
            columns="Organism",
        ).reset_index()

        resultpath = os.path.join(filepath, dataset[0])

        # References special case
        if dataset[0] == "Reference":
            references_mic_results(df, resultpath, thresholds=thresholds)
            continue  # skip for references

        pathlib.Path(resultpath).mkdir(parents=True, exist_ok=True)
        for threshold in thresholds:
            organisms_thresholded_mics = pivot_multiindex_df[
                ["Internal ID", "External ID", f"MIC{threshold} in µM"]
            ]
            cols = list(organisms_thresholded_mics.columns.droplevel())
            cols[0] = "Internal ID"
            cols[1] = "External ID"
            organisms_thresholded_mics.columns = cols
            organisms_thresholded_mics = organisms_thresholded_mics.sort_values(
                by=list(organisms_thresholded_mics.columns)[2:],
                na_position="last",
            )
            # organisms_thresholded_mics.dropna(
            #     subset=list(organisms_thresholded_mics.columns)[2:],
            #     how="all",
            #     inplace=True,
            # )
            organisms_thresholded_mics.fillna("NA", inplace=True)
            organisms_thresholded_mics.to_excel(
                os.path.join(resultpath, f"{dataset[0]}_MIC{threshold}_results.xlsx"),
                index=False,
            )


def references_mic_results(
    preprocessed_data,
    resultpath,
    thresholds=[20, 50],
):
    """
    This function saves an excel file for the reference substances.
    Since reference substances have duplicate Internal IDs
    (since they are used multiple times), they would be meaned between duplicates.
    To circumvent this, this function exists which gets the MIC
    for each reference **per (AcD) plate** instead of per Internal ID.
    """
    only_references = preprocessed_data[preprocessed_data["Dataset"] == "Reference"]
    mic_records = []
    for group_names, grp in only_references.groupby(
        [
            "Internal ID",
            "External ID",
            "Organism",
            "Dataset",
            "AcD Barcode 384",
        ]
    ):
        internal_id, external_id, organism, dataset, acd_barcode = group_names
        grp = grp.copy().sort_values(by=["Concentration"])
        record = {
            "Internal ID": internal_id,
            "External ID": external_id,
            "Organism": organism,
            "Dataset": dataset,
            "AcD Barcode 384": acd_barcode,
            "Z-Factor": list(grp["Z-Factor"])[0],
        }
        for threshold in thresholds:
            values_below_threshold = grp[grp["Relative Optical Density"] < threshold]
            # thx to jonathan - check if the OD at maximum concentration is below threshold (instead of any concentration)
            max_conc_below_threshold = list(
                grp[grp["Concentration"] == max(grp["Concentration"])][
                    "Relative Optical Density"
                ]
                < threshold
            )[0]
            if not max_conc_below_threshold:
                mic = None
            else:
                mic = values_below_threshold.iloc[0]["Concentration"]
            record[f"MIC{threshold} in µM"] = mic
        mic_records.append(record)
    mic_df = pd.DataFrame.from_records(mic_records)
    mic_df.sort_values(by=["External ID", "Organism"]).to_excel(
        os.path.join(resultpath, "References_MIC_results_eachRefID.xlsx"), index=False
    )


def primary_process_inputs(
    inputfile_path,
    mappingfile_path,
    rawfiles_path,
    map_rowname="Row 96",
    map_colname="Col 96",
    q_name="Quadrant",
):
    substances, organisms, dilutions, controls = read_inputfile(inputfile_path)
    rawdata = parse_readerfiles(rawfiles_path)
    mapapply_96_to_384(
        substances, rowname=map_rowname, colname=map_colname, q_name=q_name
    )

    mapping_df = parse_mappingfile(
        mappingfile_path,
        # motherplate_column="Origin Barcode",
        motherplate_column="AsT Barcode 384",
        childplate_column="AcD Barcode 384",
    )

    control_wbarcodes = []
    # multiply controls with number of AsT plates to later merge them with substances df
    for origin_barcode in list(substances["AsT Barcode 384"].unique()):
        controls_subdf = controls.copy()
        controls_subdf["AsT Barcode 384"] = origin_barcode
        control_wbarcodes.append(controls_subdf)
    controls_n_barcodes = pd.concat(control_wbarcodes)

    ast_plate_df = pd.merge(
        pd.concat([substances, controls_n_barcodes]), dilutions, how="outer"
    )

    mapped_organisms = pd.merge(mapping_df, organisms)

    result_df = pd.concat(
        [
            pd.merge(org_df, ast_plate_df)
            for _, org_df in pd.merge(mapped_organisms, rawdata).groupby("Organism")
        ]
    )

    for ast_barcode, ast_plate in result_df.groupby("AsT Barcode 384"):
        print(
            f"AsT Plate {ast_barcode} has size: {len(ast_plate)//len(ast_plate['AcD Barcode 384'].unique())}"
        )
        print(f"{ast_barcode} -> {ast_plate["AcD Barcode 384"].unique()}")

    return result_df


def primary_results(
    df: pd.DataFrame,
    substance_id,
    filepath="../data/results/",
    thresholds: list[float] = [50],
):
    """
    Expects the results from rda.preprocess() function.
    """
    df = df[
        (df["Dataset"] != "Reference")
        & (df["Dataset"] != "Positive Control")
        & (df["Dataset"] != "Blank")
    ].dropna(subset=["Concentration"])

    pivot_df = pd.pivot_table(
        df,
        values=["Relative Optical Density", "Replicate", "Z-Factor"],
        index=[
            substance_id,
            "Organism",
            "Concentration",
            "Dataset",
        ],
        aggfunc={
            "Relative Optical Density": ["mean"],
            "Replicate": ["count"],
        },
    ).reset_index()
    pivot_df.columns = [" ".join(x).strip() for x in pivot_df.columns.ravel()]

    for threshold in thresholds:
        pivot_df[f"Relative Growth < {threshold}"] = pivot_df.groupby(
            [substance_id, "Organism", "Dataset"]
        )["Relative Optical Density mean"].transform(lambda x: x < threshold)

        for dataset, dataset_grp in pivot_df.groupby(["Dataset"]):
            dataset = dataset[0]
            resultpath = os.path.join(filepath, dataset)
            pathlib.Path(resultpath).mkdir(parents=True, exist_ok=True)

            print(
                "Saving",
                os.path.join(resultpath, f"{dataset}_all_results.xlsx"),
            )
            dataset_grp.to_excel(
                os.path.join(resultpath, f"{dataset}_all_results.xlsx"),
                index=False,
            )
            print(
                "Saving",
                os.path.join(resultpath, f"{dataset}_all_results.csv"),
            )
            dataset_grp.to_csv(
                os.path.join(resultpath, f"{dataset}_all_results.csv"),
                index=False,
            )

            pivot_multiindex_df = pd.pivot_table(
                dataset_grp,
                values=[f"Relative Optical Density mean"],
                index=[substance_id, "Dataset", "Concentration"],
                columns="Organism",
            ).reset_index()
            cols = list(pivot_multiindex_df.columns.droplevel())
            cols[:3] = list(map(lambda x: x[0], pivot_multiindex_df.columns[:3]))
            pivot_multiindex_df.columns = cols

            # Apply threshold (active in any organism)
            thresholded_pivot = pivot_multiindex_df.iloc[
                list(
                    pivot_multiindex_df.iloc[:, 3:].apply(
                        lambda x: any(list(map(lambda i: i < threshold, x))), axis=1
                    )
                )
            ]

            # Sort by columns each organism after the other
            # return pivot_multiindex_df.sort_values(by=cols[3:])

            # Sort rows by mean between the organisms (lowest mean activity first)
            results_sorted_by_mean_activity = thresholded_pivot.iloc[
                thresholded_pivot.iloc[:, 3:].mean(axis=1).argsort()
            ]
            print(
                "Saving",
                os.path.join(
                    resultpath, f"{dataset}_threshold{threshold}_results.xlsx"
                ),
            )
            results_sorted_by_mean_activity.to_excel(
                os.path.join(
                    resultpath, f"{dataset}_threshold{threshold}_results.xlsx"
                ),
                index=False,
            )
            print(
                "Saving",
                os.path.join(resultpath, f"{dataset}_threshold{threshold}_results.csv"),
            )
            results_sorted_by_mean_activity.to_csv(
                os.path.join(resultpath, f"{dataset}_threshold{threshold}_results.csv"),
                index=False,
            )
