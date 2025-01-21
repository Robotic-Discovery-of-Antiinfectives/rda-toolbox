#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import pathlib

from .parser import (
        read_platemapping,
        parse_mappingfile,
        parse_readerfiles,
        read_inputfile,
        )
from .utility import (
        mic_assaytransfer_mapping,
        mapapply_96_to_384
    )


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
    df: pd.DataFrame, # mapped inputs
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
    # df = pd.merge(raw_df, input_df, how="outer")
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

    df[substance_id] = df[substance_id].astype(str)
    return df.round(
        {
            "Denoised Optical Density": 2,
            "Relative Optical Density": 2,
            "Z-Factor": 2,
            "Concentration": 2,
        }
    )


# def process(df: pd.DataFrame) -> pd.DataFrame:
#     processed = df.groupby(
#         ["Internal ID", "External ID", "Organism", "Concentration", "Dataset"],
#         as_index=False,
#     ).agg(
#         {
#             "Internal ID": ["first", "size"],
#             "Relative Optical Density": ["mean", "std"],
#         }
#     )
#     processed.columns = [
#         "External ID",
#         "Organism",
#         "Concentration",
#         "Dataset",
#         "Internal ID",
#         "Num Replicates",
#         "Relative Optical Density mean",
#         "Relative Optical Density std",
#     ]
#     return processed


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
        :
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


# def apply_threshold(
#     df: pd.DataFrame,
#     id_column="Internal ID",
#     negative_controls: str = "Bacteria + Medium",
#     blanks: str = "Medium",
#     measurement: str = "Relative Optical Density mean",
#     threshold=None,
# ) -> pd.DataFrame:
#     """
#     Applies provided threshold to processed data.
#     Expects a DataFrame with columns:
#     - External ID
#     - Organism
#     - Concentration
#     - Dataset
#     - Internal ID
#     - Num Replicates
#     - Relative Optical Density mean
#     - Relative Optical Density std

#     Else provide Cutoff via 'threshold' keyword argument.

#     Returns a smaller DataFrame than was given via input.
#     """

#     # Use only substance entries, no controls, no blanks etc.:
#     substance_df = df.loc[
#         (df[id_column] != blanks) & (df[id_column] != negative_controls), :
#     ].copy()
#     # Apply threshold:
#     if threshold: # overwrite possibly provided cutoff via input df
#         substance_df["Cutoff"] = threshold
#     if not threshold:
#         if "Cutoff" not in substance_df:
#             raise KeyError(
#                 "No threshold argument provided and no 'Cutoff' column in Input.xlsx" +
#                 " E.g.: apply_threshold(processed_data, threshold=50)"
#             )
#     # highest conc. needs to be below the threshold
#     # measurement at any conc. below threshold

#     # filter for groups where the measurement at max. conc. is below the given threshold
#     selection = substance_df.groupby(["Internal ID", "Organism"]).filter(
#         lambda grp: grp[grp["Concentration"] == grp["Concentration"].max()][
#             "Relative Optical Density mean"
#         ]
#         < grp["Cutoff"].mean()
#     )
#     mic_dfs = []
#     non_grouping_columns = [
#         "Concentration",
#         "Num Replicates",
#         "Relative Optical Density mean",
#         "Relative Optical Density std",
#     ]
#     grouping_columns = list(
#         filter(lambda x: x not in non_grouping_columns, selection.columns)
#     )
#     for grp_columns, grp in selection.groupby(grouping_columns):
#         mic_df = pd.DataFrame(
#             {key: [value] for key, value in zip(grouping_columns, grp_columns)}
#         )
#         mic_df[f"MIC {threshold}"] = grp.iloc[
#             (
#                 grp.sort_values(by=["Concentration"])[
#                     "Relative Optical Density mean"
#                 ]
#                 < list(grp["Cutoff"].unique())[0]
#             ).argmax()
#         ]["Concentration"]
#         mic_df[f"Relative Optical Density mean"] = grp.iloc[
#             (
#                 grp.sort_values(by=["Concentration"])[
#                     "Relative Optical Density mean"
#                 ]
#                 < list(grp["Cutoff"].unique())[0]
#             ).argmax()
#         ]["Relative Optical Density mean"]
#         mic_dfs.append(mic_df)

#     return pd.concat(mic_dfs)


def mic_process_inputs(
    substances_file: str,
    ast_mapping_file: str,
    acd_mapping_file: str,
    rawfiles_path: str,
):
    substances, organisms, dilutions, controls = read_inputfile(substances_file)
    # substances = pd.read_excel(substances_file, sheet_name="Substances")
    # organisms = pd.read_excel(substances_file, sheet_name="Organisms")
    # dilutions = pd.read_excel(substances_file, sheet_name="Dilutions")
    # controls = pd.read_excel(substances_file, sheet_name="Controls")

    rawdata = parse_readerfiles(rawfiles_path)

    # Split control position:
    # controls["Row_384"] = controls["Position"].apply(lambda x: x[0])
    # controls["Col_384"] = controls["Position"].apply(lambda x: x[1:])

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
        for col_i, conc in enumerate(list(dilutions["Concentration"].unique())):
            # Add concentration:
            subst_row["Concentration"] = conc
            # Add corresponding column:
            subst_row["Col_384"] = int(col_positions_384[init_pos][col_i])
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

    # merge rawdata with input specifications
    df = pd.merge(rawdata, acd_single_concentrations_df, how="outer")

    return df


def mic_results(df, filepath, thresholds=[20, 50]):
    """
    Expects the results from rda.preprocess() function.
    Means measurements between replicates and obtains the MIC values per substance and organism.
    Saves excel files per dataset and sheets per organism with Minimum Inhibitory Concentrations (MICs)
    at the given thresholds.
    """

    df = df[
        (df["Dataset"] != "Negative Control")
        & (df["Dataset"] != "Blank")
    ].dropna(subset=["Concentration"])
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
            ["Concentration", "Relative Optical Density mean", "Z-Factor mean", "Z-Factor std"]
        ].sort_values(by=["Concentration"])
        #print(grp)
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
            if not max_conc_below_threshold:
                mic = None
            else:
                mic = values_below_threshold.iloc[0]["Concentration"]
            record[f"MIC{threshold} in µM"] = mic
        mic_records.append(record)
    # Drop entries where no MIC could be determined
    mic_df = pd.DataFrame.from_records(mic_records)
    # mic_df.dropna(
    #     subset=[f"MIC{threshold} in µM" for threshold in thresholds],
    #     how="all",
    #     inplace=True,
    # )
    mic_df.round(2).to_excel(os.path.join(filepath, "MIC_Results_AllDatasets_longformat.xlsx"), index=False)
    for dataset, dataset_grp in mic_df.groupby(["Dataset"]):
        pivot_multiindex_df = pd.pivot_table(
            dataset_grp,
            values=[f"MIC{threshold} in µM" for threshold in thresholds] + ["Z-Factor mean", "Z-Factor std"],
            index=["Internal ID", "External ID", "Dataset"],
            columns="Organism",
        ).reset_index()

        resultpath = os.path.join(filepath, dataset[0])

        # References special case
        if dataset[0] == "Reference":
            references_mic_results(df, resultpath, thresholds=thresholds)
            continue # skip for references

        pathlib.Path(resultpath).mkdir(parents=True, exist_ok=True)
        for threshold in thresholds:
            organisms_thresholded_mics = pivot_multiindex_df[
                ["Internal ID", "External ID", f"MIC{threshold} in µM"]
            ]
            cols = list(organisms_thresholded_mics.columns.droplevel())
            cols[0] = "Internal ID"
            cols[1] = "External ID"
            organisms_thresholded_mics.columns = cols
            organisms_thresholded_mics = (
                organisms_thresholded_mics.sort_values(
                    by=list(organisms_thresholded_mics.columns)[2:],
                    na_position="last",
                )
            )
            # organisms_thresholded_mics.dropna(
            #     subset=list(organisms_thresholded_mics.columns)[2:],
            #     how="all",
            #     inplace=True,
            # )
            organisms_thresholded_mics.fillna("NA", inplace=True)
            organisms_thresholded_mics.to_excel(
                os.path.join(
                    resultpath, f"{dataset[0]}_MIC{threshold}_results.xlsx"
                ),
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
    only_references = preprocessed_data[
        preprocessed_data["Dataset"] == "Reference"
    ]
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
            values_below_threshold = grp[
                grp["Relative Optical Density"] < threshold
            ]
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
            for _, org_df in pd.merge(
                mapped_organisms, rawdata
            ).groupby("Organism")
        ]
    )

    for ast_barcode, ast_plate in result_df.groupby("AsT Barcode 384"):
        print(f"AsT Plate {ast_barcode} has size: {len(ast_plate)//len(ast_plate['AcD Barcode 384'].unique())}")
        print(f"{ast_barcode} -> {ast_plate["AcD Barcode 384"].unique()}")

    return result_df


def primary_results(
    df: pd.DataFrame,
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
            "Internal ID",
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
            ["Internal ID", "Organism", "Dataset"]
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
                index=["Internal ID", "Dataset", "Concentration"],
                columns="Organism",
            ).reset_index()
            cols = list(pivot_multiindex_df.columns.droplevel())
            cols[:3] = list(
                map(lambda x: x[0], pivot_multiindex_df.columns[:3])
            )
            pivot_multiindex_df.columns = cols

            # Apply threshold (active in any organism)
            thresholded_pivot = pivot_multiindex_df.iloc[
                  list(pivot_multiindex_df.iloc[:, 3:].apply(lambda x: any(list(map(lambda i: i < threshold, x))), axis=1))
            ]

            # Sort by columns each organism after the other
            # return pivot_multiindex_df.sort_values(by=cols[3:])

            # Sort rows by mean between the organisms (lowest mean activity first)
            results_sorted_by_mean_activity = thresholded_pivot.iloc[
                thresholded_pivot.iloc[:, 3:].mean(axis=1).argsort()
            ]
            print("Saving", os.path.join(resultpath, f"{dataset}_threshold{threshold}_results.xlsx"))
            results_sorted_by_mean_activity.to_excel(
                os.path.join(resultpath, f"{dataset}_threshold{threshold}_results.xlsx"),
                index=False,
            )
            print("Saving", os.path.join(resultpath, f"{dataset}_threshold{threshold}_results.csv"))
            results_sorted_by_mean_activity.to_csv(
                os.path.join(resultpath, f"{dataset}_threshold{threshold}_results.csv"),
                index=False,
            )
