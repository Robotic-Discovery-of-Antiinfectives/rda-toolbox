#!/usr/bin/env python3

import numpy as np
import pandas as pd
import string

import base64
import io
import rdkit
from rdkit.Chem import Draw
from rdkit import Chem
import math
import pathlib
import os


def get_rows_cols(platetype: int) -> tuple[int, int]:
    """
    Obtain number of rows and columns as tuple for corresponding plate type.
    """
    match platetype:
        case 96:
            return 8, 12
        case 384:
            return 16, 24
        case _:
            raise ValueError("Not a valid plate type")


def generate_inputtable(readout_df=None, platetype: int = 384):
    """
    Generates an input table for the corresponding readout dataframe.
    If not readout df is provided, create a minimal input df.
    """
    if readout_df is None:
        barcodes = ["001PrS01001"]
    else:
        barcodes = readout_df["Barcode"].unique()

    substance_df = pd.DataFrame(
        {
            "ID": [f"Substance {i}" for i in range(1, platetype + 1)],
            f"Row_{platetype}": [*list(string.ascii_uppercase[:16]) * 24],
            f"Col_{platetype}": sum([[i] * 16 for i in range(1, 24 + 1)], []),
            "Concentration in mg/mL": 1,
        }
    )
    layout_df = pd.DataFrame(
        {
            "Barcode": barcodes,
            "Replicate": [1] * len(barcodes),
            "Organism": [
                f"Placeholder Organism {letter}"
                for letter in string.ascii_uppercase[: len(barcodes)]
            ],
        }
    )
    df = pd.merge(layout_df, substance_df, how="cross")
    return df


def map_96_to_384(
    df_row: pd.Series,
    rowname: str,
    colname: str,
    q_name: str,
) -> tuple[pd.Series, pd.Series]:
    """
    Maps the rows and columns of 4 96-Well plates into a single 384-Well plate.

    - Maps in Z order.
    - Takes row, column and quadrant (each of the 96-well plates is one quadrant) of a well from 4 96-well plates and maps it to the corresponding well in a 384-well plate

    Returns the 384-Well plate row and column.
    Example: `df["Row_384"], df["Col_384"] = zip(*df.apply(map_96_to_384, axis=1))`
    """
    # TODO: Write tests for this mapping function

    row = df_row[rowname]  # 96-well plate row
    col = df_row[colname]  # 96-well plate column
    quadrant = df_row[q_name]  # which of the 4 96-well plate

    rowmapping = dict(
        zip(
            string.ascii_uppercase[0:8],
            np.array_split(list(string.ascii_uppercase)[0:16], 8),
        )
    )
    colmapping = dict(zip(list(range(1, 13)), np.array_split(list(range(1, 25)), 12)))
    row_384 = rowmapping[row][0 if quadrant in [1, 2] else 1]
    col_384 = colmapping[col][0 if quadrant in [1, 3] else 1]
    return row_384, col_384


def mapapply_96_to_384(
    df: pd.DataFrame,
    rowname: str = "Row_96",
    colname: str = "Column_96",
    q_name: str = "Quadrant",
) -> pd.DataFrame:
    """Apply to a DataFrame the mapping of 96-well positions to 384-well positions.

    - Maps in Z order.
    The DataFrame has to have columns with:
    - 96-well plate row positions
    - 96-well plate column positions
    - 96-well plate to 384-well plate quadrants
    *(4 96-well plates fit into 1 384-well plate)*
    """
    df["Row_384"], df["Col_384"] = zip(
        *df.apply(
            lambda row: map_96_to_384(
                row, rowname=rowname, colname=colname, q_name=q_name
            ),
            axis=1,
        )
    )
    return df


def position_to_rowcol(pos: str) -> tuple[str, int]:
    """
    Splits a position like "A1" into row and col e.g. ("A", "1").
    """
    return pos[0], int(pos[1:])


def split_position(
    df: pd.DataFrame,
    position: str = "Position",
    row: str = "Row_384",
    col: str = "Col_384",
) -> pd.DataFrame:
    """
    Split a position like "A1" into row and column positions ("A", 1) and adds them as columns to the DataFrame.
    Hint: Remove NAs before applying this function. E.g. `split_position(df.dropna(subset="Position"))`
    """
    df[row] = df[position].apply(lambda x: str(x[0]))
    df[col] = df[position].apply(lambda x: str(x[1:]))
    return df


def get_selection(df, threshold_value, x_column="Relative Optical Density"):
    """
    Apply this ahead of get_upsetplot_df (to obtain dummies df).
    After all the above, apply UpSetAltair.
    """
    selection_results = df[df[x_column] < threshold_value].copy()
    selection_results["threshold"] = f"<{threshold_value}"

    return selection_results


def get_upsetplot_df(df, set_column="Organism", counts_column="ID"):
    """
    Function to obtain a correctly formatted DataFrame.
    According to [UpSetR-shiny](https://github.com/hms-dbmi/UpSetR-shiny)
    this table is supposed to be encoded in binary and set up so that each column represents a set, and each row represents an element.
    If an element is in the set it is represented as a 1 in that position. If an element is not in the set it is represented as a 0.

    *Thanks to: https://stackoverflow.com/questions/37381862/get-dummies-for-pandas-column-containing-list*
    """
    tmp_df = (
        df.groupby(counts_column)[set_column].apply(lambda x: x.unique()).reset_index()
    )
    dummies_df = (
        pd.get_dummies(
            tmp_df.join(
                pd.Series(
                    tmp_df[set_column]
                    .apply(pd.Series)
                    .stack()
                    .reset_index(1, drop=True),
                    name=set_column + "1",
                )
            )
            .drop(set_column, axis=1)
            .rename(columns={set_column + "1": set_column}),
            columns=[set_column],
        )
        .groupby(counts_column, as_index=False)
        .sum()
    )
    # remove "{set_column}_" from set column labels
    dummies_df.columns = list(
        map(
            lambda x: "".join(x.split("_")[1:]) if x.startswith(set_column) else x,
            dummies_df.columns,
        )
    )
    # remove any dots as they interfere with altairs plotting.
    dummies_df.columns = dummies_df.columns.str.replace(".", "")
    return dummies_df.drop(columns=["count"], errors="ignore")


def chunks(l, n):
    """
    Useful function if you want to put a certain amount
    of observations into one plot.
    Yield n number of striped chunks from l.
    """
    for i in range(0, n):
        yield l[i::n]


def mic_assaytransfer_mapping(position, orig_barcode, ast_platemapping):
    """
    This is a rather unfinished function to map 96 well motherplates to 384 well assay transfer (AsT) plates.

    """
    row = position[0]
    col = int(position[1:])
    orig_barcode = str(orig_barcode)
    rowmapping = dict(
        zip(
            string.ascii_uppercase[0:8],
            np.array_split(list(string.ascii_uppercase)[0:16], 8),
        )
    )
    colmapping = dict(zip(list(range(1, 13)), [1, 2] * 13))
    mapping = {
        1: 0,
        2: 0,
        3: 1,
        4: 1,
        5: 0,
        6: 0,
        7: 1,
        8: 1,
        9: 0,
        10: 0,
        11: 1,
        12: 1,
    }

    row_384 = rowmapping[row][mapping[col]]
    col_384 = colmapping[col]

    # TODO: sometimes only the middle (5-8) or the last part of a motherplate is taken...
    # Currently this would lead to an index error since A5 -> ast_of_3 = 1 but if this is the only AsT plate, the first AsT plate is "missing"...
    # Possible solution would be to try - except decreasing the position by 4 iteratively...
    # try A5 except IndexError: A5 - 4 -> try A1 success etc.
    # or A12 -> ast_of_3 = 2 -> IndexError. A12 - 4 -> A8 -> IndexError. A8 - 4 -> A4 works...
    # seems like a bad solution though.
    if col in [1, 2, 3, 4]:
        ast_of_3 = 0
    elif col in [5, 6, 7, 8]:
        ast_of_3 = 1
    else:
        ast_of_3 = 2
    barcode_384_ast = ast_platemapping[orig_barcode][0][ast_of_3]
    return str(row_384), str(col_384), barcode_384_ast


def mol_to_bytes(mol, format="png"):
    img = Draw.MolToImage(mol)
    buffer = io.BytesIO()
    img.save(buffer, format=format)
    return buffer


def imgbuffer_to_imgstr(imgbuffer, prefix="data:image/png;base64,", suffix=""):
    """
    Encode imagebuffer to string (default base64-encoded string).
    Example: imgbuffer_to_imgstr(mol_to_bytes(mol)), prefix="<img src='data:image/png;base64,", suffix="'/>'")
    """
    str_equivalent_image = base64.b64encode(imgbuffer.getvalue()).decode()
    img_tag = prefix + str_equivalent_image + suffix
    return img_tag


def smiles_to_imgstr(smiles):
    """
    Converts a smiles string to a base64 encoded image string (e.g. for plotting in altair).
    It's a convenience function consisting of rda.utility.mol_to_bytes() and rda.utility.imgbuffer_to_imgstr(),
    use these if you want more fine grained control over the format of the returned string.
    Example: df["image"] = df["smiles"].apply(lambda smiles: smiles_to_imgstr(smiles))
    """
    return imgbuffer_to_imgstr(mol_to_bytes(Chem.MolFromSmiles(smiles)))


def prepare_visualization(df, by_id="Internal ID", whisker_width=1, exclude_negative_zfactors=True):
    """
    Does formatting for the facet lineplots.
    """
    df = df.copy()
    if exclude_negative_zfactors:
        df = df[df["Z-Factor"] > 0]
    df.loc[:, "Used Replicates"] = df.groupby([by_id, "Concentration", "Organism"])[
        ["Replicate"]
    ].transform("count")
    df.loc[:, "Mean Relative Optical Density"] = (
        df.groupby([by_id, "Concentration", "Organism"])[["Relative Optical Density"]]
        .transform("mean")
        .round(2)
    )
    df.loc[:, "Std. Relative Optical Density"] = (
        df.groupby([by_id, "Concentration", "Organism"])[["Relative Optical Density"]]
        .transform("std")
        .round(2)
    )
    df.loc[:, "uerror"] = (
        df["Mean Relative Optical Density"] + df["Std. Relative Optical Density"]
    )
    df.loc[:, "lerror"] = (
        df["Mean Relative Optical Density"] - df["Std. Relative Optical Density"]
    )

    tmp_list = []
    for _, grp in df.groupby([by_id, "Organism"]):
        # use replicate == 1 as the meaned OD is the same in all 3 replicates anyways
        maxconc_below_threshold = (
            grp[(grp["Replicate"] == 1) & (grp["Concentration"] == 50)][
                "Mean Relative Optical Density"
            ]
            < 50
        )
        # print(list(maxconc_below_threshold)[0])
        grp["max_conc_below_threshold"] = list(maxconc_below_threshold)[0]
        tmp_list.append(grp)
        # .sort_values(by=["Concentration"], ascending=False)
        # grp_sorted[grp_sorted["Concentration"] == 50]["Mean Relative Optical Density"]
        # print(grp.aggregate())
    df = pd.concat(tmp_list)
    # print(df)
    # df["highest_conc_bigger_50"] = df.groupby([by_id, "Organism"])[
    #     ["Mean Relative Optical Density"]
    # ].transform(
    #     lambda meas_per_conc: list(meas_per_conc)[0] > 50
    # )
    # print(df)

    df["at_all_conc_bigger_50"] = df.groupby([by_id, "Organism"])[
        ["Mean Relative Optical Density"]
    ].transform(lambda meas_per_conc: all([x > 50 for x in list(meas_per_conc)]))
    # Bin observations into artificial categories for plotting later:
    plot_groups = pd.DataFrame()
    for _, grp in df.groupby(["AsT Barcode 384"]):
        # divide the observations per plate into chunks
        # number of chunks is defined by using a maximum of 10 colors/observations per plot
        num_chunks = math.ceil(len(grp[by_id].unique()) / 10)
        for nr, chunk in enumerate(list(chunks(grp[by_id].unique(), num_chunks))):
            plot_groups = pd.concat(
                [
                    plot_groups,
                    pd.DataFrame(
                        {
                            by_id: chunk,
                            "AsT Plate Subgroup": sum([[nr] * len(chunk)], []),
                        }
                    ),
                ]
            ).reset_index(drop=True)
    df = pd.merge(df, plot_groups)
    return df


def save_plot_per_dataset(
    data: pd.DataFrame,
    plotfunc,
    location: str,
    plotname: str = None,
    saveformats: list[str] = ["svg", "html"],
) -> None:
    """
    This is a convenience function which splits a dataframe into each dataset given in the 'Dataset' column.
    Then it applies the given plotting function 'plotfunc' to these splits,
    automatically creates folders if non existent and saves the plots to the corresponding dataset folders.
    Examples:
        save_plot_per_dataset(preprocessed_data, rda.lineplots_facet, "../figures/")

        # if an anonymous function (lambda) is used, a plotname has to be provided:
        save_plot_per_dataset(mic_results_long, lambda x: rda.mic_hitstogram(x, "MIC50 in µM"), "../figures/", plotname="MIC_Hits_Distribution")
    """
    if plotname is None:
        plotname = plotfunc.__name__
        if plotname == "<lambda>":
            raise TypeError("Please provide a plotname when using a lambda function.")

    data = data.loc[
        (data["Dataset"] != "Negative Control") & (data["Dataset"] != "Blank")
    ]
    reference_df = data.loc[data["Dataset"] == "Reference"]
    for dataset in filter(lambda x: x != "Reference", data["Dataset"].unique()):
        dataset_data = data.loc[data["Dataset"] == dataset]
        if "AcD Barcode 384" in dataset_data:
            dataset_barcodes = list(dataset_data["AcD Barcode 384"].unique())
            dataset_references = reference_df.loc[
                (reference_df["AcD Barcode 384"].isin(dataset_barcodes)),
                :,
            ]
        else:
            dataset_references = pd.DataFrame()
        set_plot = plotfunc(pd.concat([dataset_data, dataset_references]))
        folder_location = os.path.join(location, dataset)
        pathlib.Path(folder_location).mkdir(parents=True, exist_ok=True)
        for fformat in saveformats:
            print(
                "Saving: ",
                os.path.join(
                    folder_location,
                    f"{dataset}_{plotname}.{fformat}",
                ),
            )
            set_plot.save(
                os.path.join(
                    folder_location,
                    f"{dataset}_{plotname}.{fformat}",
                )
            )


def get_mapping_dict(
    mapping_df: pd.DataFrame,
    mother_column: str = "AsT Barcode 384",
    child_column: str = "AcD Barcode 384",
) -> dict[str, list[str]]:
    mapping_dict = {}
    # TODO: add possibility for MP -> AsT mapping
    for mother_barcode, grp in mapping_df.groupby(mother_column):
        mapping_dict[mother_barcode] = list(grp[child_column].unique())

    return mapping_dict


def lowest_level_dict(mapping_dict):
    # Check if the dict is already on the lowest level (dict of lists)
    if type(list(mapping_dict.values())[0]) is list:
        return mapping_dict
    lowest_mapping_dict = {}
    for d in mapping_dict.values():
        lowest_mapping_dict.update(d)
    # Else recursively reduce the dict
    return lowest_level_dict(lowest_mapping_dict)


def add_precipitation(rawdata, precipitation, mapping_dict):
    if precipitation.empty:
        return rawdata
    precip_all_acd_barcodes = []
    mapping_dict = lowest_level_dict(mapping_dict)
    for acd_barcode, precip_grp in precipitation.groupby("AcD Barcode 384"):
        for parent_barcode, child_barcodes in mapping_dict.items():
            if acd_barcode in child_barcodes:
                for child_barcode in child_barcodes:
                    acd_precip = precip_grp.copy()
                    acd_precip["AcD Barcode 384"] = child_barcode
                    precip_all_acd_barcodes.append(acd_precip)
                # print(acd_barcode, child_barcodes)
    mapped_precipitation = pd.concat(precip_all_acd_barcodes).drop(columns=["Raw Optical Density", "Layout", "Limit of Quantification"])
    return pd.merge(rawdata, mapped_precipitation)


def get_minimum_precipitation_conc(
    grp: pd.DataFrame,
    precip_conc_multiplicator
) -> float | None:
    grp = grp.sort_values("Concentration")
    if grp["Precipitated"].any():
        min_precip_conc = grp["Concentration"][grp["Precipitated"].idxmax()]
        return min_precip_conc * precip_conc_multiplicator
    else:
        return None


def _save_tables(
    resultpath: str, resulttables, fileformats: list[str] = ["xlsx", "csv"]
):
    """
    Save result tables to "<resultpath>".
    Creates corresponding folders for each dataset.
    """
    for result in resulttables:  # cached property of subclasses
        filedir = os.path.join(resultpath, result.dataset)
        pathlib.Path(filedir).mkdir(parents=True, exist_ok=True)
        if "xlsx" or "excel" in fileformats:
            result.table.to_excel(
                os.path.join(filedir, f"{result.file_basename}.xlsx"), index=False
            )
        if "csv" in fileformats:
            result.table.to_csv(
                os.path.join(filedir, f"{result.file_basename}.csv"), index=False
            )


def _save_figures(
    resultpath: str, resultfigures, fileformats: list[str] = ["svg", "html"]
):
    for result in resultfigures:  # cached property of subclasses
        filedir = os.path.join(resultpath, result.dataset)
        pathlib.Path(filedir).mkdir(parents=True, exist_ok=True)
        for file_format in fileformats:
            result.figure.save(
                os.path.join(filedir, f"{result.file_basename}.{file_format}"),
            )
