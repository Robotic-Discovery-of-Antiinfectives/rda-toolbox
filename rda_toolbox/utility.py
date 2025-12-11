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

from io import BytesIO
import xlsxwriter
import re

import altair as alt

from collections.abc import Sequence
from typing import Tuple, Any


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
    df[col] = df[position].apply(lambda x: int(str(x[1:])))
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

def mic_assaytransfer_mapping(
    position: str,
    orig_barcode: Any,
    ast_platemapping: dict,
    *,
    strict: bool = False,
) -> Tuple[str, str, str]:
    """
    Map a 96-well motherplate position to a 384-well AsT plate position.

    Parameters
    ----------
    position : str
        Well position on the 96-well plate in the form 'A1'..'H12'.
    orig_barcode : Any
        Identifier of the 96-well motherplate. Will be cast to str.
    ast_platemapping : mapping
        Mapping from motherplate barcode -> sequence of AsT plate barcodes.
        It should support `ast_platemapping[orig_barcode]`.

        Common layouts:
        - dict[str, Sequence[str]]
        - pandas.Series/row where `row[0]` holds a Sequence[str]

    strict : bool, default False
        If True, raise a ValueError when the requested third of the plate
        (1/3, 2/3, 3/3) does not exist in `ast_platemapping`.
        If False, fall back to the last available AsT plate and continue.

    Returns
    -------
    (row_384, col_384, ast_barcode) : (str, str, str)
        384-well plate row (A–P), column (1–24, here only 1–2 used),
        and the corresponding AsT plate barcode.

    Notes
    -----
    - 96-well plate: rows A–H, cols 1–12.
    - 384-well plate: rows A–P, cols 1–24.
      This function maps each 96-well to one of two 384 rows (2x upscaling)
      and to column 1 or 2 (alternating).
    - The 96-well plate is conceptually split into three vertical thirds:
        cols  1–4  -> AsT plate index 0
        cols  5–8  -> AsT plate index 1
        cols  9–12 -> AsT plate index 2
    """

    # ---- Normalize and validate input ----
    if not isinstance(position, str):
        raise TypeError(f"position must be a string like 'A1', got {type(position)}")

    position = position.strip().upper()
    if len(position) < 2:
        raise ValueError(f"Invalid well position {position!r}")

    row = position[0]
    col_str = position[1:]

    if row not in "ABCDEFGH":
        raise ValueError(f"Row {row!r} out of range for 96-well plate (A–H).")

    try:
        col = int(col_str)
    except ValueError as exc:
        raise ValueError(f"Column {col_str!r} is not an integer in position {position!r}.") from exc

    if not (1 <= col <= 12):
        raise ValueError(f"Column {col} out of range for 96-well plate (1–12).")

    orig_barcode = str(orig_barcode)

    # ---- Build row mapping (A–H -> pairs like [A,B], [C,D], ..., [O,P]) ----
    # 384 rows we use: A–P (16 rows), grouped into 8 pairs
    rows_384 = list(string.ascii_uppercase[:16])  # ['A', ..., 'P']
    row_pairs = [rows_384[i : i + 2] for i in range(0, 16, 2)]  # [['A','B'], ['C','D'], ..., ['O','P']]

    row_index_96 = "ABCDEFGH".index(row)

    # Equivalent to your mapping dict:
    # mapping = {1:0, 2:0, 3:1, 4:1, 5:0, 6:0, 7:1, 8:1, 9:0, 10:0, 11:1, 12:1}
    # This can be expressed as:
    mapping_idx = ((col - 1) // 2) % 2  # gives 0,0,1,1,0,0,1,1,0,0,1,1 for col 1..12

    row_384 = row_pairs[row_index_96][mapping_idx]

    # ---- Column mapping (1–12 -> 1 or 2 alternating) ----
    # Your original colmapping: 1->1, 2->2, 3->1, 4->2, ...
    col_384 = 1 if (col % 2) == 1 else 2

    # ---- Determine which 1/3 of the motherplate this is in (0,1,2) ----
    # 1–4 -> 0; 5–8 -> 1; 9–12 -> 2
    ast_of_3 = (col - 1) // 4

    # ---- Normalize access to ast_platemapping ----
    def _get_ast_plates(mapping_obj, barcode: str) -> list[str]:
        """Return a list of AsT barcodes for a given motherplate barcode."""
        try:
            entry = mapping_obj[barcode]
        except KeyError as exc:
            raise KeyError(
                f"No entry for motherplate barcode {barcode!r} in ast_platemapping."
            ) from exc

        # For pandas row / Series where first column holds the list
        # try to access [0]; if that fails, use entry directly.
        try:
            candidate = entry[0]
        except Exception:
            candidate = entry

        # Normalize to list[str]
        if isinstance(candidate, (str, bytes)):
            plates = [candidate]
        elif isinstance(candidate, Sequence):
            plates = list(candidate)
        else:
            # Last resort: wrap in list
            plates = [str(candidate)]

        # Filter out obvious empties
        plates = [str(p) for p in plates if p not in (None, "", "nan")]
        if not plates:
            raise ValueError(f"ast_platemapping[{barcode!r}] contains no valid AsT barcodes.")
        return plates

    ast_plates = _get_ast_plates(ast_platemapping, orig_barcode)

    # ---- Choose the correct AsT plate, safely ----
    if ast_of_3 >= len(ast_plates):
        # We are requesting e.g. the 3rd third (index 2) but only 1 or 2 AsT plates exist.
        if strict:
            raise ValueError(
                f"Motherplate {orig_barcode!r} has only {len(ast_plates)} AsT plate(s), "
                f"cannot select segment index {ast_of_3} for column {col}."
            )
        # Non-strict mode: fall back to the last available AsT plate
        # and keep going, but at least make it explicit.
        ast_index = len(ast_plates) - 1
    else:
        ast_index = ast_of_3

    barcode_384_ast = ast_plates[ast_index]

    return str(row_384), str(col_384), str(barcode_384_ast)


# def mic_assaytransfer_mapping(position, orig_barcode, ast_platemapping):
#     """
#     This is a rather unfinished function to map 96 well motherplates to 384 well assay transfer (AsT) plates.

#     """
#     row = position[0]
#     col = int(position[1:])
#     orig_barcode = str(orig_barcode)
#     rowmapping = dict(
#         zip(
#             string.ascii_uppercase[0:8],
#             np.array_split(list(string.ascii_uppercase)[0:16], 8),
#         )
#     )
#     colmapping = dict(zip(list(range(1, 13)), [1, 2] * 13))
#     mapping = {
#         1: 0,
#         2: 0,
#         3: 1,
#         4: 1,
#         5: 0,
#         6: 0,
#         7: 1,
#         8: 1,
#         9: 0,
#         10: 0,
#         11: 1,
#         12: 1,
#     }

#     row_384 = rowmapping[row][mapping[col]]
#     col_384 = colmapping[col]

#     # TODO: sometimes only the middle (5-8) or the last part of a motherplate is taken...
#     # Currently this would lead to an index error since A5 -> ast_of_3 = 1 but if this is the only AsT plate, the first AsT plate is "missing"...
#     # Possible solution would be to try - except decreasing the position by 4 iteratively...
#     # try A5 except IndexError: A5 - 4 -> try A1 success etc.
#     # or A12 -> ast_of_3 = 2 -> IndexError. A12 - 4 -> A8 -> IndexError. A8 - 4 -> A4 works...
#     # seems like a bad solution though.
#     # num_of_ast_plates = len(ast_platemapping[orig_barcode][0])
#     # if num_of_ast_plates < 3:
#     #     print("Did someone just take one third of a motherplate (Only one or two AsT plates for one MP)?")

#     if col in [1, 2, 3, 4]:
#         ast_of_3 = 0
#     elif col in [5, 6, 7, 8]:
#         ast_of_3 = 1
#     else:
#         ast_of_3 = 2
#     # print("orig_barcode", orig_barcode, "ast_of_3", ast_of_3)
#     print(position, orig_barcode, ast_platemapping)
#     print("ast_platemapping", ast_platemapping)
#     print("col", col)
#     # print(ast_platemapping[orig_barcode][0])
#     barcode_384_ast = ast_platemapping[orig_barcode][0][ast_of_3]
#     return str(row_384), str(col_384), barcode_384_ast


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


def write_excel_MolImages(df: pd.DataFrame, filename: str, molcol_header: str):
    """
    Writes images (.png) of molecules structures into an excel file derived from the given dataframe.
    """
    if molcol_header not in df:
        raise ValueError(
            f"Missing {molcol_header} column in df."
        )

    df["img_buf"] = df[molcol_header].apply(lambda x: mol_to_bytes(x) if x else None)
    writer = pd.ExcelWriter(filename, engine="xlsxwriter")

    workbook = writer.book
    worksheet = workbook.add_worksheet()
    for i, imgbuf in enumerate(list(df["img_buf"]), start=1):
        worksheet.set_column(0, 0, 20)
        worksheet.set_row(i, 120)
        worksheet.insert_image(
            f"A{i+1}", "img.png", {"image_data": imgbuf, "x_scale": 0.5, "y_scale": 0.5}
        )

    df.drop(columns=["img_buf", molcol_header]).to_excel(writer, startcol=1, index=False)
    workbook.close()


def smiles_to_imgstr(smiles):
    """
    Converts a smiles string to a base64 encoded image string (e.g. for plotting in altair).
    It's a convenience function consisting of rda.utility.mol_to_bytes() and rda.utility.imgbuffer_to_imgstr(),
    use these if you want more fine grained control over the format of the returned string.
    Example: df["image"] = df["smiles"].apply(lambda smiles: smiles_to_imgstr(smiles))
    """
    return imgbuffer_to_imgstr(mol_to_bytes(Chem.MolFromSmiles(smiles)))


def inchi_to_imgstr(inchi):
    """
    Converts a inchi string to a base64 encoded image string (e.g. for plotting in altair).
    It's a convenience function consisting of rda.utility.mol_to_bytes() and rda.utility.imgbuffer_to_imgstr(),
    use these if you want more fine grained control over the format of the returned string.
    Example: df["image"] = df["inchi"].apply(lambda inchi: inchi_to_imgstr(inchi))
    """
    return imgbuffer_to_imgstr(mol_to_bytes(Chem.MolFromInchi(inchi)))


def prepare_visualization(
    df,
    by_id="Internal ID",
    whisker_width=1,
    exclude_negative_zfactors=True,
    threshold=50.0,
):
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
        # print(grp)
        maxconc_below_threshold = (
            grp[
                (grp["Replicate"] == 1)
                & (grp["Concentration"] == grp["Concentration"].max())
            ]["Mean Relative Optical Density"]
            < threshold
        )
        grp["max_conc_below_threshold"] = list(maxconc_below_threshold)[0]
        tmp_list.append(grp)
        # .sort_values(by=["Concentration"], ascending=False)
        # grp_sorted[grp_sorted["Concentration"] == 50]["Mean Relative Optical Density"]
        # print(grp.aggregate())
    df = pd.concat(tmp_list)
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
    plotname: str | None = None,
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
    mapped_precipitation = pd.concat(precip_all_acd_barcodes).drop(
        columns=["Raw Optical Density", "Layout", "Limit of Quantification"]
    )
    return pd.merge(rawdata, mapped_precipitation, how="outer")


def get_minimum_precipitation_conc(
    grp: pd.DataFrame, precip_conc_multiplicator
) -> float | None:
    grp = grp.sort_values("Concentration")
    if grp["Precipitated"].any():
        # Get the first (minimum) concentration at which precipitation happened:
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
        multiindex = False
        if result.table.columns.nlevels > 1:  # Check whether multiple levels are on the columns
            multiindex = True
        if "xlsx" or "excel" in fileformats:
            result.table.to_excel(
                os.path.join(filedir, f"{result.file_basename}.xlsx"), index=multiindex
            )
        if "csv" in fileformats:
            result.table.to_csv(
                os.path.join(filedir, f"{result.file_basename}.csv"), index=multiindex
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


def to_excel_molimages(
    df: pd.DataFrame, filename: str, desired_columns: list[str], mol_col: str = "mol"
):
    """
    Writes a dataframe containing RDKit molecule objects to an excel file containing the molecular structures as PNG images.
    Needs a column in df with RDKit mol object (e.g. rdkit.Chem.MolFromInchi, MolFromMolBlock, MolFromSmiles etc.)
    """
    writer = pd.ExcelWriter(filename, engine="xlsxwriter")
    workbook = writer.book
    # workbook = xlsxwriter.Workbook("images_bytesio.xlsx")
    worksheet = workbook.add_worksheet()
    df["PIL_img"] = df[mol_col].map(Draw.MolToImage)

    def get_imgbuffer(image):
        stream = BytesIO()
        image.save(stream, format="PNG")
        return stream

    for i, imgbuf in enumerate(df["PIL_img"].map(get_imgbuffer), start=1):
        worksheet.set_column(0, 0, 20)
        worksheet.set_row(i, 120)
        worksheet.insert_image(
            f"A{i+1}", "img.png", {"image_data": imgbuf, "x_scale": 0.5, "y_scale": 0.5}
        )

    df.loc[:, desired_columns].to_excel(writer, startcol=1, index=False)
    workbook.close()


def check_activity_conditions(
        organism_measurements,
        organism_bscores,
        measurement_threshold,
        bscore_threshold,
        ):
    check_conditions = []
    for i, value in enumerate(organism_measurements):
        if value < measurement_threshold:  # If below growth threshold:
            if (
                    # Check if the sample which seems 'active' for ONE organism has a viable B-Score for this exact organism
                    list(organism_bscores)[i] < bscore_threshold
                    ):  # AND if below B-Score threshold
                check_conditions.append(True)
            else:
                check_conditions.append(False)
        else:
            check_conditions.append(False)
    # print(check_conditions)
    return any(
            check_conditions
            )  # Any because we want the values for other organisms if sample is active in one


def format_organism_name(raw_organism_name: str) -> str:
    """
    Create internal, formatted orgnames but keep external for plots...

    """
    # I dont know how many spaces people introduce...
    normalized_spaces = re.sub(r'\s+', ' ', raw_organism_name).strip()
    # I dont know which organism names people think of... (MRSA ST033793 vs. Acinetobacter Baumannii ATCC 17978)
    # Or which other evil things people come up with
    return normalized_spaces.lower()


def read_sdf_withproperties(sdf_filepath: str) -> pd.DataFrame:
    """
    Reads a SDF file and returns a DataFrame containing the molecules as rdkit molobjects
    as well as all the encoded properties in the SDF block.
    """
    suppl = Chem.SDMolSupplier(sdf_filepath)
    mols = []
    nonmol_counter = 0
    for mol in suppl:
        if not mol:
            nonmol_counter += 1
            continue
        mol_props = {"mol": mol}
        propnames = mol.GetPropNames()
        for prop in propnames:
            mol_props[prop] = mol.GetProp(prop)
        mols.append(mol_props)
    if nonmol_counter > 0:
        print(f"Ignored molecules: {nonmol_counter}")
    return pd.DataFrame(mols)



def smiles_grid_altair(
        df: pd.DataFrame,
        smiles_col: str = "smiles",
        n_cols: int = 6,
        img_size: int = 128,
        tooltip_cols: list[str] | None = None,
        drop_invalid: bool = True,
        background: str | None = '#ffffff',
        gridtitle: str | None = "Molecule Grid",
        ):
    """
    Render a grid of molecule images from a DataFrame using Altair.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain a SMILES column (default 'smiles'). Other columns are arbitrary.
    smiles_col : str
        Column name with SMILES strings.
    n_cols : int
        Number of grid columns.
    img_size : int
        Size (pixels) of each PNG image (square).
    tooltip_cols : list[str] | None
        If provided, use only these columns as tooltips. Otherwise include all non-internal columns.
    drop_invalid : bool
        If True, rows with invalid/unparsable SMILES are removed. If False, they'll be kept with empty images.
    background : str | None
        Optional CSS color (e.g., '#ffffff') for chart background.

    Returns
    -------
    alt.Chart
    """

    if smiles_col not in df.columns:
        raise ValueError(f"Column '{smiles_col}' not found in DataFrame.")

    data = df.copy()
    # --- helper: SMILES -> data URL (PNG in-memory) ---
    def _smiles_to_data_url(smi: str) -> str | None:
        try:
            mol = Chem.MolFromSmiles(str(smi))
            if mol is None:
                return None
            # RDKit renders to PIL Image; save to BytesIO as PNG
            img = Draw.MolToImage(mol)
            buf = BytesIO()
            img.save(buf, format="png")
            encoded = base64.b64encode(buf.getvalue()).decode()
            return f"data:image/png;base64,{encoded}"
        except Exception:
            print("Exception: ", smi)
            return None

    # Create image column
    data["_image_url"] = data[smiles_col].apply(_smiles_to_data_url)

    # Handle invalids
    if drop_invalid:
        data = data.loc[data["_image_url"].notna()].copy()

    if len(data) == 0:
        raise ValueError("No valid molecules to render (all SMILES failed to parse?).")

    # Grid coordinates
    data = data.reset_index(drop=True)
    data["_idx"] = np.arange(len(data))
    data["_col"] = (data["_idx"] % n_cols).astype(int)
    data["_row"] = (data["_idx"] // n_cols).astype(int)

    # Tooltips: include all user columns by default (exclude internal helpers)
    internal_cols = {"_image_url", "_idx", "_col", "_row", "smiles"}
    if tooltip_cols is None:
        tooltip_cols = [c for c in data.columns if c not in internal_cols]

    # Build chart
    chart = (
            alt.Chart(data)
            .mark_image(width=img_size, height=img_size)
            .encode(
                x=alt.X("_col:O", axis=None),
                # reverse rows so row 0 is at the top
                y=alt.Y("_row:O", axis=None, sort="ascending"),
                url="_image_url:N",
                tooltip=tooltip_cols,
                )
            .properties(
                width=n_cols * img_size,
                height=(int(np.ceil(len(data) / n_cols))) * img_size,
                background=background,
                title=gridtitle,
                )
            )

    return chart
