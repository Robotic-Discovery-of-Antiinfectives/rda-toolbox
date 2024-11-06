#!/usr/bin/env python3

import numpy as np
import pandas as pd
import string


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


def generate_inputtable(readout_df = None, platetype: int = 384):
    """
    Generates an input table for the corresponding readout dataframe.
    If not readout df is provided, create a minimal input df.
    """
    if readout_df is None:
        barcodes = ["001PrS01001"]
    else:
        barcodes = readout_df["Barcode"].unique()

    substance_df = pd.DataFrame({
        "ID": [f"Substance {i}" for i in range(1, platetype+1)],
        f"Row_{platetype}": [*list(string.ascii_uppercase[:16]) * 24],
        f"Col_{platetype}": sum([[i] * 16 for i in range(1, 24+1)], []),
        "Concentration in mg/mL": 1,
    })
    layout_df = pd.DataFrame({
        "Barcode": barcodes,
        "Replicate": [1] * len(barcodes),
        "Organism": [f"Placeholder Organism {letter}" for letter in string.ascii_uppercase[:len(barcodes)]]
    })
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
    colmapping = dict(
        zip(list(range(1, 13)), np.array_split(list(range(1, 25)), 12))
    )
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
        df
        .groupby(counts_column)[set_column]
        .apply(lambda x: x.unique())
        .reset_index()
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
            lambda x: "".join(x.split("_")[1:])
            if x.startswith(set_column)
            else x,
            dummies_df.columns,
        )
    )
    # remove any dots as they interfere with altairs plotting.
    dummies_df.columns = dummies_df.columns.str.replace(".", "")
    return dummies_df.drop(columns=["count"])


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
    mapping = {1: 0, 2: 0, 3: 1, 4: 1, 5: 0, 6: 0, 7: 1, 8: 1, 9: 0, 10: 0, 11: 1, 12: 1}

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
