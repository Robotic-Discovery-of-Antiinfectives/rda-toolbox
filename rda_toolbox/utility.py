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

    Takes row, column and quadrant (each of the 96-well plates is one quadrant) of a well from 4 96-well plates and maps it to the corresponding well in a 384-well plate.
    Returns the 384-Well plate row and column.
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
    """
    Apply to a DataFrame the mapping of 96-well positions to 384-well positions.
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
