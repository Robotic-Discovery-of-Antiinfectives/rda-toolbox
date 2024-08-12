#!/usr/bin/env python3

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
