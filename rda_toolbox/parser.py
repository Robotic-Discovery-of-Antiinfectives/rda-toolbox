#!/usr/bin/env python3
# Data Handling
# Strings
import re

# from pathlib import Path
from io import BytesIO, StringIO

# System
# from os import listdir, makedirs
import os
from os.path import basename  # , exists, isfile, join

# import openpyxl
import numpy as np
import pandas as pd

from .utility import get_rows_cols

# from functools import reduce


# import string


def readerfile_parser(
    filename: str, file_object: StringIO, resulttable_header: str = "Results"
) -> dict:
    """
    Parser for files created by the BioTek Cytation C10 Confocal Imaging Reader.
    """
    lines = file_object.readlines()
    lines = list(filter(None, map(lambda x: x.strip("\n").strip("\r"), lines)))
    if len(lines) == 0:
        raise ValueError(f"Empty raw file {filename}.")

    # search the file for plate type definition and use it to derive number of rows and columns
    found_plate_type = re.findall(r"Plate Type;[A-z ]*([0-9]*)", "".join(lines))
    plate_type = 96  # define default plate type and let it be 96-well plate as this is what we started with
    if found_plate_type:
        plate_type = int(found_plate_type[0])

    num_rows, num_columns = get_rows_cols(plate_type)

    filedict = dict()
    metadata = dict()
    filedict["Reader Filename"] = filename
    filedict["plate_type"] = plate_type
    # TODO: get barcode via regex
    barcode_found = re.findall(
        r"\d{3}[A-Z][a-z][A-Z]\d{2}\d{3}", filedict["Reader Filename"]
    )
    if not barcode_found:
        filedict["Barcode"] = filedict["Reader Filename"]
    else:
        filedict["Barcode"] = barcode_found[0]
    # filedict["Barcode"] = Path(filedict["Reader Filename"]).stem.split("_")[-1]

    results = np.empty([num_rows, num_columns], dtype=float)
    # using dtype=str results in unicode strings of length 1 ('U1'), therefore we use 'U25'
    layout = np.empty([num_rows, num_columns], dtype="U25")
    concentrations = np.empty([num_rows, num_columns], dtype=float)

    metadata_regex = r";?([a-zA-Z0-9 \/]*)[;:]+([a-zA-Z0-9 \/\\:_.-]*),?"
    line_num = 0
    while line_num < len(lines):
        if lines[line_num] == resulttable_header:
            line_num += 1
            header = map(
                int, lines[line_num].strip("\n").split(";")[1:]
            )  # get the header
            index = [""] * num_rows
            for _row_num in range(num_rows):  # for the next num_rows, read result data
                line_num += 1
                res_line = lines[line_num].split(";")
                # Split at ; and slice off rowlabel and excitation/emission value:
                index[_row_num] = res_line[0]
                results[_row_num] = res_line[1:-1]
            # Initialize DataFrame from results and add it to filedict
            filedict["Raw Optical Density"] = pd.DataFrame(
                data=results, index=index, columns=header
            )
            line_num += 1
        elif lines[line_num] == "Layout":  # For the next num_rows, read layout data
            line_num += 1
            header = list(
                map(int, lines[line_num].strip("\n").split(";")[1:])
            )  # Because we use header twice here, we collect it via list()
            index = [""] * num_rows
            for _row_num in range(num_rows):
                line_num += 1
                layout_line = lines[line_num].split(";")
                index[_row_num] = layout_line[0]
                layout[_row_num] = layout_line[1:-1]
                # Each second line yields a concentration layout line
                line_num += 1
                conc_line = lines[line_num].split(";")
                concentrations[_row_num] = [
                    None if not x else float(x) for x in conc_line[1:-1]
                ]
            # Add layouts to filedict
            filedict["Layout"] = pd.DataFrame(data=layout, index=index, columns=header)
            filedict["Concentration"] = pd.DataFrame(
                data=concentrations, index=index, columns=header
            )
            line_num += 1
        else:
            metadata_pairs = re.findall(metadata_regex, lines[line_num])
            line_num += 1
            if not metadata_pairs:
                continue
            else:
                for key, value in metadata_pairs:
                    if not all(
                        [key, value]
                    ):  # if any of the keys or values are empty, skip
                        continue
                    else:
                        metadata[key.strip(" :")] = value.strip(" ")
    filedict["metadata"] = metadata
    return filedict


def filepaths_to_filedicts(filepaths: list[str]) -> list[dict]:
    """
    Wrapper function to obtain a list of dictionaries which contain the raw files information like

    - different entries of metadata
        - Plate Type
        - Barcode
        - Date
        - Time
        - etc.
    - Raw Optical Density (DataFrame)
    - Concentration (DataFrame)
    - Layout (DataFrame)
    """
    filedicts = []
    for path in filepaths:
        file = open(path)
        contents = StringIO(file.read())
        filedicts.append(readerfile_parser(basename(path), contents))
        file.close()
    return filedicts


def collect_metadata(filedicts: list[dict]) -> pd.DataFrame:
    """
    Helperfunction to collect the metadata from all reader files into a dataframe.
    """
    allmetadata_df = pd.DataFrame()
    for filedict in filedicts:
        meta_df = pd.DataFrame(filedict["metadata"], index=[0])
        meta_df["Barcode"] = filedict["Barcode"]
        allmetadata_df = pd.concat([allmetadata_df, meta_df], ignore_index=True)
    return allmetadata_df


def collect_results(filedicts: list[dict]) -> pd.DataFrame:
    """
    Collect and merge results from the readerfiles.
    """
    allresults_df = pd.DataFrame(
        {"Row": [], "Column": [], "Raw Optical Density": []}
    )  # , "Layout": [], "Concentration": []})
    platetype_s = list(set(fd["plate_type"] for fd in filedicts))
    if len(platetype_s) == 1:
        platetype = platetype_s[0]
    else:
        raise Exception(f"Different plate types used {platetype_s}")

    for filedict in filedicts:
        # long_layout_df = get_long_df("Layout")
        # long_concentrations_df = get_long_df("Concentration")
        # long_rawdata_df = get_long_df("Raw Optical Density")

        long_rawdata_df = pd.melt(
            filedict["Raw Optical Density"].reset_index(names="Row"),
            id_vars=["Row"],
            var_name="Column",
            value_name="Raw Optical Density",
        )

        long_rawdata_df["Barcode"] = filedict["Barcode"]
        # df_merged = reduce(
        #     lambda  left,right: pd.merge(left,right,on=['Row', 'Column'], how='outer'),
        #     [long_rawdata_df, long_layout_df, long_concentrations_df]
        # )
        allresults_df = pd.concat([allresults_df, long_rawdata_df], axis=0)
        platetype = filedict["plate_type"]

    allresults_df.rename(
        columns={"Row": f"Row_{platetype}", "Column": f"Col_{platetype}"}, inplace=True
    )
    return allresults_df.reset_index(drop=True)


def parse_readerfiles(path: str | None) -> pd.DataFrame | None:
    """
    Reads CytationC10 readerfiles (plain text files) and merges the results into a DataFrame which is returned.
    Wrapper for readerfiles_rawdf to keep backwards compatibility.
    Improves readerfiles_rawdf, provide a single path for convenience.
    """
    if not path:
        return None
    paths = [
            os.path.join(path, f)
            for f in os.listdir(path)
            if os.path.isfile(os.path.join(path, f))
    ]
    df = readerfiles_rawdf(paths)
    df["Col_384"] = df["Col_384"].astype(int)
    return df

def readerfiles_rawdf(paths: list[str]) -> pd.DataFrame:
    """Parses data from files declared by filepaths and merges the results into a DataFrame
    :param paths: A list of filepaths corresponding to the raw reader files generated by Cytation10
    :type paths: list[str]
    :return: A DataFrame in tidy and long format with the raw readerfile contents
    :rtype: pd.DataFrame

    :Example:

        ```Python
        import glob

        rawdata_df = readerfiles_rawdf(glob.glob("path/to/raw/files/*"))
        ```
    """
    filedicts = filepaths_to_filedicts(paths)
    rawdata = collect_results(filedicts)
    rawdata["Col_384"] = rawdata["Col_384"].astype(str)
    rawdata.rename(columns={"Barcode": "AcD Barcode 384"}, inplace=True)
    return rawdata


def readerfiles_metadf(paths: list[str]) -> pd.DataFrame:
    """
    Parses metadata from files declared by filepaths and merges the results into a DataFrame.
    """
    filedicts = filepaths_to_filedicts(paths)
    return collect_metadata(filedicts)


def process_inputfile(file_object):
    """
    Read Input excel file which should have the following columns:
        - Barcode
        - Organism
        - Row_384
        - Col_384
        - ID
    Optional columns:
        - Concentration in mg/mL (or other units)
        - Cutoff
    """
    if not file_object:
        return None
    excel_file = pd.ExcelFile(file_object)
    substance_df = pd.read_excel(excel_file, "substances")
    layout_df = pd.read_excel(excel_file, "layout")
    df = pd.merge(layout_df, substance_df, how="cross")
    # df.rename(columns={
    #     "barcode": "Barcode",
    #     "replicate": "Replicate",
    #     "organism": "Organism",
    #     "plate_row": "Row_384",
    #     "plate_column": "Col_384",
    #     "id": "ID",
    #     "concentration": "Concentration in mg/mL",
    # }, inplace=True)
    df["ID"] = df["ID"].astype(str)
    return df


def read_platemapping(filepath: str, orig_barcodes: list[str]):
    """
    Reads a mappingfile generated by the barcode reader.

    """
    filedict = dict()
    orig_barcodes = list(map(str, orig_barcodes))
    with open(filepath) as file:
        filecontents = file.read().splitlines()
        origin_barcode = ""
        origin_replicates = []
        for line in filecontents:
            line = line.split(";")
            if len(line) == 1 and line[0] in orig_barcodes:
                origin_barcode = line[0]
                origin_replicates.append(origin_barcode)
                # print("Origin barcode: ", origin_barcode)
                if origin_barcode not in filedict:
                    filedict[origin_barcode] = []
            else:
                filedict[origin_barcode].append(line)
        replicates_dict = {i:origin_replicates.count(i) for i in origin_replicates}
        if sorted(list(filedict.keys())) != sorted(orig_barcodes):
            raise ValueError(
                f"The origin barcodes from the mappingfile '{os.path.basename(filepath)}' and MP barcodes in MIC_input.xlsx do not coincide."
            )
        return filedict, replicates_dict


def parse_mappingfile(
    filepath: str,
    motherplate_column: str = "Origin Plate",
    childplate_column: str = "AcD Barcode 384",
):
    """
    Simple mappingfile parser function.
    Expects to start with a "Motherplate" line followed by corresponding "Childplates" in a single line.
    """
    filedict = dict()
    with open(filepath) as file:
        filecontents = file.read().splitlines()
        key = None
        for line in filecontents:
            line = line.split(";")
            if len(line) == 1:
                key = line[0]
            else:
                if not key:
                    raise ValueError(
                        "Motherplate barcode expected on first line."
                    )
                if key in filedict:
                    filedict[key].append(line)
                else:
                    filedict[key] = [line]
    mapping_df = pd.DataFrame(
        [
            (motherplate, childplate, rep_num, rack_nr)
            for motherplate, replicates in filedict.items()
            for rep_num, childplates in enumerate(replicates, start=1)
            for rack_nr, childplate in enumerate(childplates, start=1)
        ],
        columns=[motherplate_column, childplate_column, "Replicate", "Rack"],
    )
    return mapping_df


def read_inputfile(inputfile_path: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:

    dtypes = { # define type dict to read the correct types from excel
        'Internal ID': str,
        'PlateNr 96': int,
        'MP Barcode 96': str,
        'Position 96': str,
        'Row 96': str,
        'Col 96': int,
        'PlateNr 384': int,
        'AsT Barcode 384': str,
        'Quadrant': int,
        'Dataset': str,
        'Row 384': str,
        'Col 384': int,
        'Rack': int,
        'Organism': str,
        # 'Row_384': str,
        # 'Col_384': int,
    }
    substances = pd.read_excel(
        inputfile_path,
        sheet_name="Substances",
        dtype=dtypes,
    )
    organisms = pd.read_excel(inputfile_path, sheet_name="Organisms", dtype=dtypes)
    dilutions = pd.read_excel(inputfile_path, sheet_name="Dilutions", dtype=dtypes)
    controls = pd.read_excel(inputfile_path, sheet_name="Controls", dtype=dtypes)

    controls["Row_384"] = controls["Position"].apply(lambda x: str(x[0]))
    controls["Col_384"] = controls["Position"].apply(lambda x: int(x[1:]))
    controls.drop(columns="Position", inplace=True)

    return substances, organisms, dilutions, controls
