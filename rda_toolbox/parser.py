#!/usr/bin/env python3
# Data Handling
# Strings
import re

# from pathlib import Path
from typing import IO

# System
# from os import listdir, makedirs
import os
from os.path import basename  # , exists, isfile, join

# import openpyxl
import numpy as np
import pandas as pd

from .utility import get_rows_cols, format_organism_name, position_to_rowcol

# from functools import reduce


# import string


def readerfile_parser(
    filename: str, file_object: IO[str], resulttable_header: str = "Results"
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
        r"\d{3}[A-Z][a-z]?[a-zA-Z]\d{2}\d{3}", filedict["Reader Filename"]
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
            header = list(
                map(int, lines[line_num].strip("\n").split(";")[1:])
            )  # get the header as a concrete list
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
        try:
            with open(path, encoding="utf-8", errors="ignore") as fh:
                filedicts.append(readerfile_parser(basename(path), fh))
        except OSError as exc:
            raise OSError(f"Failed to read {path!r}: {exc}") from exc
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


def parse_readerfiles(path: str | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reads CytationC10 readerfiles (plain text files) and merges the results into
    two DataFrames (rawdata and metadata) which is returned.
    Wrapper for readerfiles_rawdf to keep backwards compatibility.
    Improves readerfiles_rawdf, provide a single path for convenience.
    """
    if not path:
        return pd.DataFrame(), pd.DataFrame()
    paths = [
            os.path.join(path, f)
            for f in os.listdir(path)
            if os.path.isfile(os.path.join(path, f))
    ]
    df_raw = readerfiles_rawdf(paths)
    df_raw["Col_384"] = df_raw["Col_384"].astype(int)
    df_meta = readerfiles_metadf(paths)
    return df_raw, df_meta

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
        for i, line in enumerate(filecontents):
            line = line.split(";")
            if i % 2 == 0:  # if i is even (expect MPs on even lines, alternating with childplates)
            # if len(line) == 1:
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


def _validate_inputfile_structure(inputfile_path: str, substance_id: str) -> None:
    """
    Validate that the input Excel file exists and contains the expected sheets
    and minimal required columns. Raises ValueError with a human-readable list
    of issues if anything important is missing.

    Checks performed:
    - file exists
    - required sheets: Substances, Organisms, Dilutions, Controls
    - Substances sheet contains the column given by `substance_id`
    - Organisms sheet contains 'Organism'
    - Controls sheet contains a column starting with 'Position' and all non-null
      values in that column look like a letter followed by digits (e.g. 'A1', 'P24')
    """
    issues: list[str] = []

    if not inputfile_path or not os.path.exists(inputfile_path):
        raise FileNotFoundError(f"Input file not found: {inputfile_path!r}")

    try:
        xls = pd.ExcelFile(inputfile_path)
    except Exception as exc:
        raise ValueError(f"Unable to open Excel file {inputfile_path!r}: {exc}") from exc

    available_sheets = [str(s).strip() for s in xls.sheet_names]
    required_sheets = ["Substances", "Organisms", "Dilutions", "Controls"]
    missing_sheets = [s for s in required_sheets if s not in available_sheets]
    if missing_sheets:
        issues.append(
            "Missing sheets: " + ", ".join(missing_sheets) + ". "
            "Expected sheets are: Substances, Organisms, Dilutions, Controls."
        )

    # If Substances sheet exists, check for substance_id column
    if "Substances" in available_sheets:
        try:
            subs_df = pd.read_excel(xls, "Substances", nrows=0)
        except Exception:
            issues.append("Could not read the 'Substances' sheet.")
        else:
            if substance_id not in subs_df.columns:
                issues.append(
                    f"Substances sheet does not contain the required column {substance_id!r}."
                    " This column identifies each substance (e.g. an internal ID)."
                )
            else:
                # attempt to read some rows to validate IDs
                try:
                    subs_df = pd.read_excel(xls, "Substances")
                except Exception:
                    issues.append("Could not read the contents of the 'Substances' sheet.")
                else:
                    if subs_df[substance_id].isnull().any():
                        issues.append(
                            f"Substances.{substance_id} contains empty values. Every substance must have an ID."
                        )
                    dup_mask = subs_df[substance_id].duplicated(keep=False)
                    if dup_mask.any():
                        example_dups = subs_df.loc[dup_mask, substance_id].astype(str).unique()[:5]
                        issues.append(
                            f"Substances.{substance_id} contains duplicate IDs. Example duplicates: {', '.join(map(str, example_dups))}."
                        )
                    

    # Organisms sheet -> must have 'Organism'
    if "Organisms" in available_sheets:
        try:
            org_df = pd.read_excel(xls, "Organisms", nrows=0)
        except Exception:
            issues.append("Could not read the 'Organisms' sheet.")
        else:
            if "Organism" not in org_df.columns:
                issues.append(
                    "Organisms sheet does not contain the required column 'Organism'."
                )
            else:
                try:
                    org_df = pd.read_excel(xls, "Organisms", nrows=50)
                except Exception:
                    issues.append("Could not read the contents of the 'Organisms' sheet.")
                else:
                    if org_df["Organism"].dropna().empty:
                        issues.append("Organisms sheet appears to be empty — at least one Organism entry is required.")
    # --- Dilutions sheet checks ---
    if "Dilutions" in available_sheets:
        try:
            dil_header = pd.read_excel(xls, "Dilutions", nrows=0)
        except Exception:
            issues.append("Could not read the 'Dilutions' sheet header.")
        else:
            try:
                dil_df = pd.read_excel(xls, "Dilutions", nrows=50)
            except Exception:
                issues.append("Could not read the contents of the 'Dilutions' sheet.")
            else:
                # look for at least one numeric concentration-like column name
                name_like = [
                    c for c in dil_df.columns
                    if re.search(r"conc|concent|dilut|dose", str(c), re.I)
                ]
                # detect numeric columns by dtype as fallback
                numeric_cols = [c for c in dil_df.columns if pd.api.types.is_numeric_dtype(dil_df[c])]
                if not name_like and not numeric_cols:
                    issues.append(
                        "Dilutions sheet does not appear to contain any concentration/dilution columns. "
                        "Expected a column with concentration/dilution values (name containing 'conc'/'dilut' or numeric values)."
                    )
                else:
                    # decide which columns to treat as concentration columns (prefer name_like)
                    conc_cols = name_like if name_like else numeric_cols[:1]

                    # check for units: either in column header (e.g. "Concentration (mg/mL)" or "conc_mg_per_ml")
                    unit_header_pattern = re.compile(
                        r"\b(mg/?ml|ug/?ml|µg/?ml|ng/?ml|g/?ml|mM|uM|µM|M|mol/?L|mmol/?L|%|per ?ml)\b",
                        re.I,
                    )
                    header_has_unit = any(bool(unit_header_pattern.search(str(c))) for c in conc_cols)

                    # or a dedicated Unit/Units column with at least one non-empty value
                    unit_col_candidates = [c for c in dil_df.columns if str(c).strip().lower() in ("unit", "units", "concentration unit")]
                    unit_values_present = False
                    if unit_col_candidates:
                        uc = unit_col_candidates[0]
                        unit_values_present = dil_df[uc].dropna().astype(str).str.strip().any()

                    if not (header_has_unit or unit_values_present):
                        issues.append(
                            "Dilutions sheet: no concentration unit detected for the concentration column(s). "
                            "Please include units (e.g. 'mg/mL' or 'mM') either in the column header "
                            "(e.g. 'Concentration (mg/mL)') or add a 'Unit' column with values. "
                        )
    # Controls sheet -> must have a 'Position*' column
    if "Controls" in available_sheets:
        try:
            ctrl_df = pd.read_excel(xls, "Controls")
        except Exception:
            issues.append("Could not read the 'Controls' sheet.")
        else:
            poscols = [c for c in ctrl_df.columns if str(c).startswith("Position")]
            if not poscols:
                issues.append(
                    "Controls sheet must contain a column that starts with 'Position' "
                    "(e.g. 'Position 96' or 'Position 384')."
                )
            else:
                poscol = poscols[0]
                # Validate format of position entries: letter(s) followed by digits
                invalid_positions = []
                for i, val in enumerate(ctrl_df[poscol].dropna().astype(str), start=1):
                    if len(val) < 2:
                        invalid_positions.append((i, val))
                        continue
                    row_letter = val[0]
                    col_part = val[1:]
                    if not row_letter.isalpha() or not col_part.isdigit():
                        invalid_positions.append((i, val))
                if invalid_positions:
                    sample = ", ".join(f"{idx}:{v!r}" for idx, v in invalid_positions[:5])
                    issues.append(
                        f"Controls.{poscol} contains entries that are not in the expected "
                        f"format (letter + digits). Sample invalid entries (row:index:value): {sample}. "
                        "Positions should look like 'A1' or 'P24'."
                    )

    if issues:
        raise ValueError(
            "Input file validation failed. Please fix the following issues in the Excel file:\n- "
            + "\n- ".join(issues)
        )
    

def read_inputfile(inputfile_path: str, substance_id) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
# Validate file structure and content before attempting to parse
    _validate_inputfile_structure(inputfile_path, substance_id)

    dtypes = { # define type dict to read the correct types from excel
        substance_id: str,
        'PlateNr 96': str,  # This could be Int, but lab members chose alphabetic platenumbers (in addition)
        'MP Barcode 96': str,
        'Position 96': str,
        'Row 96': str,
        'Col 96': int,
        'PlateNr 384': str,
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
    ).rename(columns={substance_id: "Internal ID"})

    organisms = pd.read_excel(inputfile_path, sheet_name="Organisms", dtype=dtypes)
    organisms["Organism formatted"] = organisms["Organism"].apply(format_organism_name)
    dilutions = pd.read_excel(inputfile_path, sheet_name="Dilutions", dtype=dtypes)
    controls = pd.read_excel(inputfile_path, sheet_name="Controls", dtype=dtypes)
    def _split_unit_columns(dilutions: pd.DataFrame) -> pd.DataFrame:
        """
        Detect columns with units embedded in the header (e.g. "Concentration (mg/mL)"
        or "Concentration in mg/mL") and normalize them by:
        - renaming the value column to a sensible base name (e.g. "Concentration")
        - adding a unit column (e.g. "Unit" or "<Base> Unit") with the detected unit
        """
        unit_header_re = re.compile(r"^(?P<name>.+?)\s*\((?P<unit>.+?)\)\s*$")
        in_header_re = re.compile(r"^(?P<name>.+?)\s+in\s+(?P<unit>.+?)\s*$", re.I)

        rename_map: dict[str, str] = {}
        add_unit_cols: list[tuple[str, str]] = []

        for col in list(dilutions.columns):
            col_str = str(col).strip()

            # handle "Name in Unit" first
            m_in = in_header_re.match(col_str)
            if m_in:
                base = m_in.group("name").strip()
                unit = m_in.group("unit").strip()
                if base.lower() == "concentration":
                    rename_map[col] = "Concentration"
                    unit_col_name = "Unit"
                else:
                    new_base = base
                    if new_base in dilutions.columns and new_base != col:
                        rename_map[col] = col  # keep original
                        unit_col_name = f"{col} Unit"
                    else:
                        rename_map[col] = new_base
                        unit_col_name = f"{new_base} Unit"
                add_unit_cols.append((unit_col_name, unit))
                continue

            # handle "Name (Unit)"
            m = unit_header_re.match(col_str)
            if m:
                base = m.group("name").strip()
                unit = m.group("unit").strip()
                if base.lower() == "concentration":
                    rename_map[col] = "Concentration"
                    unit_col_name = "Unit"
                else:
                    new_base = base
                    if new_base in dilutions.columns and new_base != col:
                        rename_map[col] = col  # keep original
                        unit_col_name = f"{col} Unit"
                    else:
                        rename_map[col] = new_base
                        unit_col_name = f"{new_base} Unit"
                add_unit_cols.append((unit_col_name, unit))

        if rename_map:
            dilutions = dilutions.rename(columns=rename_map)
            for (unit_col_name, unit) in add_unit_cols:
                if unit_col_name not in dilutions.columns:
                    dilutions[unit_col_name] = unit

        return dilutions

    # normalize dilutions column headers and extract unit columns if present
    dilutions = _split_unit_columns(dilutions)


    # Allow endings like 'Position 96', 'Position 384' etc.
    poscol = controls.columns[controls.columns.str.startswith("Position")][0]
    controls["Row_384"] = controls[poscol].apply(lambda x: position_to_rowcol(x)[0])
    controls["Col_384"] = controls[poscol].apply(lambda x: position_to_rowcol(x)[1])
    controls.drop(columns=poscol, inplace=True)

    return substances, organisms, dilutions, controls
