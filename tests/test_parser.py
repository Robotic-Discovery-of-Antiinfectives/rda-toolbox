import os
import pandas as pd
import pytest
from pathlib import Path

from rda_toolbox.parser import _validate_inputfile_structure, read_inputfile, read_platemapping


def _write_excel(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name, index=False)


def test_missing_sheets(tmp_path):
    path = tmp_path / "input_missing_sheets.xlsx"
    # only Substances sheet, others missing -> should report missing sheets
    subs = pd.DataFrame({"MyID": ["s1", "s2"], "Dataset": ["d1", "d2"]})
    _write_excel(path, {"Substances": subs})
    with pytest.raises(ValueError) as exc:
        _validate_inputfile_structure(str(path), "MyID")
    assert "Missing sheets" in str(exc.value)


def test_missing_substance_id_column(tmp_path):
    path = tmp_path / "input_missing_id.xlsx"
    # create all required sheets but Substances lacks the requested id column
    subs = pd.DataFrame({"Internal ID": ["s1", "s2"], "Dataset": ["d1", "d2"]})
    orgs = pd.DataFrame({"Organism": ["E. coli"]})
    dil = pd.DataFrame({"conc": [1.0]})
    ctrl = pd.DataFrame({"Position 96": ["A1"]})
    _write_excel(path, {"Substances": subs, "Organisms": orgs, "Dilutions": dil, "Controls": ctrl})
    with pytest.raises(ValueError) as exc:
        _validate_inputfile_structure(str(path), "MyID")
    assert "Substances sheet does not contain the required column" in str(exc.value)


def test_missing_concentration_unit(tmp_path):
    path = tmp_path / "input_no_unit.xlsx"
    # minimal valid sheets but Dilutions has numeric concentration without unit and no Unit column
    subs = pd.DataFrame({"MyID": ["s1"], "Dataset": ["d1"]})
    orgs = pd.DataFrame({"Organism": ["E. coli"]})
    # concentration column present but header has no unit and no separate Unit column
    dil = pd.DataFrame({"MyID": ["s1"], "concentration": [1.0]})
    ctrl = pd.DataFrame({"Position 96": ["A1"]})
    _write_excel(path, {"Substances": subs, "Organisms": orgs, "Dilutions": dil, "Controls": ctrl})
    with pytest.raises(ValueError) as exc:
        _validate_inputfile_structure(str(path), "MyID")
    assert "no concentration unit" in str(exc.value).lower()


def test_invalid_control_positions(tmp_path):
    path = tmp_path / "input_bad_controls.xlsx"
    subs = pd.DataFrame({"MyID": ["s1"], "Dataset": ["d1"]})
    orgs = pd.DataFrame({"Organism": ["E. coli"]})
    dil = pd.DataFrame({"MyID": ["s1"], "conc": [1.0]})
    # invalid positions (not letter+digits)
    ctrl = pd.DataFrame({"Position 96": ["1A", "X", ""]})
    _write_excel(path, {"Substances": subs, "Organisms": orgs, "Dilutions": dil, "Controls": ctrl})
    with pytest.raises(ValueError) as exc:
        _validate_inputfile_structure(str(path), "MyID")
    assert "Controls.Position" in str(exc.value) or "Controls" in str(exc.value)


def test_read_inputfile_success(tmp_path):
    path = tmp_path / "input_good.xlsx"
    # complete minimal valid input with explicit Unit column
    subs = pd.DataFrame({"MyID": ["s1", "s2"], "Name": ["alpha", "beta"], "Dataset": ["d1", "d2"]})
    orgs = pd.DataFrame({"Organism": ["E. coli", "S. aureus"]})
    dil = pd.DataFrame({"Concentration": ["10", "10"], "Unit": ["uM", "uM"]})
    ctrl = pd.DataFrame({"Position 96": ["A1", "B2"]})
    _write_excel(path, {"Substances": subs, "Organisms": orgs, "Dilutions": dil, "Controls": ctrl})

    substances, organisms, dilutions, controls = read_inputfile(str(path), "MyID")

    # substances should have been renamed to Internal ID
    assert "Internal ID" in substances.columns
    assert substances["Internal ID"].astype(str).tolist() == ["s1", "s2"]

    # organisms should have Organism and formatted column
    assert "Organism" in organisms.columns
    assert "Organism formatted" in organisms.columns

    # controls should now have Row_384 and Col_384 columns, and original Position column removed
    assert "Row_384" in controls.columns
    assert "Col_384" in controls.columns
    assert controls["Row_384"].iloc[0] == "A"
    assert controls["Col_384"].iloc[0] == 1

    # dilutions must contain a Unit column and not be empty
    assert "Unit" in dilutions.columns
    assert dilutions["Unit"].astype(str).tolist() == ["uM", "uM"]


def test_read_inputfile_unit_in_header(tmp_path):
    path = tmp_path / "input_unit_in_header.xlsx"
    # unit provided in the concentration header instead of a separate Unit column
    subs = pd.DataFrame({"MyID": ["s1"], "Dataset": ["d1"]})
    orgs = pd.DataFrame({"Organism": ["E. coli"]})
    dil = pd.DataFrame({"Concentration in mg/mL": ["5.0"]})
    ctrl = pd.DataFrame({"Position 96": ["A1"]})
    _write_excel(path, {"Substances": subs, "Organisms": orgs, "Dilutions": dil, "Controls": ctrl})

    _, _, dilutions, _ = read_inputfile(str(path), "MyID")

    # unit must have been extracted into a Unit column
    assert "Unit" in dilutions.columns
    assert dilutions["Unit"].astype(str).iloc[0] == "mg/mL"

    # there should still be a concentration column (either original or normalized)
    assert any("Concentration" in col for col in dilutions.columns)


def test_read_platemapping_parses_lines_and_counts_replicates():
    filecontents_a = ["5500271031", "011AsT01005;011AsT01006", "5500271056", "007AsT02001"]
    orig_barcodes_a = ["5500271031", "5500271056"]
    assert read_platemapping(filecontents_a, orig_barcodes_a) == (
        {
            "5500271031": [["011AsT01005", "011AsT01006"]],
            "5500271056": [["007AsT02001"]],
        },
        {"5500271031": 1, "5500271056": 1},
    )

    filecontents_b = [
        "011AsT01001",
        "011AcD01064;011AcD01022;011AcD01043;011AcD01085;011AcD01001",
        "007AsT02001",
        "011AcD01065;011AcD01023;011AcD01044;011AcD01086;011AcD01002",
        "011AsT01003",
        "011AcD01066;011AcD01024;011AcD01045;011AcD01087;011AcD01003",
        "011AsT01004",
        "011AcD01067;011AcD01025;011AcD01046;011AcD01088;011AcD01004",
        "011AsT01005",
        "011AcD01068;011AcD01026;011AcD01047;011AcD01089;011AcD01005",
        "011AsT01006",
        "011AcD01069;011AcD01027;011AcD01048;011AcD01090;011AcD01006",
        "011AsT01007",
        "011AcD01070;011AcD01028;011AcD01049;011AcD01091;011AcD01007",
        "011AsT01008",
        "011AcD01071;011AcD01029;011AcD01050;011AcD01092;011AcD01008",
        "011AsT01009",
        "011AcD01072;011AcD01030;011AcD01051;011AcD01093;011AcD01009",
        "011AsT01010",
        "011AcD01073;011AcD01031;011AcD01052;011AcD01094;011AcD01010",
        "002AsT01010",
        "011AcD01074;011AcD01032;011AcD01053;011AcD01095;011AcD01011",
        "002AsT01011",
        "011AcD01075;011AcD01033;011AcD01054;011AcD01096;011AcD01012",
        "002AsT01012",
        "011AcD01076;011AcD01034;011AcD01055;011AcD01097;011AcD01013",
        "002AsT01013",
        "011AcD01077;011AcD01035;011AcD01056;011AcD01098;011AcD01014",
        "002AsT01014",
        "011AcD01078;011AcD01036;011AcD01057;011AcD01099;011AcD01015",
        "002AsT01015",
        "011AcD01079;011AcD01037;011AcD01058;011AcD01100;011AcD01016",
        "002AsT01016",
        "011AcD01080;011AcD01038;011AcD01059;011AcD01101;011AcD01017",
        "002AsT01017",
        "011AcD01081;011AcD01039;011AcD01060;011AcD01102;011AcD01018",
        "002AsT01018",
        "011AcD01082;011AcD01040;011AcD01061;011AcD01103;011AcD01019",
        "002AsT01019",
        "011AcD01083;011AcD01041;011AcD01062;011AcD01104;011AcD01020",
        "002AsT01020",
        "011AcD01084;011AcD01042;011AcD01063;011AcD01105;011AcD01021",
    ]
    orig_barcodes_b = ["011AsT01005", "011AsT01006", "007AsT02001"]
    assert read_platemapping(filecontents_b, orig_barcodes_b) == (
        {
            "011AsT01005": [["011AcD01068", "011AcD01026", "011AcD01047", "011AcD01089", "011AcD01005"]],
            "011AsT01006": [["011AcD01069", "011AcD01027", "011AcD01048", "011AcD01090", "011AcD01006"]],
            "007AsT02001": [["011AcD01065", "011AcD01023", "011AcD01044", "011AcD01086", "011AcD01002"]],
        },
        {"007AsT02001": 1, "011AsT01005": 1, "011AsT01006": 1},
    )

    filecontents_d = [
        "001AsT03001",
        "001AcD03001;001AcD03019;001AcD03037;001AcD03005;001AcD03073",
        "001AsT03001",
        "001AcD03002;001AcD03020;001AcD03038;001AcD03056;001AcD03074",
    ]
    orig_barcodes_d = ["001AsT03001"]
    assert read_platemapping(filecontents_d, orig_barcodes_d) == (
        {
            "001AsT03001": [
                ["001AcD03001", "001AcD03019", "001AcD03037", "001AcD03005", "001AcD03073"],
                ["001AcD03002", "001AcD03020", "001AcD03038", "001AcD03056", "001AcD03074"],
            ],
        },
        {"001AsT03001": 2},
    )
