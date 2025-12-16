import os
import pandas as pd
import pytest
from pathlib import Path

from rda_toolbox.parser import _validate_inputfile_structure, read_inputfile


def _write_excel(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name, index=False)


def test_missing_sheets(tmp_path):
    path = tmp_path / "input_missing_sheets.xlsx"
    # only Substances sheet, others missing -> should report missing sheets
    subs = pd.DataFrame({"MyID": ["s1", "s2"]})
    _write_excel(path, {"Substances": subs})
    with pytest.raises(ValueError) as exc:
        _validate_inputfile_structure(str(path), "MyID")
    assert "Missing sheets" in str(exc.value)


def test_missing_substance_id_column(tmp_path):
    path = tmp_path / "input_missing_id.xlsx"
    # create all required sheets but Substances lacks the requested id column
    subs = pd.DataFrame({"Internal ID": ["s1", "s2"]})
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
    subs = pd.DataFrame({"MyID": ["s1"]})
    orgs = pd.DataFrame({"Organism": ["E. coli"]})
    # concentration column present but header has no unit and no separate Unit column
    dil = pd.DataFrame({"MyID": ["s1"], "concentration": [1.0]})
    ctrl = pd.DataFrame({"Position 96": ["A1"]})
    _write_excel(path, {"Substances": subs, "Organisms": orgs, "Dilutions": dil, "Controls": ctrl})
    with pytest.raises(ValueError) as exc:
        _validate_inputfile_structure(str(path), "MyID")


def test_invalid_control_positions(tmp_path):
    path = tmp_path / "input_bad_controls.xlsx"
    subs = pd.DataFrame({"MyID": ["s1"]})
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
    subs = pd.DataFrame({"MyID": ["s1", "s2"], "Name": ["alpha", "beta"]})
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
    subs = pd.DataFrame({"MyID": ["s1"]})
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