import pytest

from rda_toolbox.experiment_classes import MIC


def test_validate_mapping_dicts_accepts_consistent_mapping():
    mic = MIC.__new__(MIC)
    mp_ast_mapping_dict = {"MP-1": ["AST-1", "AST-2"]}
    ast_acd_mapping_dict = {
        "AST-1": ["ACD-1", "ACD-2"],
        "AST-2": ["ACD-3", "ACD-4"],
    }

    mic._validate_mapping_dicts(mp_ast_mapping_dict, ast_acd_mapping_dict)


def test_validate_mapping_dicts_raises_for_inconsistent_mapping():
    mic = MIC.__new__(MIC)
    mp_ast_mapping_dict = {"MP-1": ["AST-1", "AST-MISSING", ""]}
    ast_acd_mapping_dict = {
        "AST-1": ["ACD-1", ""],
        "": ["ACD-2"],
    }

    with pytest.raises(ValueError) as exc:
        mic._validate_mapping_dicts(mp_ast_mapping_dict, ast_acd_mapping_dict)

    message = str(exc.value)
    assert "Please check the mapping .txt files." in message
    assert "AsT barcodes missing in AsT -> AcD mapping" in message
    assert "Invalid AsT barcodes in MP -> AsT mapping" in message
    assert "Invalid AsT barcodes in AsT -> AcD mapping" in message
    assert "Invalid AcD barcodes in AsT -> AcD mapping" in message


def test_validate_mapping_dicts_accepts_switched_order():
    mic = MIC.__new__(MIC)
    mp_ast_mapping_dict = {"MP-1": ["AST-2", "AST-1"]}
    ast_acd_mapping_dict = {
        "AST-1": ["ACD-2", "ACD-1"],
        "AST-2": ["ACD-4", "ACD-3"],
    }

    mic._validate_mapping_dicts(mp_ast_mapping_dict, ast_acd_mapping_dict)


def test_validate_mapping_dicts_raises_for_badly_switched_entries():
    mic = MIC.__new__(MIC)
    mp_ast_mapping_dict = {"MP-1": ["AST-1", "AST-2"]}
    # AST-2 got "switched" to AST-22 by mistake -> AST-2 is now missing
    ast_acd_mapping_dict = {
        "AST-1": ["ACD-1", "ACD-2"],
        "AST-22": ["ACD-3", "ACD-4"],
    }

    with pytest.raises(ValueError) as exc:
        mic._validate_mapping_dicts(mp_ast_mapping_dict, ast_acd_mapping_dict)

    message = str(exc.value)
    assert "Please check the mapping .txt files." in message
    assert "AsT barcodes missing in AsT -> AcD mapping" in message
    assert "AST-2" in message
