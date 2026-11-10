"""Unit tests for dijetfit.common: category selection, labels, file-name parsing, message counting."""
import pytest

from dijetfit.common import category_labels, count_messages, mass_from_filename, select_categories

NAMES = ["Jet_dY0p0to0p5_mass", "Jet_dY0p5to1p0_mass", "Jet_dY0p0to1p0_mass", "Jet_dY0p0to0p5_pt",
         "eventsProcessed"]


def test_select_with_capture_group_and_exclude():
    cats, used = select_categories(NAMES, [r"(Jet_dY\d+p\d+to\d+p\d+)_mass"], [r"Jet_dY0p0to1p0_mass"])
    assert cats == {"Jet_dY0p0to0p5": "Jet_dY0p0to0p5_mass", "Jet_dY0p5to1p0": "Jet_dY0p5to1p0_mass"}
    assert used == {r"(Jet_dY\d+p\d+to\d+p\d+)_mass": 2}


def test_select_without_group_uses_full_name_and_whole_match():
    cats, _ = select_categories(NAMES, [r"Jet_dY0p0to0p5_.*"])
    assert cats == {"Jet_dY0p0to0p5_mass": "Jet_dY0p0to0p5_mass", "Jet_dY0p0to0p5_pt": "Jet_dY0p0to0p5_pt"}
    # a partial match is not enough
    cats, used = select_categories(NAMES, [r"Jet_dY0p0"])
    assert cats == {} and used == {r"Jet_dY0p0": 0}


def test_select_several_patterns_and_unused_pattern():
    cats, used = select_categories(NAMES, [r"(Jet_dY0p0to0p5)_mass", r"(Jet_dY0p5to1p0)_mass", "nomatch"])
    assert sorted(cats) == ["Jet_dY0p0to0p5", "Jet_dY0p5to1p0"]
    assert used["nomatch"] == 0


def test_select_duplicate_category_names_is_an_error():
    with pytest.raises(ValueError):
        select_categories(NAMES, [r"(Jet_dY0p0to0p5)_.*"])


def test_category_labels_strip_common_prefix_at_word_boundary():
    lab = category_labels(["Ana_Sel_dY0p0to0p5", "Ana_Sel_dY0p5to1p0"])
    assert lab == {"Ana_Sel_dY0p0to0p5": "dY0p0to0p5", "Ana_Sel_dY0p5to1p0": "dY0p5to1p0"}
    assert category_labels(["only_one"]) == {"only_one": "only_one"}


def test_mass_from_filename():
    assert mass_from_filename("/a/b/ZprimeToQQ_M500_TuneCP5.root") == 500
    assert mass_from_filename("x-ZprimeQQ_M1250_synthetic.root") == 1250
    assert mass_from_filename("data.root") is None


def test_count_messages():
    text = "\n".join([
        "[#0] ERROR:Minimization -- RooMinimizer::calculateHessErrors() Error when calculating Hessian",
        "[#1] WARNING:InputArguments -- something",
        "Warning in <TCanvas::Constructor>: Deleting canvas with same name",
        "Error in <TFormula::Eval>: bad",
        "[#1] INFO:Fitting -- RooAddition::defaultErrorLevel(...) using its error level",
        "plain line with errors and warnings in lower case",
    ])
    assert count_messages(text) == (2, 2)
