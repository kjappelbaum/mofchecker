# -*- coding: utf-8 -*-
"""Test for the command line interface"""
import json
import pathlib

from click.testing import CliRunner

from mofchecker import cli

from .conftest import THIS_DIR

TEST_DIR = pathlib.Path(THIS_DIR) / "test_files"


def test_json_output():
    """Test that output is valid JSON"""
    runner = CliRunner()
    result = runner.invoke(
        cli.run,
        [str(TEST_DIR / "ABAVIJ_clean.cif"), str(TEST_DIR / "ABAVIJ_clean.cif")],
    )
    assert result.exit_code == 0

    json_list = json.loads(result.output)
    assert len(json_list) == 2, json_list
    assert json_list[0]["name"] == "ABAVIJ_clean", json_list[0]


def test_select_descriptors():
    """Test that one can select specific descriptors"""
    runner = CliRunner()
    result = runner.invoke(
        cli.run,
        [str(TEST_DIR / "ABAVIJ_clean.cif"), "-d", "has_metal"],
    )
    assert result.exit_code == 0

    json_list = json.loads(result.output)
    assert len(json_list) == 1, json_list
    assert list(json_list[0].keys()) == ["has_metal"]
