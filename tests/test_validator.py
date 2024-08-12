#!/usr/bin/env python3

import os

import pytest

from rda_toolbox.parser import filepaths_to_df

rawfiles_testfolder = "./tests/rawfiles_testfolder"
rawfile_paths = os.listdir(rawfiles_testfolder)


# def validator():
#     return
def test_filepaths_to_filedicts():
    assert filepaths_to_df(rawfile_paths) is None
