#!/usr/bin/env python3

import os

import pytest

from rda_toolbox.parser import readerfiles_to_df

testfile = os.path.abspath("./tests/rawfiles_testfolder/test.txt")
# os.path.abspath(rawfiles_testfolder)
# rawfile_paths = os.listdir(rawfiles_testfolder)


# def validator():
#     return
# def test_filepaths_to_filedicts():
#     assert readerfiles_to_df(rawfile_paths) is None

# readerfiles_to_df([testfile])
