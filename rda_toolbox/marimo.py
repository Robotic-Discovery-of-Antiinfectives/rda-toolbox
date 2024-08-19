#!/usr/bin/env python3

from io import StringIO
from .parser import readerfile_parser
from .parser import collect_results
from .parser import collect_metadata
# from marimo import FileUploadResults


def readeruploads_rawdf(uploadfiles):
    """
    Wrapper function to collect the results from marimos file upload button.
    This function is similar to readerfiles_rawdf() except that the input
    is a list of marimo.ui.file() outputs instead of a list of filepaths.
    """
    if not uploadfiles:
        return None
    filedicts = [
        readerfile_parser(upload.name, StringIO(str(upload.contents, 'utf-8')))
        for upload in uploadfiles
    ]

    result_df = collect_results(filedicts)
    metadata_df = collect_metadata(filedicts)

    return result_df
