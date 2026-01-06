# TODO: Write tests for every function in rda_toolbox/utility.py

# There are zero automated checks that mapapply_96_to_384,
# background_normalize_zfactor, Altair plot builders,
# or MIC result formatting behave as expected.
# Adding lightweight tests
# (e.g., fixtures that feed synthetic plate data through background_normalize_zfactor
#  and verify normalization/Z‑factor math, or verifying plateheatmaps leaves input untouched)
# will catch the stability issues above before they regress again.

# TODO: write a test for rda_toolbox/utility.py:mapapply_96_to_384
# to verify that it correctly maps 96-well data into 384-well format.

# TODO: write a test for rda_toolbox/utility.py:background_normalize_zfactor
# to verify that it correctly normalizes plate data and computes Z-factors.
