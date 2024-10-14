Example of a MIC Assay Analysis
=============

Collecting the inputs:
*******************

The first step in the analysis of the MIC assay is to collect all variable inputs
.. Link MIC_Input.xlsx template in the repo here


.. code-block:: Python
   :linenos:

    import rda_toolbox as rda

    rda.mic_process_inputs(
        # Input specifications excel:
        "../data/input/MIC_Input.xlsx",
        # Barcode reader file which shows Motherplate to AsT plate mapping
        "../data/input/DiS_MP_AsT_2024-10-08.txt",
        # Barcode reader file which shows Ast plate to AcD plate mapping
        "../data/input/AmA_AsT_AcD_20241009.txt",
    ).save("../data/processed/prepared_input_mapping_table.csv", index=False)
