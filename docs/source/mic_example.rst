Example of a MIC Assay Analysis
=============

Collecting the inputs:
*******************

The first step in the analysis of the MIC assay is to collect all variable inputs
.. Link MIC_Input.xlsx template in the repo here
MIC_Input.xlsx consists of multiple sheets

**Substances**

+-----------------------+-------------+-------------+-------------+--------------------+---------------+----------------+
| Dataset               | External ID | Internal ID | Origin Rack | Origin Position 96 | MP Barcode 96 | MP Position 96 |
+=======================+=============+=============+=============+====================+===============+================+
| Cooperation Partner A | 12345       | 1           |             |                    |               |                |
+-----------------------+-------------+-------------+-------------+--------------------+---------------+----------------+

**Organisms**

+-------------------+------+
| Organism          | Rack |
+===================+======+
| Organism A ST1234 | 1    |
+-------------------+------+

**Dilutions**

+---------------+
| Concentration |
+===============+
| 50            |
| 25            |
| 12.5          |
| 6.25          |
| 3,125         |
| 1.56          |
| 0.78          |
| 0.38          |
| 0.2           |
| 0.1           |
| 0.05          |
+---------------+

**Controls**

+------------------+-------------------+----------+
| Dataset          | Internal ID       | Position |
+==================+===================+==========+
| Negative Control | Bacteria + Medium | A23      |
| Negative Control | Bacteria + Medium | C23      |
| Negative Control | Bacteria + Medium | E23      |
| Negative Control | Bacteria + Medium | G23      |
| Negative Control | Bacteria + Medium | I23      |
| Negative Control | Bacteria + Medium | K23      |
| Negative Control | Bacteria + Medium | M23      |
| Negative Control | Bacteria + Medium | O23      |
| Blank            | Medium            | A24      |
| Blank            | Medium            | B24      |
| Blank            | Medium            | C24      |
| Blank            | Medium            | D24      |
| Blank            | Medium            | E24      |
| Blank            | Medium            | F24      |
| Blank            | Medium            | G24      |
| Blank            | Medium            | H24      |
| ...              | ...               | ...      |
| Reference        | Rifampicin        | N23      |
| Reference        | Vancomycin        | P23      |
+------------------+-------------------+----------+

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
