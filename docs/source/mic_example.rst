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
| Organism B ST5678 | 2    |
+-------------------+------+
| Organism C ST9101 | 3    |
+-------------------+------+
| ...               | ...  |
+-------------------+------+

**Dilutions**

+---------------+
| Concentration |
+===============+
| 50            |
+---------------+
| 25            |
+---------------+
| 12.5          |
+---------------+
| 6.25          |
+---------------+
| 3,125         |
+---------------+
| 1.56          |
+---------------+
| 0.78          |
+---------------+
| 0.38          |
+---------------+
| 0.2           |
+---------------+
| 0.1           |
+---------------+
| 0.05          |
+---------------+

**Controls**

+------------------+-------------------+----------+
| Dataset          | Internal ID       | Position |
+==================+===================+==========+
| Negative Control | Bacteria + Medium | A23      |
+------------------+-------------------+----------+
| Negative Control | Bacteria + Medium | C23      |
+------------------+-------------------+----------+
| Negative Control | Bacteria + Medium | E23      |
+------------------+-------------------+----------+
| Negative Control | Bacteria + Medium | G23      |
+------------------+-------------------+----------+
| Negative Control | Bacteria + Medium | I23      |
+------------------+-------------------+----------+
| Negative Control | Bacteria + Medium | K23      |
+------------------+-------------------+----------+
| Negative Control | Bacteria + Medium | M23      |
+------------------+-------------------+----------+
| Negative Control | Bacteria + Medium | O23      |
+------------------+-------------------+----------+
| Blank            | Medium            | A24      |
+------------------+-------------------+----------+
| Blank            | Medium            | B24      |
+------------------+-------------------+----------+
| Blank            | Medium            | C24      |
+------------------+-------------------+----------+
| Blank            | Medium            | D24      |
+------------------+-------------------+----------+
| Blank            | Medium            | E24      |
+------------------+-------------------+----------+
| Blank            | Medium            | F24      |
+------------------+-------------------+----------+
| Blank            | Medium            | G24      |
+------------------+-------------------+----------+
| Blank            | Medium            | H24      |
+------------------+-------------------+----------+
| ...              | ...               | ...      |
+------------------+-------------------+----------+
| Reference        | Rifampicin        | N23      |
+------------------+-------------------+----------+
| Reference        | Vancomycin        | P23      |
+------------------+-------------------+----------+

.. code-block:: Python
   :linenos:

   import rda_toolbox as rda

   input_mapping = rda.mic_process_inputs(
       # Input specifications excel:
       "../data/input/MIC_Input.xlsx",
       # Barcode reader file which shows Motherplate to AsT plate mapping
       "../data/input/DiS_MP_AsT_2024-10-08.txt",
       # Barcode reader file which shows Ast plate to AcD plate mapping
       "../data/input/AmA_AsT_AcD_20241009.txt",
       # Raw files path
       "../data/raw/",
   )
   input_mapping.to_csv("../data/processed/prepared_input_mapping_table.csv", index=False)


+----------------+-------------+---------------+-----------------+-------------+---------------+----------------+---------+---------+-----------------+---------------+----------+-----------------+-----------+----------+
| Dataset        | External ID | Original Rack | Origin Position | Internal ID | MP Barcode 96 | MP Position 96 | Row_384 | Col_384 | AsT Barcode 384 | Concentration | Position | AcD Barcode 384 | Replicate | Organism |
+----------------+-------------+---------------+-----------------+-------------+---------------+----------------+---------+---------+-----------------+---------------+----------+-----------------+-----------+----------+
| Coop Partner A | 123456789   | 2233445566    | A1              | 987654321   | 3456677889    | A1             |         |         |                 |               |          |                 |           |          |
+----------------+-------------+---------------+-----------------+-------------+---------------+----------------+---------+---------+-----------------+---------------+----------+-----------------+-----------+----------+
|                |             |               |                 |             |               |                |         |         |                 |               |          |                 |           |          |
+----------------+-------------+---------------+-----------------+-------------+---------------+----------------+---------+---------+-----------------+---------------+----------+-----------------+-----------+----------+

.. code-block:: Python
   :linenos:

   rawfiles = rda.parse_readerfiles("../data/raw/")
   rawfiles.to_csv("../data/processed/rawdata.csv", index=False)

+---------+---------+---------------------+-----------------+
| Row_384 | Col_384 | Raw Optical Density | AcD Barcode 384 |
+---------+---------+---------------------+-----------------+
|         |         |                     |                 |
+---------+---------+---------------------+-----------------+
|         |         |                     |                 |
+---------+---------+---------------------+-----------------+


.. code-block:: Python
   :linenos:

   preprocessed_data = rda.preprocess(
       input_mapping,
       substance_id="Internal ID",
       measurement="Optical Density",
       negative_controls="Bacteria + Medium",
       blanks="Medium",
       norm_by_barcode="AcD Barcode 384"
   )
   preprocessed_data.to_csv("../data/processed/preprocessed_data.csv", index=False)


First Visualizations for Quality Control:
*******************

After processing all the inputs and preprocessing the rawdata like subtracting background noise etc., we can do some first visualizations.

Create faceted heatmaps for each raw plate for quality control.

.. code-block:: Python
   :linenos:

    plate_heatmaps = rda.plateheatmaps(preprocessed_data, substance_id="Internal ID", barcode="AsT Barcode 384", negative_control="Bacteria + Medium", blank="Medium")
    plate_heatmaps.save("../figures/plateheatmaps.svg")
    plate_heatmaps.save("../figures/plateheatmaps.html")


Create lineplots with a horizontal rule at 50% relative growth.

.. code-block:: Python
   :linenos:

    lineplots = rda.lineplots_facet(
        preprocessed_data
        hline_y=50
    )
    lineplots.save(f"../figures/MIC_Lineplots_AllDatasets.svg")
    lineplots.save(f"../figures/MIC_Lineplots_AllDatasets.html")


Obtain MIC results:
*******************

Calculate the mean between replicates, apply MIC threshold, save the results for each dataset in its corresponding folders under <Projectfolder>/data/results/{Datasets}.

.. code-block:: Python
   :linenos:

   mic_results(preprocessed_data, "../data/results/")
