# Minimal Inhibitory Concentration (MIC) Assay


### Read the inputs, initialize the assay class

```Python
import rda_toolbox as rda

mic = rda.MIC(
    "../data/raw/",  # Folderpath for rawfiles
    "../data/input/MIC_Input.xlsx",  # Input excel table
    "../data/input/DiS_MP_AsT_2024-12-02.txt",  # Mapping file from Motherplates to AssayTransfer plates
    "../data/input/AmA_AsT_AcD_20241204.txt",  # Mapping file from AssayTransfer to ActivityDetermination plates
    plate_type=384,
    measurement_label="Raw Optical Density",
    negative_controls="Organism + Medium",  # Label of negative controls ('Bacteria + Medium', 'Fungi + Medium', 'Organism + Medium', 'Negative Controls' etc.)
    precip_exclude_outlier=True,  # Exclude outliers from the precipitation
    precipitation_rawfilepath = "../data/raw/Precipitation_measurements/",  # Folderpath for precipitation rawfiles
)
```

### Cytation 10 readout table header
In the Cytation C10 reader software you can define things like table headers.
To detect the result matrix and be flexible in naming this table, the keyword `cyt10_matrixheader_mapping` was introduced.
The default value for this keyword is `cyt10_matrixheader_mapping = {"Results": "Raw Optical Density"}` (Usually its better to leave out the "Raw").

```Python
import rda_toolbox as rda

mic = rda.MIC(
    rawfiles_folderpath="../data/raw/24 h/",
    inputfile_path="../data/input/MIC_Input.xlsx",
    mp_ast_mapping_filepath="../data/input/DiS_MP_AsT_2024-12-02.txt",
    ast_acd_mapping_filepath="../data/input/AmA_AsT_AcD_20260317.txt",
    precipitation_rawfilepath="../data/raw/Precipitation/",
    cyt10_matrixheader_mapping={"Read 1:554,593": "Fluorescence", "Read 2:450": "Optical Density"},
    negative_controls="Organism + Medium",
)
```

### Save the results

```Python
mic.save_results(<tables path>, <figures path>, <processed data path>, figureformats=["svg", "html"], tableformats=["xlsx", "csv"])
```

If everything went well you can stop now.

If errors occured you may inspect in-between results and debug from there.

(**Check your inputs!**)

---

## In-between inspection:

### View in-between results (e.g. in a notebook)

Its possible to inspect the assay object:

```Python
mic.__dict__
```

Show how plates are related to each other (hirarchical dictionary):
```Python
mic._mapping_dict
```

### Tables

```Python
mic.mapped_input_df
```

```Python
mic.processed
```

```Python
mic.results
```

### Visualizations

```Python
mic.plateheatmap
```

### Save the results separately

```Python
mic.save_tables("../data/results/")
```

```Python
mic.save_figures("../figures/")
```
