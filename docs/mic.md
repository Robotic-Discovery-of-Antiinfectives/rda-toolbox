# Minimal Inhibitory Concentration (MIC) Assay


## Read the inputs, initialize the assay class

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


## View in-between results (e.g. in a notebook)

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

## Save the results

```Python
mic.save_tables("../data/results/")
```

```Python
mic.save_figures("../figures/")
```

```Python
mic.save_results("../") # project root as directory path
```
