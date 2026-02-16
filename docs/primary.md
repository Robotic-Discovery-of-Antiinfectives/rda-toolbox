# Primary Screen Assay



## Read the inputs, initialize the assay class

```Python
import rda_toolbox as rda

primary = rda.PrimaryScreen(
    "../data/raw/",  # Folder where the raw readerfiles are located
    "../data/input/PrS_Input.xlsx",  # Assay specific Input sheet
    "../data/input/AmA_AsT_AcD_20241204.txt",  # Mapping file (MP -> AcD plates)
    map_rowname="Row 96",
    map_colname="Col 96",
    # Folder where the raw readerfiles for precipitation test are located
    precipitation_rawfilepath = "../data/raw/Precipitation_measurements/"
)
```

### Compounds already in 384-well plates
Example without mapping from 96 to 384 (motherplate barcode is reused as AsT Barcode):

```Python
import rda_toolbox as rda

primary = rda.PrimaryScreen(
    "../data/raw/",
    "../data/input/PrS_Input.xlsx",
    mappingfile_path = "../data/input/AmA_AsT_AcD_20251202.txt",  # Contains mapping of AsT -> AcD plates
    needs_mapping=False,
    precipitation_rawfilepath="../data/raw/Precipitation/",
    substance_id="Internal ID",
    ast_barcode_header="MP Barcode 384"  # Motherplate is reused as AsT Barcode
)
```

Its possible to inspect the assay object:

```Python
primary.__dict__
```

## View in-between results (e.g. in a notebook)

### Tables

```Python
primary.mapped_input_df
```

```Python
primary.processed
```

```Python
primary.results
```

### Visualizations

```Python
primary.plateheatmap
```

## Save the results

```Python
# Save all tables
primary.save_tables("../data/results/")
# Save all figures
primary.save_figures("../figures/")
# Save results (figures and tables)
primary.save_results(<tables path>, <figures path>, <processed data path>, figureformats=["svg", "html"], tableformats=["xlsx", "csv"])
```

# Primary Screen Results
