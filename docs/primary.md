# Primary Screen Assay



```Python
import rda_toolbox as rda

primary = rda.PrimaryScreen(
    "../data/raw/",  # Folder where the raw readerfiles are located
    "../data/input/PrS_Input.xlsx",  # Assay specific Input sheet
    "../data/input/AmA_AsT_AcD_20241204.txt",  # Mapping file
    map_rowname="Row 96",
    map_colname="Col 96",
    # Folder where the raw readerfiles for precipitation test are located
    precipitation_rawfilepath = "../data/raw/Precipitation_measurements/"
)
```


Its possible to inspect the assay:

```Python
primary.__dict__
```

```Python
# Save all tables
primary.save_tables("../data/results/")
# Save all figures
primary.save_figures("../figures/")
# Save results (figures and tables)
primary.save()
```
