# Minimal Inhibitory Concentration (MIC) Assay


## Read the inputs, initialize the assay class

```Python
import rda_toolbox as rda

mic = rda.MIC(
    "../data/raw/",
    "../data/input/MIC_Input.xlsx",
    "../data/input/DiS_MP_AsT_2024-12-02.txt",
    "../data/input/AmA_AsT_AcD_20241204.txt",
    precipitation_rawfilepath = "../data/raw/Precipitation_measurements/",
)
```

Its possible to inspect the assay object:

```Python
mic.__dict__
```

## View in-between results (e.g. in a notebook)

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
mic.save("../") # project root as directory path
```
