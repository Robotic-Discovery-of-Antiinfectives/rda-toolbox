# Robotic-assisted Discovery of Antiinfectives

This package aims to provide a toolbox for data analysis in the field of drug discovery.

- **[Docs](https://robotic-discovery-of-antiinfectives.github.io/rda-toolbox/)**
- **[PyPi](https://pypi.org/project/rda-toolbox/)**


---

The aim is to provide functions to help evaluate the following assays:
- Primary Screen
- MIC (Minimum Inhibitory Concentration) Assay
- Cellviability

### Usage Example
`pip install rda-toolbox`

or

`pip install "git+https://github.com/Robotic-Discovery-of-Antiinfectives/rda-toolbox.git"`


```Python
#!/usr/bin/env python3

import rda_toolbox as rda
import glob

rda.readerfiles_rawdf(glob.glob("path/to/raw/readerfiles/*"))
```


### File Parsing
- Read output files and return readouts in a [tidy](https://r4ds.had.co.nz/tidy-data.html), [long](https://towardsdatascience.com/long-and-wide-formats-in-data-explained-e48d7c9a06cb) DataFrame

#### **Supported readers:**
- Cytation C10

### Plotting
This package uses [Vega-Altair](https://altair-viz.github.io/index.html) for creating interactive (or static) visualizations.

- Plate Heatmaps
- Upset plots
  - `UpSetAltair` plotting function is taken from https://github.com/hms-dbmi/upset-altair-notebook and modified
  - This part of this package is licensed under MIT license.
<!-- https://testdriven.io/blog/python-project-workflow/ -->


### New Release
1) Update `pyproject.toml` release version
2) Update `docs/source/conf.py` release version
3) On GitHub go to *releases* and `Draft a new release`

### This package is managed via [Poetry](https://python-poetry.org/docs/basic-usage/)
- `poetry build`
- `poetry publish`
