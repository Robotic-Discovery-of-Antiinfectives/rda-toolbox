# Robotic-assisted Discovery of Antiinfectives

This package aims to provide a toolbox for data analysis in the field of drug discovery.

---

The aim is to provide functions to help evaluate the following assays:
- Primary Screen
- MIC (Minimum Inhibitory Concentration) Assay
- Cellviability

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
