# Setting Up An Experiment

In this example I will use marimo notebooks. `rda-toolbox` should work with any other notebook or virtual environment utility as well as in plain python files.

## Initialize a Virtual Environment

```Sh
uv venv
```

```Sh
source .venv/bin/activate
```

## Marimo Notebooks

Install marimo:
```Sh
uv pip install marimo
```

Start a notebook using a sandboxed environment:
```Sh
marimo edit --sandbox Analysis.py
```
