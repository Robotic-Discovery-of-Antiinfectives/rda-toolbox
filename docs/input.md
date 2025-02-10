# Input Specifications

Some functions tightly adhere to these specifications and require a certain setup to work.<br>
Most experiments require a specific input excel sheet containing information about assay parameters.

---

Because reproducibility is a **very** important topic, especially in science:

- It is highly recommended to use [uv](https://docs.astral.sh/uv/pip/environments/) for environments
- Create a virtual environment for each project (Screen)
- Always at least provide a requirements.txt (`uv pip freeze > requirements.txt`) for each finished project


## Example Project Folder Structure

```Sh
<YYYYMMDD_Assay_OrganismType>
├── code/
│   └── .venv/ # uv venv .venv
├── data/
│   ├── input/
│   │   └── <Assay>_Input.xlsx
│   ├── meta/
│   │   └── logs/
│   ├── processed/
│   ├── raw/
│   └── results/
├── figures/
├── methods/
│   └── module_01/
├── readme.md
└── report/
```

- [MIC Input File](https://github.com/Robotic-Discovery-of-Antiinfectives/rda-toolbox/blob/main/Projectfolder_TEMPLATE/data/input/MIC_Input.xlsx)
- [Primary Input File]()
- []()
