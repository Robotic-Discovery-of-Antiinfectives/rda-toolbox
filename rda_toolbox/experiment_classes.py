#!/usr/bin/env python

import pandas as pd
import altair as alt
from functools import cached_property

import string


from .utility import get_rows_cols, mapapply_96_to_384
from .parser import parse_readerfiles, read_inputfile, parse_mappingfile
from .process import preprocess


class Experiment:
    """
    Superclass for all experiments.
    Reads rawdata into a DataFrame.

    Attributes
    ----------
    rawdata : pd.DataFrame
        DataFrame containing the rawdata

    Methods
    ----------
    save_plots
        Save all the resulting plots to figuredir
    save_tables
        Save all the resulting tables to tabledir
    save
        Save all plots and tables to resultdir
    """

    def __init__(self, rawfiles_folderpath: str, plate_type: int):
        self._plate_type = plate_type
        self._rows, self._columns = get_rows_cols(plate_type)
        self._rawfiles_folderpath = rawfiles_folderpath
        self.rawdata = parse_readerfiles(rawfiles_folderpath)

    # def save_plots(self, figuredir: str):
    # def save_tables(self, tabledir: str):
    # def save(self, resultdir: str):
    #     self.save_plots()
    #     self.save_tables()


class Precipitation(Experiment):
    def __init__(
        self,
        rawfiles_folderpath: str,
        plate_type: int = 384,  # Define default plate_type for experiment
        measurement_label: str = "Raw Optical Density",
        background_locations: type[None | pd.DataFrame] = None,
    ):
        super().__init__(rawfiles_folderpath, plate_type)
        self._measurement_label = measurement_label

        # Set default background location as 24th column
        if background_locations is None:
            self.background_locations = pd.DataFrame(
                {
                    f"Row_{self._plate_type}": list(
                        string.ascii_uppercase[: self._rows]
                    ),
                    f"Col_{self._plate_type}": [
                        self._columns for i in range(self._rows)
                    ],
                    "Layout": "Background",
                }
            )
        else:
            if not "Layout" in background_locations:
                background_locations["Layout"] = "Background"
            self.background_locations = background_locations.rename(
                columns={
                    "Row": f"Row_{self._plate_type}",
                    "Column": f"Col_{self._plate_type}",
                }
            )

        self.rawdata_w_layout = pd.merge(
            self.rawdata, self.background_locations, how="outer"
        ).fillna({"Layout": "Substance"})
        # print(self.get_results())
        # self.results = self.get_results()

    @property
    def limit_of_quantification(self):  # "Bestimmungsmaß"
        background = self.rawdata_w_layout[
            self.rawdata_w_layout["Layout"] == "Background"
        ][self._measurement_label]
        loq = round(background.mean() + 10 * background.std(), 3)
        self.rawdata_w_layout["Limit of Quantification"] = loq
        return loq

    @cached_property
    def results(self):
        self.rawdata_w_layout["Precipitated"] = self.rawdata_w_layout[
            self._measurement_label
        ].apply(lambda x: x > self.limit_of_quantification)
        return self.rawdata_w_layout

    # let it have its own heatmap function for now:
    def plateheatmap(self):
        base = alt.Chart(
            self.results,
        ).encode(
            alt.X("Col_384:O").axis(labelAngle=0, orient="top").title(None),
            alt.Y("Row_384:O").title(None),
            tooltip=list(self.results.columns),
        )

        heatmap = base.mark_rect().encode(
            alt.Color(
                self._measurement_label,
                scale=alt.Scale(
                    scheme="redyellowblue",
                    domain=[0, self.limit_of_quantification, 1],
                    reverse=True,
                ),
            ).title(self._measurement_label)
        )
        text = base.mark_text(baseline="middle", align="center", fontSize=7).encode(
            alt.Text(f"{self._measurement_label}:Q", format=".1f"),
            color=alt.condition(
                alt.datum[self._measurement_label]
                < max(self.results[self._measurement_label]) / 2,
                alt.value("black"),
                alt.value("white"),
            ),
        )
        return alt.layer(heatmap, text).facet(
            column="AcD Barcode 384",
            title=alt.Title(
                "Precipitation Test",
                subtitle=[f"Limit of Quantification: {self.limit_of_quantification}"],
            ),
        )


class PrimaryScreen(Experiment):
    def __init__(
        self,
        rawfiles_folderpath: str,
        inputfile_path: str,
        mappingfile_path: str,
        plate_type: int = 384,  # Define default plate_type for experiment
        measurement_label: str = "Raw Optical Density",
        map_rowname: str = "Row_96",
        map_colname: str = "Col_96",
        q_name: str = "Quadrant",
        substance_id: str = "Internal ID",
        negative_controls: str = "Bacteria + Medium",
        blanks: str = "Medium",
        norm_by_barcode: str = "AcD Barcode 384",
        thresholds : list[float] = [50.0],
    ):
        # TODO: inherit the save_* functions from Experiment superclass
        super().__init__(rawfiles_folderpath, plate_type)
        self._measurement_label = measurement_label
        self._mappingfile_path = mappingfile_path
        self._inputfile_path = inputfile_path
        self._substances_unmapped, self._organisms, self._dilutions, self._controls = (
            read_inputfile(inputfile_path)
        )
        self.substances = mapapply_96_to_384(
            self._substances_unmapped,
            rowname=map_rowname,
            colname=map_colname,
            q_name=q_name,
        )
        self._mapping_df = parse_mappingfile(
            mappingfile_path,
            motherplate_column="AsT Barcode 384",
            childplate_column="AcD Barcode 384",
        )
        self._substance_id = substance_id
        self._negative_controls = negative_controls
        self._blanks = blanks
        self._norm_by_barcode = norm_by_barcode

    @cached_property
    def mapped_input_df(self):
        control_wbarcodes = []
        # multiply controls with number of AsT plates to later merge them with substances df
        for origin_barcode in list(self.substances["AsT Barcode 384"].unique()):
            controls_subdf = self._controls.copy()
            controls_subdf["AsT Barcode 384"] = origin_barcode
            control_wbarcodes.append(controls_subdf)
        controls_n_barcodes = pd.concat(control_wbarcodes)

        ast_plate_df = pd.merge(
            pd.concat([self.substances, controls_n_barcodes]),
            self._dilutions,
            how="outer",
        )

        mapped_organisms = pd.merge(self._mapping_df, self._organisms)

        result_df = pd.concat(
            [
                pd.merge(org_df, ast_plate_df)
                for _, org_df in pd.merge(mapped_organisms, self.rawdata).groupby(
                    "Organism"
                )
            ]
        )

        for ast_barcode, ast_plate in result_df.groupby("AsT Barcode 384"):
            print(
                f"AsT Plate {ast_barcode} has size: {len(ast_plate)//len(ast_plate['AcD Barcode 384'].unique())}"
            )
            print(f"{ast_barcode} -> {ast_plate["AcD Barcode 384"].unique()}")
        return result_df

    @cached_property
    def processed(self):
        return preprocess(
            self.mapped_input_df,
            substance_id = self._substance_id,
            measurement = self._measurement_label.strip(
                "Raw "
            ),  # I know this is weird, its because of how background_normalize_zfactor works,
            negative_controls = self._negative_controls,
            blanks = self._blanks,
            norm_by_barcode = self._norm_by_barcode,
        )

    @cached_property
    def results(self):
        df = self.processed.copy()
        df = df[
            (df["Dataset"] != "Reference")
            & (df["Dataset"] != "Positive Control")
            & (df["Dataset"] != "Blank")
        ].dropna(subset=["Concentration"])

        pivot_df = pd.pivot_table(
            df,
            values=["Relative Optical Density", "Replicate", "Z-Factor"],
            index=[
                "Internal ID",
                "Organism",
                "Concentration",
                "Dataset",
            ],
            aggfunc={
                "Relative Optical Density": ["mean"],
                "Replicate": ["count"],
            },
        ).reset_index()
        pivot_df.columns = [" ".join(x).strip() for x in pivot_df.columns.ravel()]

    # def plateheatmap(self):


class MIC(Experiment):  # Minimum Inhibitory Concentration
    def __init__(
        self,
        rawfiles_folderpath,
        plate_type=384,  # Define default plate_type for experiment
        measurement_label="Raw Optical Density",
    ):
        super().__init__(rawfiles_folderpath, plate_type)
        self._measurement_label = measurement_label

    @cached_property
    def results(self):
        pass
