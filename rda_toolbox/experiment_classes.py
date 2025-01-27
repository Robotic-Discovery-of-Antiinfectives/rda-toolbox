#!/usr/bin/env python

import pandas as pd
import altair as alt
from functools import cached_property
from dataclasses import dataclass

import os
import pathlib

import string


from .utility import (
    get_rows_cols,
    mapapply_96_to_384,
    get_upsetplot_df,
    get_mapping_dict,
    add_precipitation,
)
from .parser import parse_readerfiles, read_inputfile, parse_mappingfile
from .process import preprocess, get_thresholded_subset
from .plot import plateheatmaps, UpSetAltair


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

    # def save_plots(self, figuredir: str, result_plots: list[Resultfigure]):
    #     pass

    def save_tables(self, resultpath: str):
        """
        Save result tables to "<resultpath>".
        Creates corresponding folders for each dataset.
        """
        for result in self._resulttables:  # cached property of subclasses
            filedir = os.path.join(resultpath, result.dataset)
            pathlib.Path(filedir).mkdir(parents=True, exist_ok=True)
            result.table.to_excel(
                os.path.join(filedir, f"{result.table_basename}.xlsx"), index=False
            )
            result.table.to_csv(
                os.path.join(filedir, f"{result.table_basename}.csv"), index=False
            )

    def save_plots(self, resultpath: str, fileformats: list[str] = ["xlsx", "csv"]):
        for result in self._resultfigures:  # cached property of subclasses
            filedir = os.path.join(resultpath, result.dataset)
            pathlib.Path(filedir).mkdir(parents=True, exist_ok=True)
            for file_format in fileformats:
                result.figure.save(
                    os.path.join(filedir, f"{result.table_basename}.{file_format}"),
                    index=False,
                )

    # def save(self, resultpath: str):
    #     self.save_plots()
    #     self.save_tables()


@dataclass
class Result:
    dataset: str
    file_basename: str
    table: pd.DataFrame | None = None
    figure: alt.Chart | None = None


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
            if "Layout" not in background_locations:
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
        if (
            len(self.results["AcD Barcode 384"].unique()) % 2 == 0
        ):  # even amount of AcD Barcodes
            col_num = 4
        else:  # uneven amount of AcD Barcodes
            col_num = 3
        return (
            alt.layer(heatmap, text)
            .facet(
                facet="AcD Barcode 384",
                title=alt.Title(
                    "Precipitation Test",
                    subtitle=[
                        f"Limit of Quantification: {self.limit_of_quantification}"
                    ],
                ),
                columns=col_num,
            )
            .resolve_axis(x="independent", y="independent")
        )


class PrimaryScreen(Experiment):
    """
    Primary screen experiment. Usually done using only 1 concentration.
    """

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
        thresholds: list[float] = [50.0],
        precipitation_rawfilepath: str | None = None,
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
        self._mapping_dict = get_mapping_dict(self._mapping_df)
        self._substance_id = substance_id
        self._negative_controls = negative_controls
        self._blanks = blanks
        self._norm_by_barcode = norm_by_barcode
        self.thresholds = thresholds
        self.precipitation = (
            Precipitation(precipitation_rawfilepath)
            if precipitation_rawfilepath
            else None
        )
        self.rawdata = (  # Overwrite rawdata if precipitation data is available
            self.rawdata
            if self.precipitation is None
            else add_precipitation(
                self.rawdata, self.precipitation.results, self._mapping_dict
            )
        )
        self._processed_only_substances = self.processed[
            (self.processed["Dataset"] != "Reference")
            & (self.processed["Dataset"] != "Positive Control")
            & (self.processed["Dataset"] != "Blank")
        ]
        #     Precipitation(precipitation_rawfilepath)
        #     if precipitation_rawfilepath
        #     else None
        # )

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
                f"AsT Plate {ast_barcode} has size: {
                    len(ast_plate) // len(ast_plate['AcD Barcode 384'].unique())
                }"
            )
            print(f"{ast_barcode} -> {ast_plate['AcD Barcode 384'].unique()}")
        return result_df


    @cached_property
    def processed(self):
        # TODO: Add precipitation data :)
        # if self.precipitation:
        #     print(
        #         self.precipitation.results.loc[
        #             :, ["Row_384", "Col_384", "AcD Barcode 384", "Precipitated"]
        #         ]
        #     )
        return preprocess(
            self.mapped_input_df,
            substance_id=self._substance_id,
            measurement=self._measurement_label.strip(
                "Raw "
            ),  # I know this is weird, its because of how background_normalize_zfactor works,
            negative_controls=self._negative_controls,
            blanks=self._blanks,
            norm_by_barcode=self._norm_by_barcode,
        )

    @cached_property
    def plateheatmap(self):
        return plateheatmaps(
            self.processed,
            substance_id=self._substance_id,
            negative_control=self._negative_controls,
            blank=self._blanks,
            barcode=self._norm_by_barcode,
        )

    @cached_property
    def _resultfigures(self):
        result_figures = []
        result_figures.append(
            Result("QualityControl", "plateheatmaps", figure=self.plateheatmap)
        )

        for threshold in self.thresholds:
            subset = get_thresholded_subset(
                self._processed_only_substances,
                id_column=self._substance_id,
                negative_controls=self._negative_controls,
                blanks=self._blanks,
                threshold=threshold,
            )
            for dataset, sub_df in subset.groupby("Dataset"):
                dummy_df = get_upsetplot_df(sub_df, counts_column=self._substance_id)

                result_figures.append(
                    Result(
                        dataset,
                        f"UpSetPlot_{dataset}",
                        figure=UpSetAltair(dummy_df, title=dataset),
                    )
                )
        return result_figures

    @cached_property
    def _resulttables(self):
        """
        Retrieves result tables and returns them like list[Resulttable]
        where Resulttable is a dataclass collecting meta information about the plot.
        """
        # result_plots = dict() # {"filepath": plot}
        result_tables = []

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
        for threshold in self.thresholds:
            pivot_df[f"Relative Growth < {threshold}"] = pivot_df.groupby(
                ["Internal ID", "Organism", "Dataset"]
            )["Relative Optical Density mean"].transform(lambda x: x < threshold)

            for dataset, dataset_grp in pivot_df.groupby("Dataset"):
                # dataset = dataset[0]
                # resultpath = os.path.join(filepath, dataset)
                # result_tables[f"{dataset}_all_results"] = dataset_grp
                result_tables.append(
                    Result(dataset, f"{dataset}_all_results", table=dataset_grp)
                )

                pivot_multiindex_df = pd.pivot_table(
                    dataset_grp,
                    values=["Relative Optical Density mean"],
                    index=["Internal ID", "Dataset", "Concentration"],
                    columns="Organism",
                ).reset_index()
                cols = list(pivot_multiindex_df.columns.droplevel())
                cols[:3] = list(map(lambda x: x[0], pivot_multiindex_df.columns[:3]))
                pivot_multiindex_df.columns = cols

                # Apply threshold (active in any organism)
                thresholded_pivot = pivot_multiindex_df.iloc[
                    list(
                        pivot_multiindex_df.iloc[:, 3:].apply(
                            lambda x: any(list(map(lambda i: i < threshold, x))), axis=1
                        )
                    )
                ]

                thresholded_pivot = pivot_multiindex_df.iloc[
                    list(
                        pivot_multiindex_df.iloc[:, 3:].apply(
                            lambda x: any(list(map(lambda i: i < threshold, x))), axis=1
                        )
                    )
                ]

                # Sort by columns each organism after the other
                # return pivot_multiindex_df.sort_values(by=cols[3:])

                # Sort rows by mean between the organisms (lowest mean activity first)
                results_sorted_by_mean_activity = thresholded_pivot.iloc[
                    thresholded_pivot.iloc[:, 3:].mean(axis=1).argsort()
                ]
                # result_tables[f"{dataset}_threshold{round(threshold)}_results"] = results_sorted_by_mean_activity
                result_tables.append(
                    Result(
                        dataset,
                        f"{dataset}_threshold{round(threshold)}_results",
                        table=results_sorted_by_mean_activity,
                    )
                )
        return result_tables

    @cached_property
    def results(self):
        """
        Retrieves result tables (from self._resulttables)
        and returns them in a dictionary like:
            {"<filepath>": pd.DataFrame}
        """
        return {tbl.file_basename: tbl.table for tbl in self._resulttables}

    # @cached_property
    # def precipitation(self, precip_data):


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
