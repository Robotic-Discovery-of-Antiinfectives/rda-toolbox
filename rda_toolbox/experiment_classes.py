#!/usr/bin/env python

import pandas as pd
import altair as alt
from functools import cached_property
from dataclasses import dataclass
import numpy as np
import logging
from typing import Optional, List, Dict

import os
import pathlib

import string


from .utility import (
    get_rows_cols,
    position_to_rowcol,
    mapapply_96_to_384,
    get_upsetplot_df,
    get_mapping_dict,
    add_precipitation,
    _save_tables,
    _save_figures,
    mic_assaytransfer_mapping,
    get_minimum_precipitation_conc,
    check_activity_conditions,
    split_position
)
from .parser import (
    parse_readerfiles,
    read_inputfile,
    parse_mappingfile,
    read_platemapping,
)
from .process import preprocess, get_thresholded_subset, mic_process_inputs, add_b_score
from .plot import (
    plateheatmaps,
    UpSetAltair,
    lineplots_facet,
    potency_distribution,
    measurement_vs_bscore_scatter,
    ChartLike
)

logger = logging.getLogger(__name__)

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

    def __init__(self, rawfiles_folderpath: Optional[str], plate_type: int):
        self._plate_type = plate_type
        self._rows, self._columns = get_rows_cols(plate_type)
        self._rawfiles_folderpath = rawfiles_folderpath
        # If no path is provided, initialize empty placeholders instead of calling parse_readerfiles
        if not rawfiles_folderpath:
            self.rawdata = pd.DataFrame()
            self.metadata = {}
        else:
            self.rawdata, self.metadata = parse_readerfiles(
                rawfiles_folderpath
            )  # Get rawdata, this will later be overwritten by adding precipitation, if available

    # def save_results(
    #     self,
    #     save_path: str,
    #     figureformats: list[str] = ["svg", "html"],
    #     tableformats: list[str] = ["xlsx", "csv"]
    # ):
        # pathlib.Path(save_path)


# def save(self, resultpath: str):
#     self.save_plots()
#     self.save_tables()


@dataclass
class Result:
    dataset: str
    file_basename: str
    table: pd.DataFrame | None = None
    figure: ChartLike | None = None


class Precipitation(Experiment):
    def __init__(
        self,
        rawfiles_folderpath: str | None,
        background_locations: pd.DataFrame | list[str] | dict,
        plate_type: int = 384,  # Define default plate_type for experiment
        measurement_label: str = "Raw Optical Density",
        exclude_outlier: bool = False,
    ):
        super().__init__(rawfiles_folderpath, plate_type)
        self._measurement_label = measurement_label

        if type(background_locations) is list:
            self.background_locations = pd.DataFrame(
                list(map(position_to_rowcol, background_locations)),
                columns=[f"Row_{self._plate_type}", f"Col_{self._plate_type}"],
            )
            if "Layout" not in self.background_locations:
                self.background_locations["Layout"] = "Background"
        elif type(background_locations) is pd.DataFrame:
            if "Layout" not in background_locations:
                background_locations["Layout"] = "Background"
            self.background_locations = background_locations.rename(
                columns={
                    "Row": f"Row_{self._plate_type}",
                    "Column": f"Col_{self._plate_type}",
                }
            )
        elif (
            type(background_locations) is dict
        ):  # TODO: specify plate specific (w. barcode) background_locations using a dictionary
            # {"Barcode": "", "Position": }
            pass

        if not self._rawfiles_folderpath:
            # return early with placeholder results
            self.rawdata_w_layout = pd.DataFrame(
                {
                    "Row_384": [],
                    "Col_384": [],
                    "Raw Optical Density": [],
                    "AcD Barcode 384": [],
                    "Layout": [],
                    "Limit of Quantification": [],
                    "Precipitated": [],
                    f"Precipitated at {self._measurement_label}": [],
                }
            )
        else:
            self.rawdata_w_layout = pd.merge(
                self.rawdata, self.background_locations, how="outer"
            ).fillna({"Layout": "Substance"})

        self._background_median = self.rawdata_w_layout[
            self.rawdata_w_layout["Layout"] == "Background"
        ][self._measurement_label].median()
        # Determine the outlier using 2*median of background samples
        self._outlier = self.rawdata_w_layout[
            (self.rawdata_w_layout["Layout"] == "Background")
            & (
                self.rawdata_w_layout[self._measurement_label]
                > self._background_median * 2
            )
        ]
        if not self._outlier.empty:
            if exclude_outlier:
                # Print some info on exluded outliers:
                print("For precipitation test:")
                for index, row in self._outlier.iterrows():
                    print(
                        f"    Exluding outlier on plate {row['AcD Barcode 384']}, position {row['Row_384']}{row['Col_384']}"
                    )
                self.rawdata_w_layout.drop(self._outlier.index, inplace=True)
            else:
                print(
                    "Precipitation background outliers detected,\n",
                    "consider excluding them explicitly via the `background_locations` argument specification\n",
                    "or setting the `exclude_outlier` argument flag to `True`",
                    sep="",
                )
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
        self.rawdata_w_layout[f"Precipitated at {self._measurement_label}"] = (
            self.rawdata_w_layout[self._measurement_label].apply(
                lambda x: x if x > self.limit_of_quantification else None
            )
        )
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
        ast_barcode_header: str = "AsT Barcode 384",
        thresholds: list[float] | None = None,
        b_score_threshold: float = -3.0,
        precipitation_rawfilepath: str | None = None,
        background_locations: pd.DataFrame | list[str] | None = None,
        precip_exclude_outlier: bool = False,
        needs_mapping: bool = True,
    ):
        super().__init__(rawfiles_folderpath, plate_type)
        self._measurement_label = measurement_label
        self._mappingfile_path = mappingfile_path
        self._inputfile_path = inputfile_path
        self._substances_unmapped, self._organisms, self._dilutions, self._controls = (
            read_inputfile(inputfile_path, substance_id)
        )
        self.substances = mapapply_96_to_384(
            self._substances_unmapped,
            rowname=map_rowname,
            colname=map_colname,
            q_name=q_name,
        ) if needs_mapping else split_position(
            self._substances_unmapped,
            position="MP Position 384",
            row="Row_384",
            col="Col_384"
            )

        self._mapping_df = parse_mappingfile(
            mappingfile_path,
            motherplate_column=ast_barcode_header,
            childplate_column="AcD Barcode 384",
        )
        self._mapping_dict = get_mapping_dict(self._mapping_df, mother_column=ast_barcode_header)
        # self._substance_id = substance_id
        self._negative_controls = negative_controls
        self._blanks = blanks
        self._norm_by_barcode = norm_by_barcode
        self._ast_barcode_header = ast_barcode_header
        if thresholds is None:
            thresholds = [50.0]
        self.thresholds = thresholds
        self.b_score_threshold = b_score_threshold
        if background_locations is None:
            background_locations = [
                f"{row}24" for row in string.ascii_uppercase[:16]
            ]

        self.precipitation = (
            Precipitation(
                precipitation_rawfilepath,
                background_locations=background_locations,
                exclude_outlier=precip_exclude_outlier,
            )
            # if precipitation_rawfilepath
            # else None
        )
        self.rawdata = (  # Overwrite rawdata if precipitation data is available
            # self.rawdata
            # if self.precipitation is None
            # else
            add_precipitation(
                self.rawdata, self.precipitation.results, self._mapping_dict
            )
        )
        self._processed_only_substances = self.processed[
            (self.processed["Dataset"] != "Reference")
            & (self.processed["Dataset"] != "Positive Control")
            & (self.processed["Dataset"] != "Blank")
        ]
        self.substances_precipitation = (
            None
            if self.precipitation.results.empty
            else (
                self._processed_only_substances[
                    self._processed_only_substances["Dataset"] != "Negative Control"
                ]
                .drop_duplicates(
                    ["Internal ID", self._ast_barcode_header, "Row_384", "Col_384"]
                )
                .loc[
                    :,
                    [
                        "Internal ID",
                        # "AsT Barcode 384",
                        # "Row_384",
                        # "Col_384",
                        "Concentration",
                        "Precipitated",
                        f"Precipitated at {measurement_label}",
                    ],
                ]
                .reset_index(drop=True)
            )
        )

    def check_substances(self):
        """
        Do some sanity checks for the substances table.
        - Check if all necessary columns are present.
        - Check if substances contains missing values.
        - Check if there are duplicate Internal IDs (references excluded)
        """
        

    @cached_property
    def mapped_input_df(self):
        """
        Does mapping of the inputfile describing the tested substances with the
        corresponding mappingfile(s).
        *Basically replaces rda.process.primary_process_inputs() function so all the variables and intermediate results are available via the class*
        """
        control_wbarcodes = []
        # multiply controls with number of AsT plates to later merge them with substances df
        for origin_barcode in list(self.substances[self._ast_barcode_header].unique()):
            controls_subdf = self._controls.copy()
            controls_subdf[self._ast_barcode_header] = origin_barcode
            control_wbarcodes.append(controls_subdf)
        controls_n_barcodes = pd.concat(control_wbarcodes)

        ast_plate_df = pd.merge(
            pd.concat([self.substances, controls_n_barcodes]), # concatenate substances and controls
            self._dilutions,  # merge with dilutions on column "Dataset"
            how="outer",
            on="Dataset",  # Explicitly define on which column to merge
        )

        mapped_organisms = pd.merge(self._mapping_df, self._organisms, on="Rack")
        rawdata_mapped_organism = pd.merge(mapped_organisms, self.rawdata)

        result_df = pd.concat(
            [
                pd.merge(org_df, ast_plate_df, how="inner")
                for _, org_df in rawdata_mapped_organism.groupby(
                    "Organism formatted"
                )
            ]
        )

        for ast_barcode, ast_plate in result_df.groupby(self._ast_barcode_header):
            logger.info(
                f"AsT Plate {ast_barcode} has size: {
                    len(ast_plate) // len(ast_plate['AcD Barcode 384'].unique())
                }"
            )
            logger.info(f"{ast_barcode} -> {ast_plate['AcD Barcode 384'].unique()}")
        # result_df = result_df.rename({self._substance_id: "Internal ID"}) # rename whatever substance ID was given to Internal ID
        return result_df

    @cached_property
    def processed(self):
        processed = preprocess(
            self.mapped_input_df,
            substance_id="Internal ID",
            measurement=self._measurement_label.strip(
                "Raw "
            ),  # I know this is weird, its because of how background_normalize_zfactor works,
            negative_controls=self._negative_controls,
            blanks=self._blanks,
            norm_by_barcode=self._norm_by_barcode,
        )

        # Add B-Scores to plates without negative controls and blanks
        proc_wo_controls = processed[
            ~processed["Internal ID"].isin([self._negative_controls, self._blanks])
        ]
        # We add b_scores here since we only want them in a primary screen and preprocess() is used generally
        b_scores = (
            proc_wo_controls.groupby(self._norm_by_barcode)[
                [self._norm_by_barcode, "Row_384", "Col_384", self._measurement_label]
            ]
            .apply(lambda plate_grp: add_b_score(plate_grp))
            .reset_index(drop=True)
        )
        processed = pd.merge(processed, b_scores, how="outer")
        return processed

    @cached_property
    def plateheatmap(self):
        return plateheatmaps(
            self.processed.fillna(""),
            substance_id="Internal ID",
            negative_control=self._negative_controls,
            blank=self._blanks,
            barcode=self._norm_by_barcode,
        )

    @cached_property
    def _resultfigures(self):
        result_figures = []
        # Add QualityControl overview of the plates as heatmaps:
        result_figures.append(
            Result("QualityControl", "plateheatmaps", figure=self.plateheatmap)
        )
        # If precipitation testing was done, add it to QC result figures:
        if not self.precipitation.results.empty:
            result_figures.append(
                Result(
                    "QualityControl",
                    "Heatmap_Precipitation",
                    figure=self.precipitation.plateheatmap(),
                )
            )

        for threshold in self.thresholds:
            result_figures.append(
                Result(
                    "QualityControl",
                    "Scatter_Measurement_vs_BScore_Substances",
                    figure=measurement_vs_bscore_scatter(
                        self._processed_only_substances,
                        measurement_header="Relative Optical Density",
                        measurement_title="Relative Optical Density",
                        bscore_header="b_scores",
                        bscore_title="B-Score",
                        color_header="Dataset",
                        measurement_threshold=threshold,
                        b_score_threshold=self.b_score_threshold,
                    ).facet(row="Organism", column="Dataset"),
                )
            )
            result_figures.append(
                Result(
                    "QualityControl",
                    "Scatter_Measurement_vs_BScore_References",
                    figure=measurement_vs_bscore_scatter(
                        self.processed[
                            self.processed["Dataset"] == "Reference"
                        ].replace({np.nan: None}),
                        measurement_header="Relative Optical Density",
                        measurement_title="Relative Optical Density",
                        bscore_header="b_scores",
                        bscore_title="B-Score",
                        color_header="Dataset",
                        measurement_threshold=threshold,
                        b_score_threshold=self.b_score_threshold,
                    ).facet(row="Organism", column="Dataset"),
                )
            )

            subset = get_thresholded_subset(
                self._processed_only_substances,
                id_column="Internal ID",
                negative_controls=self._negative_controls,
                blanks=self._blanks,
                threshold=threshold,
            )
            for dataset, sub_df in subset.groupby("Dataset"):
                dataset_name = str(dataset)
                dummy_df = get_upsetplot_df(sub_df, counts_column="Internal ID")

                result_figures.append(
                    Result(
                        dataset_name,
                        f"UpSetPlot_{dataset_name}",
                        figure=UpSetAltair(dummy_df, title=dataset_name),
                    )
                )
                # ---
                only_actives = self.results[f"{dataset_name}_all_results"][
                    self.results[f"{dataset_name}_all_results"]
                    .groupby("Organism")["Relative Optical Density mean"]
                    .transform(lambda x: x < threshold)
                ]
                result_figures.append(
                    Result(
                        dataset_name,
                        f"Scatterplot_BScores_{dataset_name}",
                        figure=measurement_vs_bscore_scatter(
                            only_actives, show_area=False
                        ),
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
        # result_tables.append(Result("All", ))

        df = self.processed.copy()
        df = df[
            (df["Dataset"] != "Reference")
            & (df["Dataset"] != "Positive Control")
            & (df["Dataset"] != "Blank")
        ].dropna(subset=["Concentration"])

        pivot_df = pd.pivot_table(
            df,
            values=[
                "Relative Optical Density",
                "Replicate",
                "Z-Factor",
                "Robust Z-Factor",
                "b_scores",
            ],
            index=[
                "Internal ID",
                "Organism formatted",
                "Organism",
                "Concentration",
                "Dataset",
            ],
            aggfunc={
                # We need lists here for MultiIndex, otherwise the returned DataFrame is flat
                "Relative Optical Density": ["mean"],
                "Replicate": ["count"],
                "b_scores": ["mean"],
                "Z-Factor": ["mean"],
                "Robust Z-Factor": ["mean"],
            },
        ).reset_index()
        pivot_df.columns = [" ".join(x).strip() for x in pivot_df.columns.ravel()]

        for threshold in self.thresholds:
            # Apply Threshold to % Growth:
            for dataset, dataset_grp in pivot_df.groupby("Dataset"):
                dataset_name = str(dataset)
                if not self.precipitation.results.empty and not self.substances_precipitation is None:
                    dataset_grp = pd.merge(dataset_grp, self.substances_precipitation, how="outer")

                result_tables.append(
                    Result(dataset_name, f"{dataset_name}_all_results", table=dataset_grp)
                )

                # Apply threshold conditions:
                thresholded_dataset_grp = dataset_grp.groupby("Internal ID").filter(
                    lambda x: check_activity_conditions(
                        x["Relative Optical Density mean"],
                        x["b_scores mean"],
                        threshold,
                        self.b_score_threshold,
                    )
                )

                # Pivot the long table for excel viewability:
                pivot_multiindex_df = pd.pivot_table(
                    thresholded_dataset_grp,
                    values=["Relative Optical Density mean", "Z-Factor mean"],
                    index=["Internal ID", "Dataset", "Concentration"],
                    columns="Organism formatted",
                ).reset_index()


                # Sort rows by mean between the organisms (lowest mean measurement first)
                results_sorted_by_mean_activity = pivot_multiindex_df.loc[
                    pivot_multiindex_df.loc[
                        :,
                        list(
                            filter(
                                lambda x: x[0].startswith("Relative Optical Density"),
                                pivot_multiindex_df.columns,
                            )
                        ),
                    ]
                    .mean(axis=1)
                    .argsort()
                ]
                # Only try to merge if precipitation results exist and the precipitation dataframe is present
                if not self.precipitation.results.empty and self.substances_precipitation is not None:
                    # If pivot_table produced MultiIndex columns, flatten them to single level so pandas.merge works
                    if isinstance(results_sorted_by_mean_activity.columns, pd.MultiIndex):
                        results_sorted_by_mean_activity.columns = [
                            " ".join(col).strip() if isinstance(col, tuple) else col
                            for col in results_sorted_by_mean_activity.columns
                        ]
                    # Merge explicitly on Internal ID to avoid ambiguous/level-mismatch merges
                    results_sorted_by_mean_activity = pd.merge(
                        results_sorted_by_mean_activity,
                        self.substances_precipitation,
                        how="left",
                        on="Internal ID",
                    )

                # Correct "mean" header if its only one replicate (remove 'mean')
                if sum(thresholded_dataset_grp["Replicate count"].unique()) == 1:
                    results_sorted_by_mean_activity = results_sorted_by_mean_activity.rename(
                        columns={
                            "Relative Optical Density mean": "Relative Optical Density",
                            # "b_scores mean": "B-Score",
                        }
                    )

                # results_sorted_by_mean_activity = (
                #     results_sorted_by_mean_activity.rename(
                #         columns={"b_scores mean": "B-Score mean"}
                #     )
                # )

                results_sorted_by_mean_activity = (
                    results_sorted_by_mean_activity.fillna("NA")
                )  # Fill NA for better excel readability

                result_tables.append(
                    Result(
                        dataset_name,
                        f"{dataset_name}_threshold{round(threshold)}_results",
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

    def save_figures(self, resultpath, fileformats: list[str] = ["svg", "html"]):
        _save_figures(resultpath, self._resultfigures, fileformats=fileformats)

    def save_tables(
        self, result_path, processed_path, fileformats: list[str] = ["xlsx", "csv"]
    ):
        pathlib.Path(processed_path).mkdir(parents=True, exist_ok=True)
        self.processed.to_csv(os.path.join(processed_path, "processed.csv"))
        _save_tables(result_path, self._resulttables, fileformats=fileformats)

    def save_results(
        self,
        tables_path: str,
        figures_path: str,
        processed_path: str,
        figureformats: list[str] = ["svg", "html"],
        tableformats: list[str] = ["xlsx", "csv"],
    ):
        self.save_figures(figures_path, fileformats=figureformats)
        self.save_tables(tables_path, processed_path, fileformats=tableformats)


class MIC(Experiment):  # Minimum Inhibitory Concentration
    def __init__(
        self,
        rawfiles_folderpath,
        inputfile_path,
        mp_ast_mapping_filepath,
        ast_acd_mapping_filepath,
        plate_type=384,  # Define default plate_type for experiment
        measurement_label: str = "Raw Optical Density",
        map_rowname: str = "Row_96",
        map_colname: str = "Col_96",
        q_name: str = "Quadrant",
        substance_id: str = "Internal ID",
        negative_controls: str = "Bacteria + Medium",
        blanks: str = "Medium",
        norm_by_barcode: str = "AcD Barcode 384",
        thresholds: list[float] | None = None,
        exclude_negative_zfactors: bool = False,
        precipitation_rawfilepath: str | None = None,
        precip_background_locations: pd.DataFrame | list[str] = [
            f"{row}24" for row in string.ascii_uppercase[:16]
        ],
        precip_exclude_outlier: bool = False,
        precip_conc_multiplicator: float = 2.0,
    ):
        super().__init__(rawfiles_folderpath, plate_type)
        self._inputfile_path = inputfile_path
        self._mp_ast_mapping_filepath = mp_ast_mapping_filepath
        self._ast_acd_mapping_filepath = ast_acd_mapping_filepath
        self._measurement_label = measurement_label
        self.precipitation = (
            Precipitation(
                precipitation_rawfilepath,
                background_locations=precip_background_locations,
                exclude_outlier=precip_exclude_outlier,
            )
            # if precipitation_rawfilepath
            # else None
        )
        self.precip_conc_multiplicator = precip_conc_multiplicator
        self.rawdata = (  # Overwrite rawdata if precipitation data is available
            # self.rawdata
            # if self.precipitation is None
            # else
            add_precipitation(
                self.rawdata, self.precipitation.results, self._mapping_dict
            )
        )
        self._substances_unmapped, self._organisms, self._dilutions, self._controls = (
            read_inputfile(inputfile_path, substance_id)
        )

        self._negative_controls = negative_controls
        self._blanks = blanks
        self._norm_by_barcode = norm_by_barcode
        if thresholds is None:
            thresholds = [50.0]
        self.thresholds = thresholds
        self._processed_only_substances = (
            self.processed[  # Negative Control is still there!
                (self.processed["Dataset"] != "Reference")
                & (self.processed["Dataset"] != "Positive Control")
                & (self.processed["Dataset"] != "Blank")
            ]
        )
        self._references_results = self.processed.loc[
            self.processed["Dataset"] == "Reference"
        ]
        self.substances_precipitation = (
            None
            if self.precipitation.results.empty
            else (
                self._processed_only_substances[
                    self._processed_only_substances["Dataset"] != "Negative Control"
                ]
                .drop_duplicates(
                    ["Internal ID", "AsT Barcode 384", "Row_384", "Col_384"]
                )
                .loc[
                    :,
                    [
                        "Internal ID",
                        "AsT Barcode 384",
                        "Row_384",
                        "Col_384",
                        "Concentration",
                        "Precipitated",
                    ],
                ]
                .reset_index(drop=True)
            )
        )
        def get_min_precip_conc_df(self):
            if (self.precipitation.results.empty) and (not self.substances_precipitation):
                return None
            else:
                precip_grps = []
                # precip_df = self.substances_precipitation
                for (int_id, ast_barcode), grp in self.substances_precipitation.groupby(
                    ["Internal ID", "AsT Barcode 384"]
                    ):
                    grp = grp.sort_values("Concentration")
                    min_precip_conc = None
                    if grp.Precipitated.any():
                        min_precip_conc = grp["Concentration"][grp["Precipitated"].idxmax()] * self.precip_conc_multiplicator
                    grp["Minimum Precipitation Concentration"] = min_precip_conc
                    precip_grps.append(grp)
                precip_df = pd.concat(precip_grps)
                precip_df = precip_df[["Internal ID", "Minimum Precipitation Concentration"]]
                return precip_df
        self.substances_minimum_precipitation_conc = get_min_precip_conc_df(self)
        self._exclude_negative_zfactor = exclude_negative_zfactors
        self.mic_df = self.get_mic_df(
                # self.processed.copy()
            df = self.processed[
                (self.processed["Dataset"] != "Negative Control") & (self.processed["Dataset"] != "Blank")
            ].dropna(subset=["Concentration"]).copy()
        ).reset_index(drop=True)



    @property
    def _mapping_dict(self):
        mp_ast_mapping_dict = get_mapping_dict(
            parse_mappingfile(
                self._mp_ast_mapping_filepath,
                motherplate_column="MP Barcode 96",
                childplate_column="AsT Barcode 384",
            ),
            mother_column="MP Barcode 96",
            child_column="AsT Barcode 384",
        )
        ast_acd_mapping_dict = get_mapping_dict(
            parse_mappingfile(
                self._ast_acd_mapping_filepath,
                motherplate_column="AsT Barcode 384",
                childplate_column="AcD Barcode 384",
            ),
            mother_column="AsT Barcode 384",
            child_column="AcD Barcode 384",
        )
        mapping_dict = {}
        for mp_barcode, ast_barcodes in mp_ast_mapping_dict.items():
            tmp_dict = {}
            for ast_barcode in ast_barcodes:
                tmp_dict[ast_barcode] = ast_acd_mapping_dict[ast_barcode]
            mapping_dict[mp_barcode] = tmp_dict
        return mapping_dict

    @cached_property
    def mapped_input_df(self):
        """
        Does mapping of the inputfile describing the tested substances with the
        corresponding mappingfile(s).
        *Basically replaces rda.process.mic_process_inputs() function so all the variables and intermediate results are available via the class*
        """

        # Sorting of organisms via Rack is **very** important, otherwise data gets attributed to wrong organisms
        organisms = list(self._organisms.sort_values(by="Rack")["Organism"])
        formatted_organisms = list(self._organisms.sort_values(by="Rack")["Organism formatted"])


        ast_platemapping, _ = read_platemapping(
            self._mp_ast_mapping_filepath,
            list(map(str, self._substances_unmapped["MP Barcode 96"].unique())),
        )
        # Do some sanity checks:
        necessary_columns = [
            "Dataset",
            "Internal ID",
            "MP Barcode 96",
            "MP Position 96",
        ]
        # Check if all necessary column are present in the input table:
        if not all(
            column in self._substances_unmapped.columns for column in necessary_columns
        ):
            raise ValueError(
                f"Not all necessary columns are present in the input table.\n(Necessary columns: {necessary_columns})"
            )
        # Check if all of the necessary column are complete:
        # TODO: be more precise on what is missing
        if self._substances_unmapped[necessary_columns].isnull().values.any():
            raise ValueError("Input table incomplete, contains NA (missing) values.")
        # Check if there are duplicates in the internal IDs (apart from references)
        if any(
            self._substances_unmapped[
                self._substances_unmapped["Dataset"] != "Reference"
            ]["Internal ID"].duplicated()
        ):
            raise ValueError("Duplicate Internal IDs.")

        # Map AssayTransfer barcodes to the motherplate barcodes:
        (
            self._substances_unmapped["Row_384"],
            self._substances_unmapped["Col_384"],
            self._substances_unmapped["AsT Barcode 384"],
        ) = zip(
            *self._substances_unmapped.apply(
                lambda row: mic_assaytransfer_mapping(
                    row["MP Position 96"],
                    row["MP Barcode 96"],
                    ast_platemapping,
                ),
                axis=1,
            )
        )
        acd_platemapping, replicates_dict = read_platemapping(
            self._ast_acd_mapping_filepath,
            list(map(str, self._substances_unmapped["AsT Barcode 384"].unique())),
        )
        num_replicates = list(set(replicates_dict.values()))[0]
        single_subst_concentrations = []

        for substance, subst_row in self._substances_unmapped.groupby("Internal ID"):
            # Collect the concentrations each as rows for a single substance:
            single_subst_conc_rows = []
            init_pos = int(subst_row["Col_384"].iloc[0]) - 1
            col_positions_384 = [list(range(1, 23, 2)), list(range(2, 23, 2))]
            for col_i, conc in enumerate(
                list(self._dilutions["Concentration"].unique())
            ):
                # Add concentration:
                subst_row["Concentration"] = conc
                # Add corresponding column:
                subst_row["Col_384"] = int(col_positions_384[init_pos][col_i])
                single_subst_conc_rows.append(subst_row.copy())

            # Concatenate all concentrations rows for a substance in a dataframe
            single_subst_concentrations.append(pd.concat(single_subst_conc_rows))
        # Concatenate all self._substances_unmapped dataframes to one whole
        input_w_concentrations = pd.concat(single_subst_concentrations)

        acd_dfs_list = []
        # for ast_barcode, ast_plate in input_w_concentrations.groupby("AsT Barcode 384"):
        #     self._controls["AsT Barcode 384"] = list(
        #         ast_plate["AsT Barcode 384"].unique()
        #     )[0]

        #     ast_plate = pd.concat([ast_plate, self._controls.copy()])
        for ast_barcode, ast_plate in input_w_concentrations.groupby("AsT Barcode 384"):
            controls_with_barcode = self._controls.assign(**{"AsT Barcode 384": ast_barcode})
            ast_plate = pd.concat([ast_plate, controls_with_barcode], ignore_index=True)

            for org_i, organism in enumerate(organisms):
                for replicate in range(num_replicates):
                    # Add the AcD barcode
                    ast_plate["AcD Barcode 384"] = acd_platemapping[ast_barcode][
                        replicate
                    ][org_i]

                    ast_plate["Replicate"] = replicate + 1
                    # Add the scientific Organism name
                    ast_plate["Organism formatted"] = formatted_organisms[org_i]
                    ast_plate["Organism"] = organism
                    acd_dfs_list.append(ast_plate.copy())
                    # Add concentrations:
        acd_single_concentrations_df = pd.concat(acd_dfs_list)

        # merge rawdata with input specifications
        df = pd.merge(self.rawdata, acd_single_concentrations_df, how="outer")
        return df

    @cached_property
    def processed(self):
        return preprocess(
            self.mapped_input_df,
            substance_id="Internal ID",
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
            substance_id="Internal ID",
            barcode=self._norm_by_barcode,
            negative_control=self._negative_controls,
            blank=self._blanks,
        )

    # def lineplots_facet(self):
    #    return lineplots_facet(self.processed)

    @cached_property
    def _resultfigures(self) -> list[Result]:
        result_figures = []
        result_figures.append(
            Result("QualityControl", "plateheatmaps", figure=self.plateheatmap)
        )
        if (self.substances_precipitation is not None) and (
            not self.substances_precipitation.empty
        ):
            result_figures.append(
                Result(
                    "QualityControl",
                    "Precipitation_Heatmap",
                    figure=self.precipitation.plateheatmap(),
                )
            )

        # Save plots per dataset:

        processed_negative_zfactor = self._processed_only_substances[
            self._processed_only_substances["Z-Factor"] < 0
        ]
        if (
            not processed_negative_zfactor.empty
            and self._exclude_negative_zfactor == True
        ):
            print(
                f"{len(processed_negative_zfactor["AsT Barcode 384"].unique())} plate(s) with negative Z-Factor detected for organisms '{", ".join(processed_negative_zfactor["Organism formatted"].unique())}'.\n",
                "These plates will be excluded from the lineplots visualization!\n (If you want to include them, use the `exclude_negative_zfactors=False` flag of the MIC class)",
            )

        for dataset, dataset_data in self._processed_only_substances.groupby("Dataset"):
            dataset_name = str(dataset)
            # Look for and add the corresponding references for each dataset:
            if "AcD Barcode 384" in dataset_data:
                dataset_barcodes = list(dataset_data["AcD Barcode 384"].unique())
                corresponding_dataset_references = self._references_results.loc[
                    (
                        self._references_results["AcD Barcode 384"].isin(
                            dataset_barcodes
                        )
                    ),
                    :,
                ]
            else:
                corresponding_dataset_references = pd.DataFrame()

            lineplots_input_df = pd.concat(
                [dataset_data, corresponding_dataset_references]
            )
            lineplots_input_df = lineplots_input_df.dropna(
                subset=["Concentration"]
            ).loc[
                (lineplots_input_df["Dataset"] != "Negative Control")
                & (lineplots_input_df["Dataset"] != "Blank"),
                :,
            ]
            if not lineplots_input_df.empty:
                for threshold in self.thresholds:
                    result_figures.append(
                        Result(
                            dataset_name,
                            f"{dataset_name}_lineplots_facet_thrsh{threshold}_InternalID",
                            figure=lineplots_facet(
                                lineplots_input_df,
                                by_id="Internal ID",
                                exclude_negative_zfactors=self._exclude_negative_zfactor,
                                threshold=threshold,
                            ),
                        )
                    )
                    result_figures.append(
                        Result(
                            dataset_name,
                            f"{dataset_name}_lineplots_facet_thrsh{threshold}_ExternalID",
                            figure=lineplots_facet(
                                lineplots_input_df,
                                by_id="External ID",
                                exclude_negative_zfactors=self._exclude_negative_zfactor,
                                threshold=threshold,
                            ),
                        )
                    )

        # Save plots per threshold:
        for threshold in self.thresholds:
            for dataset, sub_df in self.mic_df.groupby("Dataset"):
                dataset_name = str(dataset)
                # print(sub_df)
                sub_df = sub_df.dropna(subset=f"MIC{threshold} in µM")
                if sub_df.empty:
                    print(f"No MICs for dataset: {dataset_name}, threshold: {threshold}")
                    continue
                dummy_df = get_upsetplot_df(
                    sub_df,
                    counts_column="Internal ID",
                    set_column="Organism",
                )

                result_figures.append(
                    Result(
                        dataset_name,
                        f"{dataset_name}_UpSetPlot",
                        figure=UpSetAltair(dummy_df, title=dataset_name),
                    )
                )
                result_figures.append(
                    Result(
                        dataset_name,
                        f"{dataset_name}_PotencyDistribution",
                        figure=potency_distribution(sub_df, threshold, dataset_name),
                    )
                )
        return result_figures

    def get_mic_df(self, df):

        pivot_df = pd.pivot_table(
            df,
            values=["Relative Optical Density", "Replicate", "Z-Factor", "Robust Z-Factor"],
            index=[
                "Internal ID",
                # "External ID",
                "Organism formatted",
                "Organism",
                "Concentration",
                "Dataset",
            ],
            aggfunc={
                "Relative Optical Density": ["mean"],
                "Replicate": ["count"],
                "Z-Factor": ["mean", "std"],
                "Robust Z-Factor": ["mean", "std"],
            },
            # margins=True
            fill_value=0 # This might result in confusion, if there are no replicates (1)
        ).reset_index()

        pivot_df.columns = [" ".join(x).strip() for x in pivot_df.columns.ravel()]

        mic_records = []
        for group_names, grp in pivot_df.groupby(
            ["Internal ID", "Organism formatted", "Dataset"]
        ):
            internal_id, organism_formatted, dataset = group_names
            # Sort by concentration just to be sure:
            grp = grp[
                [
                    "Concentration",
                    "Relative Optical Density mean",
                    "Z-Factor mean",
                    "Z-Factor std",
                    "Robust Z-Factor mean",
                    "Robust Z-Factor std",
                ]
            ].sort_values(by=["Concentration"])

            # Get rows where the OD is below the given threshold:
            record = {
                "Internal ID": internal_id,
                "Organism formatted": organism_formatted,
                "Dataset": dataset,
                "Z-Factor mean": list(grp["Z-Factor mean"])[0],
                "Z-Factor std": list(grp["Z-Factor std"])[0],
                "Robust Z-Factor mean": list(grp["Robust Z-Factor mean"])[0],
                "Robust Z-Factor std": list(grp["Robust Z-Factor std"])[0],
            }

            for threshold in self.thresholds:
                values_below_threshold = grp[
                    grp["Relative Optical Density mean"] < threshold
                ]
                # thx to jonathan - check if the OD at maximum concentration is below threshold (instead of any concentration)
                max_conc_below_threshold = list(
                    grp[grp["Concentration"] == max(grp["Concentration"])][
                        "Relative Optical Density mean"
                    ]
                    < threshold
                )[0]
                if not max_conc_below_threshold:
                    mic = None
                else:
                    mic = values_below_threshold.iloc[0]["Concentration"]
                record[f"MIC{threshold} in µM"] = mic
            mic_records.append(record)
        # Drop entries where no MIC could be determined
        mic_df = pd.DataFrame.from_records(mic_records)
        # Merge inconsistent (but maybe necessary) columns again
        mic_df = pd.merge(mic_df, df[["Internal ID", "External ID"]], on=["Internal ID"])
        mic_df = pd.merge(mic_df, self._organisms[["Organism", "Organism formatted"]], on=["Organism formatted"])
        unit_values = self._dilutions.get("Unit")
        mic_df["Unit"] = unit_values.dropna().iloc[0] if unit_values is not None and not unit_values.dropna().empty else None
        mic_df = mic_df.drop_duplicates()
        return mic_df


    @cached_property
    def _resulttables(self) -> list[Result]:
        """
        Retrieves result tables and returns them like list[Result]
        where Resulttable is a dataclass collecting meta information about the plot.
        """
        result_tables = []
        df = self.processed.copy()


        # mic_df = self.get_mic_df(df)
        references_mic_results = self.get_mic_df(
            self.processed[self.processed["Dataset"] == "Reference"].copy()
        ).reset_index(drop=True)

        result_tables.append(
            Result(
                "Reference",
                "References_MIC_results_eachRefID",
                table=references_mic_results,
            )
        )

        mic_df = self.mic_df
        # If precipitation has been done, merge MPC results on long mic_df
        if not self.precipitation.results.empty and not self.substances_minimum_precipitation_conc is None:
            mic_df = pd.merge(self.mic_df, self.substances_minimum_precipitation_conc, on="Internal ID", how="left")

        result_tables.append(
            Result("All", "MIC_Results_AllDatasets_longformat", table=self.mic_df)
        )

        for dataset, dataset_grp in mic_df.groupby("Dataset"):
            dataset_name = str(dataset)
            print(f"Preparing tables for dataset: {dataset_name}")
            
            # dataset_grp = dataset_grp.fillna("NA")
            values_list = [f"MIC{threshold} in µM" for threshold in self.thresholds] + [
                "Z-Factor mean",
                "Z-Factor std",
            ]
            # create pivot without resetting index so we can manipulate the MultiIndex columns
            pivot_df = pd.pivot_table(
                dataset_grp,
                values=values_list,
                index=["Internal ID", "Dataset"],
                columns="Organism",
            )

            # ensure columns are a MultiIndex of (value, organism) even if only one value was pivoted
            if not isinstance(pivot_df.columns, pd.MultiIndex):
                # single value case: pivot_df.columns are organisms, create first level from the single value name
                first_value_name = values_list[0]
                pivot_df.columns = pd.MultiIndex.from_product([[first_value_name], pivot_df.columns.tolist()])

            # determine full set of organisms we want to keep (preserve order from self._organisms if available)
            try:
                organisms = list(self._organisms.sort_values(by="Rack")["Organism"].tolist())
            except Exception:
                organisms = list(self._organisms["Organism"].unique())

            # build expected full MultiIndex and reindex to keep columns that were all-NaN
            expected_columns = pd.MultiIndex.from_product([values_list, organisms])
            pivot_df = pivot_df.reindex(columns=expected_columns)

            # finally reset index to get the same shape as before
            pivot_multiindex_df = pivot_df.reset_index()

            for threshold in self.thresholds:

                if pivot_multiindex_df.empty:
                    continue

                organisms_thresholded_mics = pivot_multiindex_df[
                    ["Internal ID", f"MIC{threshold} in µM"]
                ]
                cols = list(organisms_thresholded_mics.columns.droplevel())
                cols[0] = "Internal ID"
                organisms_thresholded_mics.columns = cols
                organisms_thresholded_mics = organisms_thresholded_mics.sort_values(
                    by=list(organisms_thresholded_mics.columns)[1:],
                    na_position="last",
                )

                # Fill with nan if not available
                organisms_thresholded_mics = organisms_thresholded_mics.round(2)
                organisms_thresholded_mics = organisms_thresholded_mics.astype(str)
                organisms_thresholded_mics = pd.merge(organisms_thresholded_mics, self.mic_df[["Internal ID", "External ID"]], on=["Internal ID"], how="left")
                # organisms_thresholded_mics.fillna("NA", inplace=True)

                if not self.precipitation.results.empty and not self.substances_minimum_precipitation_conc is None:
                    organisms_thresholded_mics = pd.merge(
                        organisms_thresholded_mics,
                        self.substances_minimum_precipitation_conc,
                        how="left"
                    )
                organisms_thresholded_mics = organisms_thresholded_mics.reset_index(drop=True)
                organisms_thresholded_mics = organisms_thresholded_mics.drop_duplicates()
                # Add unit column
                unit_values = self._dilutions.get("Unit")
                organisms_thresholded_mics["Unit"] = unit_values.dropna().iloc[0] if unit_values is not None and not unit_values.dropna().empty else None
                # Reorder columns
                desired_order = ["Internal ID", "External ID"]
                remaining_cols = [col for col in organisms_thresholded_mics.columns if col not in desired_order]
                organisms_thresholded_mics = organisms_thresholded_mics[desired_order + remaining_cols]
                organisms_thresholded_mics = organisms_thresholded_mics.replace("nan", "NA").fillna("NA")
                result_tables.append(
                    Result(
                        dataset_name,
                        f"{dataset_name}_MIC{int(round(threshold))}_results",
                        table=organisms_thresholded_mics.reset_index(drop=True)
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

    def save_figures(self, result_path, fileformats: list[str] = ["svg", "html"]):
        _save_figures(result_path, self._resultfigures, fileformats=fileformats)

    def save_tables(
        self, result_path, processed_path, fileformats: list[str] = ["xlsx", "csv"]
    ):
        # Create folder if not existent:
        pathlib.Path(processed_path).mkdir(parents=True, exist_ok=True)
        self.processed.to_csv(os.path.join(processed_path, "processed.csv"))
        _save_tables(result_path, self._resulttables, fileformats=fileformats)

    def save_results(
        self,
        tables_path: str,
        figures_path: str,
        processed_path: str,
        figureformats: list[str] = ["svg", "html"],
        tableformats: list[str] = ["xlsx", "csv"],
    ):
        self.save_figures(figures_path, fileformats=figureformats)
        self.save_tables(tables_path, processed_path, fileformats=tableformats)
