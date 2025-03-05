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
    position_to_rowcol,
    mapapply_96_to_384,
    get_upsetplot_df,
    get_mapping_dict,
    add_precipitation,
    _save_tables,
    _save_figures,
    mic_assaytransfer_mapping,
    get_minimum_precipitation_conc,
)
from .parser import parse_readerfiles, read_inputfile, parse_mappingfile, read_platemapping
from .process import preprocess, get_thresholded_subset, mic_process_inputs
from .plot import plateheatmaps, UpSetAltair, lineplots_facet


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

    def __init__(self, rawfiles_folderpath: str | None, plate_type: int):
        self._plate_type = plate_type
        self._rows, self._columns = get_rows_cols(plate_type)
        self._rawfiles_folderpath = rawfiles_folderpath
        self.rawdata = parse_readerfiles(
            rawfiles_folderpath
        )  # Get rawdata, this will later be overwritting by adding precipitation, if available


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
                    columns = [
                        f"Row_{self._plate_type}",
                        f"Col_{self._plate_type}"
                    ],
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
        elif type(background_locations) is dict:  # TODO: specify plate specific (w. barcode) background_locations using a dictionary
            # {"Barcode": "", "Position": }
            pass


        if not self._rawfiles_folderpath:
            # return early with placeholder results
            self.rawdata_w_layout = pd.DataFrame({
                'Row_384': [],
                'Col_384': [],
                'Raw Optical Density': [],
                'AcD Barcode 384': [],
                'Layout': [],
                'Limit of Quantification': [],
                'Precipitated': [],
                f'Precipitated at {self._measurement_label}': []
                })
        else:
            self.rawdata_w_layout = pd.merge(
                self.rawdata, self.background_locations, how="outer"
            ).fillna({"Layout": "Substance"})

        # if self.rawdata_w_layout[]
        self._background_median = self.rawdata_w_layout[
            self.rawdata_w_layout["Layout"] == "Background"
        ][self._measurement_label].median()
        # Determine the outlier using 2*median of background samples
        self._outlier = self.rawdata_w_layout[
            (self.rawdata_w_layout["Layout"] == "Background") &
            (self.rawdata_w_layout[self._measurement_label] > self._background_median * 2)
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
        self.rawdata_w_layout[f"Precipitated at {self._measurement_label}"] = self.rawdata_w_layout[
            self._measurement_label
        ].apply(lambda x: x if x > self.limit_of_quantification else None)
        # if not self._rawfiles_folderpath:
        #     # return early with placeholder results
        #     return pd.DataFrame({
        #         'Row_384': [],
        #         'Col_384': [],
        #         'Raw Optical Density': [],
        #         'AcD Barcode 384': [],
        #         'Layout': [],
        #         'Limit of Quantification': [],
        #         'Precipitated': [],
        #         f'Precipitated at {self._measurement_label}': []
        #         })
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
        background_locations: pd.DataFrame | list[str] = [f"{row}24" for row in string.ascii_uppercase[:16]],
        precip_exclude_outlier: bool = False,
    ):
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
                    ["Internal ID", "AsT Barcode 384", "Row_384", "Col_384"]
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
                        f"Precipitated at {measurement_label}"
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
        # if not all(column in self._substances_unmapped.columns for column in necessary_columns):
        #     raise ValueError(
        #         f"Not all necessary columns are present in the input table.\n(Necessary columns: {necessary_columns})"
        #     )
        # # Check if all of the necessary column are complete:
        # if substances[necessary_columns].isnull().values.any():
        #     raise ValueError(
        #         "Input table incomplete, contains NA (missing) values."
        #     )
        # # Check if there are duplicates in the internal IDs (apart from references)
        # if any(substances[substances["Dataset"] != "Reference"]["Internal ID"].duplicated()):
        #     raise ValueError("Duplicate Internal IDs.")
        pass

    @cached_property
    def mapped_input_df(self):
        """
        Does mapping of the inputfile describing the tested substances with the
        corresponding mappingfile(s).
        *Basically replaces rda.process.primary_process_inputs() function so all the variables and intermediate results are available via the class*
        """
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
        # Add QualityControl overview of the plates as heatmaps:
        result_figures.append(
            Result("QualityControl", "plateheatmaps", figure=self.plateheatmap)
        )
        # If precipitation testing was done, add it to QC result figures:
        if self.precipitation:
            result_figures.append(
                Result(
                    "QualityControl",
                    "Precipitation_Heatmap",
                    figure=self.precipitation.plateheatmap(),
                )
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
                if not self.precipitation.results.empty:
                    dataset_grp = pd.merge(dataset_grp, self.substances_precipitation)
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

                if not self.precipitation.results.empty:
                    results_sorted_by_mean_activity = pd.merge(results_sorted_by_mean_activity, self.substances_precipitation)
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

    def save_figures(self, resultpath, fileformats: list[str] = ["svg", "html"]):
        _save_figures(resultpath, self._resultfigures, fileformats=fileformats)

    def save_tables(self, resultpath, fileformats: list[str] = ["xlsx", "csv"]):
        _save_tables(resultpath, self._resulttables, fileformats=fileformats)

    def save_results(
        self,
        tables_path: str,
        figures_path: str,
        figureformats: list[str] = ["svg", "html"],
        tableformats: list[str] = ["xlsx", "csv"],
    ):
        self.save_figures(figures_path, fileformats=figureformats)
        self.save_tables(tables_path, fileformats=tableformats)


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
        thresholds: list[float] = [50.0],
        exclude_negative_zfactors: bool = True,
        precipitation_rawfilepath: str | None = None,
        precip_background_locations: pd.DataFrame | list[str] = [f"{row}24" for row in string.ascii_uppercase[:16]],
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
        self.rawdata = (  # Overwrite rawdata if precipitation data is available
            # self.rawdata
            # if self.precipitation is None
            # else
            add_precipitation(
                self.rawdata, self.precipitation.results, self._mapping_dict
            )
        )
        self._substances_unmapped, self._organisms, self._dilutions, self._controls = (
            read_inputfile(inputfile_path)
        )
        self._substance_id = substance_id
        self._negative_controls = negative_controls
        self._blanks = blanks
        self._norm_by_barcode = norm_by_barcode
        self.thresholds = thresholds
        self._processed_only_substances = self.processed[ # Negative Control is still there!
            (self.processed["Dataset"] != "Reference")
            & (self.processed["Dataset"] != "Positive Control")
            & (self.processed["Dataset"] != "Blank")
        ]
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
        self.substances_minimum_precipitation_conc = ( # minimum precipitation conc = MPC
            None
            if (self.precipitation.results.empty) and (self.substances_precipitation)
            else (
                self.substances_precipitation
                .groupby(["Internal ID", "AsT Barcode 384"])
                .apply(
                    lambda x: get_minimum_precipitation_conc(
                        x, precip_conc_multiplicator=precip_conc_multiplicator
                    ),
                    include_groups=False,
                )
                .reset_index(name="Minimum Precipitation Conc.")
                .drop(columns="AsT Barcode 384")
            )
        )
        self._exclude_negative_zfactor = exclude_negative_zfactors

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


        organisms = list(self._organisms["Organism"])

        ast_platemapping, _ = read_platemapping(
            self._mp_ast_mapping_filepath, self._substances_unmapped["MP Barcode 96"].unique()
        )
        # Do some sanity checks:
        necessary_columns = [
            "Dataset",
            "Internal ID",
            "MP Barcode 96",
            "MP Position 96",
        ]
        # Check if all necessary column are present in the input table:
        if not all(column in self._substances_unmapped.columns for column in necessary_columns):
            raise ValueError(
                f"Not all necessary columns are present in the input table.\n(Necessary columns: {necessary_columns})"
            )
        # Check if all of the necessary column are complete:
        if self._substances_unmapped[necessary_columns].isnull().values.any():
            raise ValueError(
                "Input table incomplete, contains NA (missing) values."
            )
        # Check if there are duplicates in the internal IDs (apart from references)
        if any(self._substances_unmapped[self._substances_unmapped["Dataset"] != "Reference"]["Internal ID"].duplicated()):
            raise ValueError("Duplicate Internal IDs.")

        # Map AssayTransfer barcodes to the motherplate barcodes:
        self._substances_unmapped["Row_384"], self._substances_unmapped["Col_384"], self._substances_unmapped["AsT Barcode 384"] = (
            zip(
                *self._substances_unmapped.apply(
                    lambda row: mic_assaytransfer_mapping(
                        row["MP Position 96"],
                        row["MP Barcode 96"],
                        ast_platemapping,
                    ),
                    axis=1,
                )
            )
        )
        acd_platemapping, replicates_dict = read_platemapping(
                self._ast_acd_mapping_filepath, self._substances_unmapped["AsT Barcode 384"].unique()
        )
        num_replicates = list(set(replicates_dict.values()))[0]
        single_subst_concentrations = []

        for substance, subst_row in self._substances_unmapped.groupby("Internal ID"):
            # Collect the concentrations each as rows for a single substance:
            single_subst_conc_rows = []
            init_pos = int(subst_row["Col_384"].iloc[0]) - 1
            col_positions_384 = [list(range(1, 23, 2)), list(range(2, 23, 2))]
            for col_i, conc in enumerate(list(self._dilutions["Concentration"].unique())):
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
        for ast_barcode, ast_plate in input_w_concentrations.groupby("AsT Barcode 384"):
            self._controls["AsT Barcode 384"] = list(ast_plate["AsT Barcode 384"].unique())[0]

            ast_plate = pd.concat([ast_plate, self._controls.copy()])
            for org_i, organism in enumerate(organisms):
                for replicate in range(num_replicates):
                    # Add the AcD barcode
                    ast_plate["AcD Barcode 384"] = acd_platemapping[ast_barcode][
                        replicate
                    ][org_i]

                    ast_plate["Replicate"] = replicate + 1
                    # Add the scientific Organism name
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
        if self.precipitation:
            result_figures.append(
                Result(
                    "QualityControl",
                    "Precipitation_Heatmap",
                    figure=self.precipitation.plateheatmap(),
                )
            )

        # Save plots per dataset:

        processed_negative_zfactor = self._processed_only_substances[self._processed_only_substances["Z-Factor"] < 0]
        if not processed_negative_zfactor.empty and self._exclude_negative_zfactor == True:
            print(
            f"{len(processed_negative_zfactor["AsT Barcode 384"].unique())} plate(s) with negative Z-Factor detected for organisms '{", ".join(processed_negative_zfactor["Organism"].unique())}'.\n",
                "These plates will be excluded from the lineplots visualization!\n (If you want to include them, use the `exclude_negative_zfactors=False` flag of the MIC class)"
            )

        for dataset, dataset_data in self._processed_only_substances.groupby("Dataset"):
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

            result_figures.append(
                Result(
                    dataset,
                    f"{dataset}_lineplots_facet",
                    figure=lineplots_facet(
                        pd.concat([dataset_data, corresponding_dataset_references]),
                        exclude_negative_zfactors=self._exclude_negative_zfactor
                    ),
                )
            )

        # Save plots per threshold:
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
                        f"{dataset}_UpSetPlot",
                        figure=UpSetAltair(dummy_df, title=dataset),
                    )
                )
        return result_figures

    @cached_property
    def _resulttables(self) -> list[Result]:
        """
        Retrieves result tables and returns them like list[Result]
        where Resulttable is a dataclass collecting meta information about the plot.
        """
        result_tables = []
        df = self.processed.copy()

        def get_mic_df(df):
            df = df[
                (df["Dataset"] != "Negative Control") & (df["Dataset"] != "Blank")
            ].dropna(subset=["Concentration"])

            pivot_df = pd.pivot_table(
                df,
                values=["Relative Optical Density", "Replicate", "Z-Factor"],
                index=[
                    "Internal ID",
                    "External ID",
                    "Organism",
                    "Concentration",
                    "Dataset",
                ],
                aggfunc={
                    "Relative Optical Density": ["mean"],
                    "Replicate": ["count"],
                    "Z-Factor": [
                        "mean",
                        "std",
                    ],  # does this make sense? with std its usable.
                    # "Z-Factor": ["std"],
                },
            ).reset_index()

            pivot_df.columns = [" ".join(x).strip() for x in pivot_df.columns.ravel()]

            mic_records = []
            for group_names, grp in pivot_df.groupby(
                ["Internal ID", "External ID", "Organism", "Dataset"]
            ):
                internal_id, external_id, organism, dataset = group_names
                # Sort by concentration just to be sure:
                grp = grp[
                    [
                        "Concentration",
                        "Relative Optical Density mean",
                        "Z-Factor mean",
                        "Z-Factor std",
                    ]
                ].sort_values(by=["Concentration"])

                # Get rows where the OD is below the given threshold:
                record = {
                    "Internal ID": internal_id,
                    "External ID": external_id,
                    "Organism": organism,
                    "Dataset": dataset,
                    "Z-Factor mean": list(grp["Z-Factor mean"])[0],
                    "Z-Factor std": list(grp["Z-Factor std"])[0],
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
            return mic_df

        def get_ref_mics(df):
            only_references = df[
                    df["Dataset"] == "Reference"
                    ]
            mic_records = []
            for group_names, grp in only_references.groupby(
                [
                    "Internal ID",
                    "External ID",
                    "Organism",
                    "Dataset",
                    "AcD Barcode 384",
                ]
            ):
                internal_id, external_id, organism, dataset, acd_barcode = group_names
                grp = grp.copy().sort_values(by=["Concentration"])
                record = {
                    "Internal ID": internal_id,
                    "External ID": external_id,
                    "Organism": organism,
                    "Dataset": dataset,
                    "AcD Barcode 384": acd_barcode,
                    "Z-Factor": list(grp["Z-Factor"])[0],
                }
                for threshold in self.thresholds:
                    values_below_threshold = grp[
                        grp["Relative Optical Density"] < threshold
                    ]
                    # thx to jonathan - check if the OD at maximum concentration is below threshold (instead of any concentration)
                    max_conc_below_threshold = list(
                        grp[grp["Concentration"] == max(grp["Concentration"])][
                            "Relative Optical Density"
                        ]
                        < threshold
                    )[0]
                    if not max_conc_below_threshold:
                        mic = None
                    else:
                        mic = values_below_threshold.iloc[0]["Concentration"]
                    record[f"MIC{threshold} in µM"] = mic
                mic_records.append(record)
            reference_mics = pd.DataFrame.from_records(mic_records)
            return reference_mics

        mic_df = get_mic_df(df)
        references_mic_results = get_ref_mics(df)

        result_tables.append(Result(
            "Reference",
            "References_MIC_results_eachRefID",
            table=references_mic_results,
        ))

        # If precipitation has been done, merge MPC results on long mic_df
        if not self.precipitation.results.empty:
            mic_df = pd.merge(mic_df, self.substances_minimum_precipitation_conc)

        result_tables.append(
            Result("All", "MIC_Results_AllDatasets_longformat", table=mic_df)
        )

        for dataset, dataset_grp in mic_df.groupby("Dataset"):
            pivot_multiindex_df = pd.pivot_table(
                dataset_grp,
                values=[f"MIC{threshold} in µM" for threshold in self.thresholds]
                + ["Z-Factor mean", "Z-Factor std"],
                index=["Internal ID", "External ID", "Dataset"],
                columns="Organism",
            ).reset_index()
            # self.pivot_multiindex_df = pivot_multiindex_df
            # print(pivot_multiindex_df)

            for threshold in self.thresholds:
                organisms_thresholded_mics = pivot_multiindex_df[
                    ["Internal ID", "External ID", f"MIC{threshold} in µM"]
                ]
                cols = list(organisms_thresholded_mics.columns.droplevel())
                cols[0] = "Internal ID"
                cols[1] = "External ID"
                organisms_thresholded_mics.columns = cols
                organisms_thresholded_mics = organisms_thresholded_mics.sort_values(
                    by=list(organisms_thresholded_mics.columns)[2:],
                    na_position="last",
                )

                # Fill with nan if not available
                organisms_thresholded_mics = organisms_thresholded_mics.astype(str)
                # organisms_thresholded_mics.fillna("NA", inplace=True)

                if not self.precipitation.results.empty:
                    organisms_thresholded_mics = pd.merge(
                        organisms_thresholded_mics, self.substances_minimum_precipitation_conc
                    )
                result_tables.append(
                    Result(
                        dataset,
                        f"{dataset}_MIC{round(threshold)}_results",
                        table=organisms_thresholded_mics,
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

    def save_tables(self, resultpath, fileformats: list[str] = ["xlsx", "csv"]):
        _save_tables(resultpath, self._resulttables, fileformats=fileformats)

    def save_results(
        self,
        tables_path: str,
        figures_path: str,
        figureformats: list[str] = ["svg", "html"],
        tableformats: list[str] = ["xlsx", "csv"],
    ):
        self.save_figures(figures_path, fileformats=figureformats)
        self.save_tables(tables_path, fileformats=tableformats)
