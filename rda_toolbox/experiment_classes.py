#!/usr/bin/env python

import pandas as pd
import altair as alt
from functools import cached_property

import string


from .utility import get_rows_cols
from .parser import parse_readerfiles


class Experiment:
    def __init__(self, rawfiles_folderpath, plate_type):
        self._plate_type = plate_type
        self._rows, self._columns = get_rows_cols(plate_type)
        self._rawfiles_folderpath = rawfiles_folderpath
        self.rawdata = parse_readerfiles(rawfiles_folderpath)


class Precipitation(Experiment):
    def __init__(
        self,
        rawfiles_folderpath,
        plate_type=384,  # Define default plate_type for experiment
        measurement_label="Raw Optical Density",
        background_locations=None,
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
        rawfiles_folderpath,
        plate_type=384,  # Define default plate_type for experiment
        measurement_label="Raw Optical Density",
        ):
        pass


class MIC(Experiment):  # Minimum Inhibitory Concentration
    def __init__(
        self,
        rawfiles_folderpath,
        plate_type=384,  # Define default plate_type for experiment
        measurement_label="Raw Optical Density",
        ):
        pass
