#!/usr/bin/env python3

import altair as alt
import pandas as pd
import pathlib
from .utility import (
    prepare_visualization,
    get_upsetplot_df,
)
from .process import (
    get_thresholded_subset,
)


def get_heatmap(subdf, substance_id, measurement, negative_controls, blanks):
    base = alt.Chart(
        subdf,
    ).encode(
        alt.X("Col_384:O").axis(labelAngle=0, orient="top").title(None),
        alt.Y("Row_384:O").title(None),
        tooltip=list(subdf.columns),
    )
    blank_mean = subdf[subdf[substance_id] == blanks][measurement].mean()
    negative_mean = subdf[subdf[substance_id] == negative_controls][measurement].mean()
    heatmap = base.mark_rect().encode(
        alt.Color(f"{measurement}:Q")
        .title("Optical Density")
        .scale(domain=[blank_mean, negative_mean]),
    )
    text = base.mark_text(baseline="middle", align="center", fontSize=7).encode(
        alt.Text(f"{measurement}:Q", format=".1f"),
        color=alt.condition(
            alt.datum[measurement]
            < max(subdf[subdf[substance_id] == negative_controls][measurement]) / 2,
            alt.value("black"),
            alt.value("white"),
        ),
    )
    return alt.layer(heatmap, text)


def plateheatmaps(
    df,
    substance_id="ID",
    measurement="Raw Optical Density",
    barcode="Barcode",
    negative_control="Negative Control",
    blank="Medium",
) -> alt.vegalite.v5.api.HConcatChart:
    """
    - df: Dataframe with relevant data
    - ID: column name in df containing the unique substance id
    - measurement: column name in df with the measurements to colorize via heatmaps
    - negative_control: controls with organism + medium
    - blank: controls with only medium (no organism and therefore no growth)
    Plots heatmaps of the plates from df in a gridlike manner.
    Exclude unwanted plates, for example Blanks from the df outside this function, like so
    `df[df["Organism"] != "Blank"]`
    before plotting, otherwise it will appear as an extra plate.
    """
    df["Col_384"] = df["Col_384"].astype(int)
    # df[substance_id] = df[substance_id].astype(str)
    plots = []
    for _, _organism_df in df.groupby("Organism"):
        plots.append(
            get_heatmap(
                _organism_df,
                substance_id,
                measurement,
                negative_control,
                blank,
            )
            .facet(
                row=alt.Row(f"{barcode}:N"),
                column=alt.Column("Replicate:N"),
                title=alt.Title(
                    _organism_df["Organism"].unique()[0],
                    orient="top",
                    anchor="middle",
                    dx=-20,
                ),
            )
            .resolve_scale(color="shared")
            .resolve_axis(x="independent", y="independent")
        )

    plate_heatmaps = (
        alt.hconcat(*plots).resolve_scale(color="independent").resolve_axis(y="shared")
    )
    return plate_heatmaps


def blank_heatmap(blank_df):
    blank_df = blank_df.copy()[["Row_384", "Col_384", "Raw Optical Density"]]
    title = "Blank (Only Medium) plate"
    base = alt.Chart(
        blank_df, title=alt.TitleParams("", subtitle=title) if title else ""
    ).encode(
        alt.X("Col_384:O").axis(labelAngle=0, orient="top").title(None),
        alt.Y("Row_384:O").title(None),
        tooltip=list(blank_df.columns),
    )
    heatmap = base.mark_rect().encode(
        alt.Color("Raw Optical Density:Q")
        .title("Optical Density")
        .scale(domain=[0, 100])
    )
    text = base.mark_text(baseline="middle", align="center", fontSize=10).encode(
        alt.Text("Raw Optical Density:Q", format=".1f"),
    )
    return alt.layer(heatmap, text)


# The following two functions are (modified versions) from https://github.com/hms-dbmi/upset-altair-notebook
# http://upset.plus/covid19symptoms/
# and were provided under the MIT licence
# by Sehi L'yi and Nils Gehlenborg


def upsetaltair_top_level_configuration(
    base, legend_orient="top-left", legend_symbol_size=30
):
    return (
        base.configure_view(stroke=None)
        .configure_title(
            fontSize=18, fontWeight=400, anchor="start", subtitlePadding=10
        )
        .configure_axis(
            labelFontSize=14,
            labelFontWeight=300,
            titleFontSize=16,
            titleFontWeight=400,
            titlePadding=10,
        )
        # .configure_legend(
        #     titleFontSize=16,
        #     titleFontWeight=400,
        #     labelFontSize=14,
        #     labelFontWeight=300,
        #     padding=20,
        #     orient=legend_orient,
        #     symbolType="circle",
        #     symbolSize=legend_symbol_size,
        # )
        .configure_legend(disable=True)
        .configure_concat(spacing=0)
    )


def UpSetAltair(
    data=None,
    title="",
    subtitle="",
    sets=None,
    abbre=None,
    sort_by="frequency",
    sort_by_order="ascending",
    inter_degree_frequency="ascending",
    width=1200,
    height=700,
    height_ratio=0.6,
    horizontal_bar_chart_width=300,
    set_colors_dict=dict(),
    highlight_color="#777777",
    glyph_size=200,
    set_label_bg_size=1000,
    line_connection_size=2,
    horizontal_bar_size=20,
    vertical_bar_label_size=16,
    vertical_bar_padding=20,
    set_labelstyle="normal",
):
    """This function generates Altair-based interactive UpSet plots.

    Parameters:
        - data (pandas.DataFrame): Tabular data containing the membership of each element (row) in
            exclusive intersecting sets (column).
        - sets (list): List of set names of interest to show in the UpSet plots.
            This list reflects the order of sets to be shown in the plots as well.
        - abbre (dict): Dictionary mapping set names to abbreviated set names.
        - sort_by (str): "frequency" or "degree"
        - sort_by_order (str): "ascending" or "descending"
        - inter_degree_frequency: "ascending" or "descending", only makes sense if sort_by="degree"
        - width (int): Vertical size of the UpSet plot.
        - height (int): Horizontal size of the UpSet plot.
        - height_ratio (float): Ratio of height between upper and under views, ranges from 0 to 1.
        - horizontal_bar_chart_width (int): Width of horizontal bar chart on the bottom-right.
        - set_colors_dict (dict): Dictionary containing the sets as keys with corresponding colors as values
        - highlight_color (str): Color to encode intersecting sets upon mouse hover.
        - glyph_size (int): Size of UpSet glyph (⬤).
        - set_label_bg_size (int): Size of label background in the horizontal bar chart.
        - line_connection_size (int): width of lines in matrix view.
        - horizontal_bar_size (int): Height of bars in the horizontal bar chart.
        - vertical_bar_label_size (int): Font size of texts in the vertical bar chart on the top.
        - vertical_bar_padding (int): Gap between a pair of bars in the vertical bar charts.
        - set_labelstyle (str): "normal" (default) or "italic"

    Run rda.utility.get_upsetplot_df() on the df before trying this function.
    """

    if data is None:
        print("No data and/or a list of sets are provided")
        return
    if sets is None:
        sets = list(data.columns[1:])

    if (height_ratio < 0) or (1 < height_ratio):
        print("height_ratio set to 0.5")
        height_ratio = 0.5
    if not abbre:
        abbre = {set: set for set in sets}
    if len(sets) != len(abbre):
        abbre = sets
        print(
            "Dropping the `abbre` list because the lengths of `sets` and `abbre` are not identical."
        )
    if not set_colors_dict:  # build default colors dict
        colors = [  # observable10
            "#4269d0",
            "#efb118",
            "#ff725c",
            "#6cc5b0",
            "#3ca951",
            "#ff8ab7",
            "#a463f2",
            "#97bbf5",
            "#9c6b4e",
            "#9498a0",
        ]
        if len(sets) > len(colors):
            colors = colors * len(sets)
        set_colors_dict = {key: value for key, value in zip(sets, colors[: len(sets)])}
    else:
        if sorted(list(set_colors_dict.keys())) != sorted(sets):
            raise ValueError(
                f"Wrong set names, correct names are:\n{dict((set, '') for set in sets)}"
            )
    # filter set_colors_dict with the sets which are actually in the data df (sets)
    # this might be needed if set_colors_dict if more comprehensive than the data
    set_colors_dict = {
        key: value for key, value in set_colors_dict.items() if key in sets
    }
    """
    Data Preprocessing
    """
    data = data.copy()
    data["count"] = 0
    data = data[sets + ["count"]]
    data = data.groupby(sets).count().reset_index()

    data["intersection_id"] = data.index
    data["degree"] = data[sets].sum(axis=1)
    data = data.sort_values(
        by=["count"],
        ascending=True if inter_degree_frequency == "ascending" else False,
    )

    data = pd.melt(data, id_vars=["intersection_id", "count", "degree"])
    data = data.rename(columns={"variable": "set", "value": "is_intersect"})

    set_to_abbre = pd.DataFrame(abbre.items(), columns=["set", "set_abbre"])

    set_to_order = (
        data[data["is_intersect"] == 1]
        .groupby("set")
        .sum()
        .reset_index()
        .sort_values(by="count", ascending=False)
        .filter(["set"])
    )
    set_to_order["set_order"] = list(range(len(sets)))

    degree_calculation = ""
    for s in sets:
        degree_calculation += f"(isDefined(datum['{s}']) ? datum['{s}'] : 0)"
        if sets[-1] != s:
            degree_calculation += "+"
    """
    Selections
    """
    legend_selection = alt.selection_point(fields=["set"], bind="legend")
    color_selection = alt.selection_point(
        fields=["intersection_id"], on="pointerover", empty=False
    )
    opacity_selection = alt.selection_point(fields=["intersection_id"])

    """
    Styles
    """
    vertical_bar_chart_height = height * height_ratio
    matrix_height = height - vertical_bar_chart_height
    matrix_width = width - horizontal_bar_chart_width

    vertical_bar_size = min(
        30,
        width / len(data["intersection_id"].unique().tolist()) - vertical_bar_padding,
    )

    main_color = "#3A3A3A"
    brush_opacity = alt.condition(~opacity_selection, alt.value(1), alt.value(0.6))
    brush_color = alt.condition(
        color_selection, alt.value(highlight_color), alt.value(main_color)
    )
    is_show_horizontal_bar_label_bg = len(list(abbre.values())[0]) <= 2
    horizontal_bar_label_bg_color = (
        "white" if is_show_horizontal_bar_label_bg else "black"
    )

    x_sort = alt.Sort(
        field="count" if sort_by == "frequency" else "degree",
        order=sort_by_order,
    )

    tooltip = [
        alt.Tooltip("max(count):Q", title="Cardinality"),
        alt.Tooltip("degree:Q", title="Degree"),
    ]
    """
    Plots
    """
    # To use native interactivity in Altair, we are using the data transformation functions
    # supported in Altair.
    base = (
        alt.Chart(data)
        .transform_pivot(
            "set",
            op="max",
            groupby=["intersection_id", "count"],
            value="is_intersect",
        )
        .transform_aggregate(
            # count, set1, set2, ...
            count="sum(count)",
            groupby=sets,
        )
        .transform_calculate(
            # count, set1, set2, ...
            degree=degree_calculation
        )
        .transform_filter(
            # count, set1, set2, ..., degree
            alt.datum["degree"]
            != 0
        )
        .transform_window(
            # count, set1, set2, ..., degree
            intersection_id="row_number()",
            frame=[None, None],
        )
        .transform_fold(
            # count, set1, set2, ..., degree, intersection_id
            sets,
            as_=["set", "is_intersect"],
        )
        .transform_lookup(
            # count, set, is_intersect, degree, intersection_id
            lookup="set",
            from_=alt.LookupData(set_to_abbre, "set", ["set_abbre"]),
        )
        .transform_lookup(
            # count, set, is_intersect, degree, intersection_id, set_abbre
            lookup="set",
            from_=alt.LookupData(set_to_order, "set", ["set_order"]),
        )
        .transform_filter(
            # Make sure to remove the filtered sets.
            legend_selection
        )
        .transform_window(
            # count, set, is_intersect, degree, intersection_id, set_abbre
            set_order="distinct(set)",
            frame=[None, 0],
            sort=[{"field": "set_order"}],
        )
        .transform_lookup(
            lookup="set",
            from_=alt.LookupData(set_to_order, "set", ["set_order"]),
        )
    )

    vertical_bar = (
        base.mark_bar(color=main_color)  # , size=vertical_bar_size)
        .encode(
            x=alt.X(
                "intersection_id:N",
                axis=alt.Axis(grid=False, labels=False, ticks=False, domain=True),
                sort=x_sort,
                title=None,
            ),
            y=alt.Y(
                "max(count):Q",
                axis=alt.Axis(grid=False, tickCount=3, orient="right"),
                title="Intersection Size",
            ),
            color=brush_color,
            tooltip=tooltip,
        )
        .properties(width=matrix_width, height=vertical_bar_chart_height)
    )
    vertical_bar_text = vertical_bar.mark_text(
        color=main_color, dy=-10, size=vertical_bar_label_size, fontSize=20
    ).encode(text=alt.Text("count:Q", format=".0f"))
    vertical_bar_chart = (vertical_bar + vertical_bar_text).add_params(
        color_selection,
    )

    circle_bg = (
        vertical_bar.mark_circle(size=glyph_size, opacity=1)
        .encode(
            x=alt.X(
                "intersection_id:N",
                axis=alt.Axis(grid=False, labels=False, ticks=False, domain=False),
                sort=x_sort,
                title=None,
            ),
            y=alt.Y(
                "set_order:N",
                axis=alt.Axis(grid=False, labels=False, ticks=False, domain=False),
                title=None,
            ),
            color=alt.value("#E6E6E6"),
        )
        .properties(height=matrix_height)
    )
    rect_bg = (
        circle_bg.mark_rect()
        .transform_filter(alt.datum["set_order"] % 2 == 1)
        .encode(color=alt.value("#F7F7F7"))
    )
    circle = circle_bg.transform_filter(alt.datum["is_intersect"] == 1).encode(
        color=brush_color
    )
    line_connection = (
        circle_bg.mark_bar(size=line_connection_size, color=main_color)
        .transform_filter(alt.datum["is_intersect"] == 1)
        .encode(
            y=alt.Y("min(set_order):N"),
            y2=alt.Y2("max(set_order):N"),
            color=brush_color,
        )
    )
    matrix_view = alt.layer(
        circle + rect_bg + circle_bg + line_connection + circle
    ).add_params(
        # Duplicate `circle` is to properly show tooltips.
        color_selection,
    )

    # Cardinality by sets (horizontal bar chart)
    horizontal_bar_label_bg = base.mark_circle(size=set_label_bg_size).encode(
        y=alt.Y(
            "set_order:N",
            axis=alt.Axis(grid=False, labels=False, ticks=False, domain=False),
            title=None,
        ),
        color=alt.Color(
            "set:N",
            scale=alt.Scale(
                domain=list(set_colors_dict.keys()),
                range=list(set_colors_dict.values()),
            ),
            title=None,
        ),
        opacity=alt.value(1),
    )
    horizontal_bar_label = horizontal_bar_label_bg.mark_text(
        align=("center" if is_show_horizontal_bar_label_bg else "center"),
        fontSize=20,
        fontStyle=set_labelstyle,
    ).encode(
        text=alt.Text("set_abbre:N"),
        color=alt.value(horizontal_bar_label_bg_color),
    )
    horizontal_bar_axis = (
        (horizontal_bar_label_bg + horizontal_bar_label)
        if is_show_horizontal_bar_label_bg
        else horizontal_bar_label
    )

    horizontal_bar = (
        horizontal_bar_label_bg.mark_bar(size=horizontal_bar_size)
        .transform_filter(alt.datum["is_intersect"] == 1)
        .encode(
            x=alt.X(
                "sum(count):Q",
                axis=alt.Axis(grid=False, tickCount=3),
                title="Set Size",
                # scale=alt.Scale(range=color_range)
            ),
            # color=alt.Color(None,legend=None), # remove interactivity, color and legend
        )
        .properties(width=horizontal_bar_chart_width)
    )
    horizontal_bar_text = horizontal_bar.mark_text(
        align="left", dx=2, fontSize=20
    ).encode(text="sum(count):Q")
    horizontal_bar_chart = alt.layer(horizontal_bar, horizontal_bar_text)
    # Concat Plots
    upsetaltair = alt.vconcat(
        vertical_bar_chart,
        alt.hconcat(
            matrix_view,
            horizontal_bar_axis,
            horizontal_bar_chart,
            spacing=5,
        ).resolve_scale(y="shared"),
        spacing=20,
    ).add_params(
        legend_selection,
    )

    # Apply top-level configuration
    upsetaltair = upsetaltair_top_level_configuration(
        upsetaltair,
        legend_orient="top",
        legend_symbol_size=set_label_bg_size / 2.0,
    ).properties(
        title={
            "text": title,
            "subtitle": subtitle,
            "fontSize": 20,
            "fontWeight": 500,
            "subtitleColor": main_color,
            "subtitleFontSize": 14,
        }
    )
    return upsetaltair


def UpSet_per_dataset(
    df: pd.DataFrame,  # processed
    save_formats=["pdf", "svg"],
    id_column="Internal ID",
):
    """
    UpsetPlot wrapper function which applies threshold to processed data (without controls, references etc.).
    For each dataset present in the given df, create a dummy_df for rda.UpSetAltair() and save the UpSetPlot.
    """
    subset = get_thresholded_subset(
        df,
        id_column="Internal ID",
        negative_controls="Bacteria + Medium",
        blanks="Medium",
        threshold=50,
    )

    for dataset, sub_df in subset.groupby("Dataset"):
        dummy_df = get_upsetplot_df(sub_df, counts_column=id_column)
        # Create dataset folder if non-existent
        pathlib.Path(f"../figures/{dataset}").mkdir(parents=True, exist_ok=True)
        for save_format in save_formats:
            filename = f"../figures/{dataset}/UpSetPlot_{dataset}.{save_format}"
            print("Saving", filename)
            dataset_upsetplot = rda.UpSetAltair(dummy_df, title=dataset).save(filename)


def lineplots_facet(
    df,
    hline_y=50,
    by_id="Internal ID",
    whisker_width=10,
    exclude_negative_zfactors=True,
) -> alt.vegalite.v5.api.HConcatChart:
    """
    Assay: MIC
    Input: processed_df
    Output: Altair Chart with faceted lineplots.
    Negative controls and blanks are dropped inside the function.
    """
    df = df.dropna(subset=["Concentration"]).loc[
        (df["Dataset"] != "Negative Control") & (df["Dataset"] != "Blank"), :
    ]
    df = prepare_visualization(
        df, by_id=by_id, exclude_negative_zfactors=exclude_negative_zfactors
    )
    hline_y = 50
    organism_columns = []

    color = alt.condition(
        # alt.datum.Concentration
        alt.datum.max_conc_below_threshold,
        alt.Color("External ID:N"),
        alt.value("lightgray"),
    )
    for organism, org_data in df.groupby(["Organism"]):
        base = alt.Chart(org_data).encode(color=color)  # , title=organism)
        lineplot = base.mark_line(point=True, size=0.8).encode(
            x=alt.X(
                "Concentration:O",
                title="Concentration in µM",
                axis=alt.Axis(labelAngle=-45),
            ),
            y=alt.Y(
                "Mean Relative Optical Density:Q",
                title="Relative Optical Density",
                scale=alt.Scale(domain=[-20, 160], clamp=True),
            ),
            # color="Internal ID:N",
            shape=alt.Shape("External ID:N", legend=None),
            # color=color,
            tooltip=[
                "Internal ID",
                "External ID",
                "Organism",
                "Dataset",
                "Concentration",
                "Used Replicates",
                "Raw Optical Density",
                "Mean Relative Optical Density",
                r"Std\. Relative Optical Density",
                "Z-Factor",
            ],
        )

        error_bars = base.mark_rule().encode(
            x="Concentration:O",
            y="uerror:Q",
            y2="lerror:Q",
        )
        uerror_whiskers = base.mark_tick(size=whisker_width).encode(
            x="Concentration:O",
            y="uerror:Q",
        )
        lerror_whiskers = base.mark_tick(size=whisker_width).encode(
            x="Concentration:O",
            y="lerror:Q",
        )

        hline = base.mark_rule(strokeDash=[3, 2]).encode(
            y=alt.datum(hline_y),
            # x=[alt.value(0), alt.value(50)],
            color=alt.value("black"),
        )

        org_column = (
            alt.layer(lineplot, error_bars, uerror_whiskers, lerror_whiskers, hline)
            .facet(
                row="AsT Barcode 384",
                column="AsT Plate Subgroup",
                title=alt.Title(organism, anchor="middle"),
            )
            .resolve_axis(x="independent")
            .resolve_scale(color="independent", shape="independent")
            # .add_params(selection)
        )

        organism_columns.append(org_column)
    return alt.hconcat(*organism_columns).configure_point(size=60)


def mic_hitstogram(
    data, mic_col, title="Count Distribution of Hits over Concentration"
):
    """
    It's a Hi(t)stogram...
    Plots distribution of hits over determined MICs.
    Example: mic_distribution_overview(mic_results_long, 'MIC50 in µM')
    """
    data = data.dropna(subset=[mic_col])
    bars = (
        alt.Chart(data, title=alt.Title(title))
        .mark_bar()
        .encode(
            x=alt.X(f"{mic_col}:O"),
            y=alt.Y("count(External ID):Q"),
            xOffset="Organism:N",
            color="Organism:N",
        )
    )
    text = (
        alt.Chart(data)
        .mark_text(dx=0, dy=-5)
        .encode(
            x=alt.X(f"{mic_col}:O"),
            y=alt.Y("count(External ID):Q"),
            text=alt.Text("count(External ID):Q"),
            xOffset="Organism:N",
            color="Organism:N",
        )
    )

    return alt.layer(bars, text)
