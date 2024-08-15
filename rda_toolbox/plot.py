#!/usr/bin/env python3

import altair as alt


def get_heatmap(subdf, negative_controls="Negative Control", blanks="Medium"):
    """
    - negative_controls are controls with organism + medium
    - blanks are controls with only medium (no organism and therefore no growth)
    """
    base = alt.Chart(
        subdf,
    ).encode(
        alt.X("Col_384:O").axis(labelAngle=0, orient="top").title(None),
        alt.Y("Row_384:O").title(None),
        tooltip=list(subdf.columns),
    )
    blank_mean = subdf[subdf["ID"] == blanks]["Raw Optical Density"].mean()
    negative_mean = subdf[subdf["ID"] == negative_controls][
        "Raw Optical Density"
    ].mean()
    heatmap = base.mark_rect().encode(
        alt.Color("Raw Optical Density:Q")
        .title("Optical Density")
        .scale(domain=[blank_mean, negative_mean]),
    )
    text = base.mark_text(baseline="middle", align="center", fontSize=10).encode(
        alt.Text("Raw Optical Density:Q", format=".1f"),
        color=alt.condition(
            alt.datum["Raw Optical Density"]
            < max(subdf[subdf["ID"] == blanks]["Raw Optical Density"]) / 2,
            alt.value("black"),
            alt.value("white"),
        ),
    )
    return alt.layer(heatmap, text)


def plateheatmaps(df):
    """
    Plots heatmaps of the plates from df in a gridlike manner.
    Exclude unwanted plates, for example Blanks from the df like so
        `df[df["Organism"] != "Blank"]`
    before plotting, otherwise it will appear as an extra plate.
    """
    plots = []
    for _, _organism_df in df.groupby("Organism"):
        plots.append(
            get_heatmap(_organism_df)
            .facet(
                # row=alt.Row("Barcode 384-AST:N", sort=ast_keys),
                row=alt.Row("Barcode:N"),
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
    text = base.mark_text(
        baseline="middle",
        align="center",
        fontSize=10
    ).encode(
        alt.Text("Raw Optical Density:Q", format=".1f"),
    )
    return alt.layer(heatmap, text)
