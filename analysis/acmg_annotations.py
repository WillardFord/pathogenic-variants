import marimo

__generated_with = "0.15.5"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    return (pl,)


@app.cell
def _():
    # Data from:
    # https://docs.google.com/spreadsheets/d/1ecSUe0bJ-2gAMXBsrf5U1wMEDOmiIjtR_A5388uj7HY/edit?gid=1946972754#gid=1946972754
    annotations_file = "data/ACMG SF Annotations Download - ACMG SF v3.3 List.tsv"
    return (annotations_file,)


@app.cell
def _(annotations_file, pl):
    annotations = pl.read_csv(annotations_file, separator="\t")
    annotations["Variants to report"].unique()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
