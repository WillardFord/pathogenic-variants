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
    pathogenic_file = "data/ACMG SF Annotations Download - ACMG SF v3.3 List.tsv"
    return (pathogenic_file,)


@app.cell
def _(pathogenic_file, pl):
    pathogenic_vars = pl.read_csv(pathogenic_file, separator="\t")
    pathogenic_vars["Variants to report"].unique()
    return (pathogenic_vars,)


@app.cell
def _(pathogenic_vars):
    # Save gene names to file with colon suffix, one per line, no header
    gene_names = pathogenic_vars["Gene"].unique().sort()
    gene_names_with_colon = gene_names + ":"

    with open("data/pathogenic_genes.txt", "w") as f:
        for gene in gene_names_with_colon:
            f.write(f"{gene}\n")

    print(f"Saved {len(gene_names)} unique gene names to data/pathogenic_genes.txt")
    return


@app.cell
def _():
    annotations_file = "data/clinvar.vcf.gz"
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
