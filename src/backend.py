"""
Backend.py
Author: Ramon Reilman
Version: 1.0
Year: BFV2

Usage:
Serves as the main backend of the website.
Will contain functions used in ui.py
Not meant to be used without ui.py

run:
panel serve src/ui.py --port=5100 --allow-websocket-origin='*' 
--num-procs 60 --reuse-sessions --global-loading-spinner
"""

import os
import re
import configparser
import asyncio
import hvplot.polars
import polars as pl
import panel as pn
import pandas as pd


@pn.cache
def parse_config() -> configparser:
    """Reads config file

    This function will read a config file
    This file contains all of the paths to the data
    ./data/config.ini

    Parameters
    ----------
    None

    Returns
    -------
    configparser.ConfigParser
                        ConfigParser object that contains config.ini information

    """
    config = configparser.ConfigParser()
    config.read("data/config.ini")
    return config


def process_groups(config: configparser) -> pl.DataFrame:
    """Processes group.csv information

    Will get the group name with barcode information
    Will also seperate duplicate group names with an index

    Parameters
    ----------
    config : ConfigParser
            Contains the paths to needed files

    Returns
    -------
    pl.DataFrame
            pl.DataFrame object that contains the group_name and the barcode

    """
    # Get file for group data
    path = config.get("PATHS", "group_data")

    try:
        barcodes_names: pl.dataframe = pl.read_csv(path)

        # Give duplicate names an index to differentiate them
        barcodes_names = barcodes_names.with_columns(
            controle_n=pl.int_range(pl.len()).over(" description") + 1)
        barcodes_names = barcodes_names.with_columns(
            group_and_n=pl.concat_str([pl.col(' description'), pl.col("controle_n")]))

        # Strip the strings of ecess of characters
        barcodes_names = barcodes_names.with_columns(pl.col(pl.Utf8).str.strip_chars()
                                                     ).drop("controle_n")

    # Raising errors
    except pl.exceptions.ColumnNotFoundError as error:
        print(
            f"Header in {path} file not correct, should be: {' description'}\n {error}")
        raise

    except FileNotFoundError:
        print(f"file: {path} not found or incorrect permissions")
        raise

    return barcodes_names


@pn.cache(max_items=10, per_session=True)
def read_data(config: configparser) -> pl.DataFrame:
    """Reads the main analysis data

    This function will read all of the analysis data and process it.

    Parameters
    ----------
    config : ConfigParser
            Contains the paths to needed files

    Returns
    -------
    pl.Dataframe
            pl.dataframe that contains all of the methylation data
            Extra column is added that will link the data to the group it belongs to

    """
    # Read paths and get the group information
    path = config.get("PATHS", "data_folder")
    barcodes_names = process_groups(config)
    resulting_df: pd.DataFrame = pl.DataFrame(
        {"chr": [],
         "start": [],
         "end": [],
         "frac": [],
         "valid": [],
         "group_name": []})

    # Get the analysis files
    files: list[str] = os.listdir(path)

    for file in files:
        # Check if file exists
        if os.path.isfile(
                f"{path}/{file}") and file.endswith("methylatie_ALL.csv"):

            # Read it and make it a polars df
            temp_df: pd.DataFrame = pd.read_csv(f"{path}/{file}", sep="\t")
            temp_df: pl.DataFrame = pl.from_pandas(temp_df)

            # Remove excess characters from chr
            temp_df = temp_df.with_columns(
                pl.col("chr").str.split("_").list.get(0))

            # Figure out with barcode the file has
            barcode_num: list[int] = re.findall(r"\d+", file)
            name_group: str = barcodes_names.filter(pl.col("barcode").cast(
                pl.String) == barcode_num[0]).select("group_and_n")

            # Name the group accoring to the barcode number
            temp_df: pl.DataFrame = temp_df.with_columns(
                pl.lit(name_group).alias("group_name"))
            resulting_df = pl.concat([temp_df, resulting_df])

    return resulting_df


def count_methylation_data(df: pl.DataFrame) -> pl.DataFrame:
    """Counts amount of methylation points in df

    This function will count all of the methylation points in a df
    it will be split on the groups in the df

    Parameters
    ----------
    df : pl.DataFrame
        Dataframe that contains the methylation data

    Returns
    -------
    pl.Dataframe
            pl.DataFrame that will contain the amount of methylated spots for every group

    """
    # Get All groups in the df
    all_groups = pl.DataFrame({"group_name": df["group_name"].unique()})

    # Count it
    return (df
            .select(["group_name", "frac"])
            .group_by("group_name")
            .agg([pl.len().alias("n methylations")])
            .join(all_groups, on="group_name", how="full")
            .with_columns(pl.col("group_name").fill_null(pl.col("group_name_right")))
            .drop("group_name_right")
            .fill_null(0))


def get_gene_info(annotated_bed: pl.DataFrame, genes: list[str]) -> pl.DataFrame:
    """Will get gene information

    This function will read all of the analysis data and process it.

    Parameters
    ----------
    annotated_bed : pl.DataFrame
            Contains promoter sites of (mostly) all human genes
    genes : list
            A list containing genes the user wants to see, comes from frontend

    Returns
    -------
    pl.Dataframe
            pl.DataFrame that contains promoter sites for the genes the user specified

    """
    print("Getting gene info")
    df_wanted = (annotated_bed
                 .filter(pl.col("gene_name").is_in(genes)))

    return df_wanted


def filter_df_gene(chromosome: str, start: int, end: int, df: pl.DataFrame) -> pl.DataFrame:
    """Filters the main data on wanted genes

    This function will filter the main data on the genes.
    These genes are specified in the frontend, by the user.

    Parameters
    ----------
    chromosome : str
        A chromosome name to filter on
    start : int
        a starting range to filter on
    end : int
        end of the range to filter on
    df : pl.DataFrame
        main analysis data

    Returns
    -------
    pl.Dataframe
            main analysis data, filtered on the filters specified above

    """
    return df.filter(
        (pl.col("chr") == chromosome) &
        (pl.col("start") >= start) &
        (pl.col("end") <= end))


def filter_genes(gene_list: list[str], df: pl.DataFrame, annotated_bed: pl.DataFrame) -> pl.DataFrame:
    """Gets list with genes and filters main df on it

    This function will filter the main analysis data based on a list of given genes.

    Parameters
    ----------
    gene_list : list
            A list containing genes the user wants to see, comes from frontend
    df : pl.DataFrame
            Main analysis data
    annotated_bed : pl.DataFrame
            Contains promoter sites of (mostly) all human genes

    Returns
    -------
    pl.Dataframe
            pl.dataframe that contains all of the methylation data
            Extra column is added that will link the data to the group it belongs to

    """
    # Get a df that contains gene promoter regions
    df_wanted = get_gene_info(annotated_bed, gene_list)
    final_subsetted_df: pd.DataFrame = pl.DataFrame(
        {"chr": [],
         "start": [],
         "end": [],
         "frac": [],
         "valid": [],
         "group_name": [],
         "gene_name": {}})

    # Loop through genes
    for row in df_wanted.iter_rows():
        # Extract promoter regions
        (chromosome, promoter_start, promoter_end, gene) = row
        # Filter main data based on thos reagions
        subsetted_df = filter_df_gene(chromosome, promoter_start, promoter_end, df)

        # Add gene name to the resulting df
        subsetted_df = subsetted_df.with_columns([
            pl.lit(gene).alias("gene_name")
        ])
        final_subsetted_df = pl.concat([subsetted_df, final_subsetted_df])

    return final_subsetted_df


def loading_indicator(label: str) -> pn.indicators.LoadingSpinner:
    """Create loading indicator

    Used when website should display a loading filter

    Parameters
    ----------
    label : str
            contains the label the loading indicator should show

    Returns
    -------
    pn.indicators.LoadingSpinner
            LoadingSpinner

    """
    return pn.indicators.LoadingSpinner(
        value=True, name=label, size=25, align="center"
    )


async def plot_barchart(df: pl.DataFrame) -> hvplot.plot:
    """Plots a barchart

    This function will plot a hvplot barchart

    Parameters
    ----------
    df : pl.DataFrame
        main analysis dataframe

    Returns
    -------
    hvplot.barchart

    """
    # Get count data and plot
    df = count_methylation_data(df)
    barplot = df.hvplot.bar(x="group_name", y="n methylations",
                            color="group_name", cmap="Category10",
                            width=1125, rot=20, height=600,
                            title="Number of methylations for every group",
                            xlabel="Group name")
    return barplot


async def plot_scatter(df: pl.DataFrame) -> hvplot.plot:
    """Plots a scatter

    This function will plot a hvplot scatter

    Parameters
    ----------
    df : pl.DataFrame
        main analysis dataframe

    Returns
    -------
    hvplot.scatter

    """
    return df.hvplot.scatter(x="start", y="chr", by="group_name", width=1125,
                             dynamic=False,
                             alpha=0.2,
                             height=600,
                             title="Methylated DNA points",
                             xlabel="Start positon of methylation",
                             ylabel="Chromosome",)


async def plot_density(df: pl.DataFrame) -> hvplot.plot:
    """Plots a density

    This function will plot a hvplot density

    Parameters
    ----------
    df : pl.DataFrame
        main analysis dataframe

    Returns
    -------
    hvplot.density

    """
    df = df.select(["group_name", "start"])

    return df.hvplot.kde(by="group_name", width=1125, height=600,
                         title="Density of methylation positions",
                         xlabel="genomic positions",
                         subplots=True)


@pn.cache(max_items=10, per_session=True)
async def plot_plots(df: pl.DataFrame, want_scatter: list[str]) -> list[tuple]:
    """Plots all wanted plots

    This function will plot all wanted plots

    Parameters
    ----------
    df : pl.DataFrame
        main analysis dataframe
    want_scatter : list
        Used a check if the scatterplot should be made

    Returns
    -------
    tuple
        This tuple contains for everyplot a list with title of the page and the plot

    """
    if df.is_empty():
        return loading_indicator("Data missing!")

    # Start plotting
    barplot_task = asyncio.create_task(plot_barchart(df))
    density_task = asyncio.create_task(plot_density(df))

    scatter_task = asyncio.create_task(
        plot_scatter(df)) if want_scatter else None

    # Await tasks
    barplot = await barplot_task
    density = await density_task

    scatter = await scatter_task if scatter_task else pn.pane.Markdown("""
                                                                       # To get a scatter plot please select 1 or more genes
                                                                       Since the scatter plot works extremely slow, you have to select a section of genes to view the start points.
                                                                       """)

    print("Plotted!")
    return [("Barplot", barplot),
            ("Density plot",density),
            ("Scatter plot", scatter)]


def load_bed_file(config: configparser) -> pl.DataFrame:
    """Loads a bed file

    Loads the annotated bed file

    Parameters
    ----------
    config : ConfigParser
            Contains the paths to needed files

    Returns
    -------
    pl.DataFrame
        Contains the gene promoter information

    """
    path = config["PATHS"]["annotated_bed"]
    try:
        return pl.read_csv(path, separator=",")
    except FileNotFoundError:
        print(f"file: {path} not found or incorrect permissions")
    return None


def filter_chr(chr_list: list[str], df: pl.DataFrame) -> pl.DataFrame:
    """Filters main df on chr

    Filters the main analysis data based on given chr values

    Parameters
    ----------
    chr_list : list
            List containing wanted chromosomes
    df : pl.DataFrame
            Main analysis data

    Returns
    -------
    pl.DataFrame
        main analysis data filtered on wanted chromosomes

    """
    return (df.filter(pl.col("chr").is_in(chr_list)))


def filter_group(group_list: list[str], df: pl.DataFrame) -> pl.DataFrame:
    """Filters main df on group

    Filters the main analysis data based on given group values

    Parameters
    ----------
    chr_list : list
            List containing wanted groups
    df : pl.DataFrame
            Main analysis data

    Returns
    -------
    pl.DataFrame
        main analysis data filtered on wanted groups

    """
    return (df.filter(pl.col("group_name").is_in(group_list)))


def filter_ranges(min_range: int, max_range: int, df: pl.DataFrame) -> pl.DataFrame:
    """Filters main df on ranges

    Filters the main analysis data based on given a given start-end range

    Parameters
    ----------
    min_range: int
            min range to filter on
    max_range: int
            max range to filter on
    df : pl.DataFrame
            Main analysis data

    Returns
    -------
    pl.DataFrame
        main analysis data filtered on range

    """
    return (df.filter((pl.col("start") >= min_range) &
                      (pl.col("end") <= max_range)))


def head_variation(df, n_amount):
    """returns n gene variation df rows

    Parameters
    ----------
    n_amount : int
            amount of rows to return
    df : pl.DataFrame
            gene variation dataframe

    Returns
    -------
    pl.DataFrame
        gene variation with n_amount rows returned

    """
    if isinstance(df, pl.DataFrame):
        df = df.head(n=n_amount)
        return pn.pane.DataFrame(df.to_pandas())
    return df


def read_variation_genes(config: configparser):
    """Reads gene variation file

    Parameters
    ----------
    config : ConfigParser
            Contains the paths to needed files

    Returns
    -------
    pl.DataFrame
        gene variation dataframe (/data/gene_variation.csv)

    """
    # Check if file exists, if not return markdown
    path = config.get("PATHS", "top_genes")
    if os.path.isfile(path):
        return pl.read_csv(path).drop_nulls()
    return pn.pane.Markdown("""
                            # No gene variation file found
                            ## Please run the count_best_genes.py script!
                            ## Or ask web-manager (Ramon) to do so!
                            """)
