import os
import hvplot.polars
import polars as pl
import panel as pn
import re
import pandas as pd
import configparser

@pn.cache
def parse_config()-> configparser:
    config = configparser.ConfigParser()
    config.read("data/config.ini")
    return config


def process_groups(config: configparser) -> pl.DataFrame:
    path = config.get("PATHS", "group_data")

    barcodes_names: pl.dataframe = pl.read_csv(path)

    barcodes_names = barcodes_names.with_columns(controle_n = pl.int_range(pl.len()).over(" description")+1)
    barcodes_names = barcodes_names.with_columns(group_and_n = pl.concat_str([pl.col(' description'), pl.col("controle_n")]))

    barcodes_names = barcodes_names.with_columns(pl.col(pl.Utf8).str.strip_chars()).drop("controle_n")


    return barcodes_names


@pn.cache
def read_data(config: configparser) -> pl.DataFrame:
    path = config.get("PATHS", "data_folder")
    barcodes_names = process_groups(config)
    resulting_df: pd.DataFrame = pl.DataFrame(
        {"chr":[],
         "start":[],
         "end":[],
         "frac":[],
         "valid":[],
         "group_name":[]}
    )
    files: list[str] = os.listdir(path)

    for file in files:
        if os.path.isfile(f"{path}/{file}") and file.endswith("methylatie_ALL.csv"):
            temp_df: pd.DataFrame = pd.read_csv(f"{path}/{file}", sep="\t")
            temp_df: pl.DataFrame = pl.from_pandas(temp_df)
            barcode_num: list[int] = re.findall(r"\d+", file)

            name_group: str = barcodes_names.filter(pl.col("barcode").cast(pl.String) == barcode_num[0]).select("group_and_n")
            temp_df: pl.DataFrame = temp_df.with_columns(pl.lit(name_group).alias("group_name"))
            resulting_df = pl.concat([temp_df, resulting_df])


    return resulting_df


def count_methylation_data(df: pl.DataFrame) -> pl.DataFrame:
    all_groups = pl.DataFrame({"group_name": df["group_name"].unique()})

    return (df
    .select(["group_name", "frac"])
    .group_by("group_name")
    .agg([pl.len().alias("n methylations")])
    .join(all_groups, on="group_name", how="full")
    .with_columns(pl.col("group_name").fill_null(pl.col("group_name_right")))
    .drop("group_name_right")
    .fill_null(0)
)


def get_gene_info(annotated_bed: pl.DataFrame, genes: list[str]) -> pl.DataFrame:
    print("Getting gene info")
    df_wanted = (annotated_bed
                 .filter(pl.col("gene_name").is_in(genes)))
    return df_wanted


def filter_genes(gene_list: list[str], df: pl.DataFrame, annotated_bed: pl.DataFrame) -> pl.DataFrame:
    df_wanted = get_gene_info(annotated_bed, gene_list)
    final_subsetted_df: pd.DataFrame = pl.DataFrame(
        {"chr":[],
         "start":[],
         "end":[],
         "frac":[],
         "valid":[],
         "group_name":[]})

    for row in df_wanted.iter_rows():
        (chromosome, promoter_start, promoter_end, gene_name) = row
        subsetted_df = df.filter(
            (pl.col("chr") == chromosome) &
            (pl.col("start") >= promoter_start) &
            (pl.col("end") <= promoter_end))

        final_subsetted_df = pl.concat([subsetted_df, final_subsetted_df])


    return final_subsetted_df


def plot_barchart(df: pl.DataFrame) -> hvplot.plot:
    df = count_methylation_data(df)
    barplot = df.hvplot.bar(x = "group_name", y="n methylations",
                            color = "group_name", cmap = "Category10",
                            width = 750, rot=20, height = 400,
                            title = "Number of methylations for every group",
                            xlabel = "Group name")
    return barplot

def plot_scatter(df: pl.DataFrame) -> hvplot.plot:
    return df.hvplot.scatter(x = "start", y="chr", by="group_name", width = 750,
                             height = 400,
                             title = "Methylated DNA points",
                             xlabel = "Start positon of methylation",
                             ylabel = "Chromosome")

def plot_density(df: pl.DataFrame) -> hvplot.plot:
    df = df.select(["group_name", "start"])
    return df.hvplot.kde(by="group_name", width = 750, height = 400,
                         title = "Density of methylation positions",
                         xlabel = "genomic positions")

def plot_plots(df: pl.DataFrame) -> pn.layout.Row[hvplot.plot]:
    barplot = plot_barchart(df)
    #scatter = plot_scatter(df)
    density = plot_density(df)
    return pn.layout.Row(barplot, density)
    


def load_bed_file(config: configparser) -> pl.DataFrame:
    return pl.read_csv(config["PATHS"]["annotated_bed"], separator=",")


def filter_chr(chr_list: list[str], df: pl.DataFrame) -> pl.DataFrame:
    return (df.filter(pl.col("chr").is_in(chr_list)))


def filter_group(group_list: list[str], df: pl.DataFrame) -> pl.DataFrame:
    return (df.filter(pl.col("group_name").is_in(group_list)))


def filter_ranges(min_range: int, max_range: int, df: pl.DataFrame) -> pl.DataFrame:
    return (df.filter((pl.col("start") >= min_range) &
                      (pl.col("end") <= max_range)))
