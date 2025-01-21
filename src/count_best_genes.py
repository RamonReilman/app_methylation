"""
count_best_genes.py
Author: Ramon Reilman
Version: 1.0
Year: BFV2

Usage:
This script will create a csv file that contains the std of the normalized positions of all methylation points

run:
python3 count_best_genes.py
"""

import backend as be
import polars as pl


def get_top_x_genes(gene_df: pl.DataFrame, df: pl.DataFrame):
    """Get all the methylation data for every gene

    This function gets all of the methylated points for every gene

    Parameters
    ----------
    gene_df: pl.DataFrame
            Contains the promoter region information for all genes
    df: pl.DataFrame
            Contains the main analysis data

    Returns
    -------
    pl.DataFrame
        Containing all of the methylated points, with the fitting gene

    """
    df_count_groups_list = []
    dict_start_end = {}

    # Loop through the gene promoter region lines
    for gene_info in gene_df.iter_rows():

        # Get the info and add it to a dict if it is not in it
        chr, start, end, _ = gene_info
        if chr not in dict_start_end:
            dict_start_end[chr] = []

        # Only continue if these ranges have not been done yet
        if [start, end] in dict_start_end[chr]:
            continue
        
        # Get all genes in this range
        gene_df_filtered = be.filter_chr([chr], gene_df)
        gene_df_filtered = be.filter_ranges(start, end, gene_df_filtered)
        gene_list = pl.Series(gene_df_filtered.select("gene_name")).to_list()

        # Filter main data for these genes
        new_df = be.filter_df_gene(chr, start, end, df)
        # Add all the gene names to this
        for gene in gene_list:
            counted_df_annot = new_df.with_columns([
                pl.lit(gene).alias("gene")
            ])
            df_count_groups_list.append(counted_df_annot)

        # Add range to dict
        # Only do this range once
        dict_start_end[chr].append([start, end])

    df_count_groups = pl.concat(df_count_groups_list)
    return df_count_groups


def calc_std(df: pl.DataFrame):
    """calculates the std

    Calculate the normalized standard deviation on position
    To get position variation in genes

    Parameters
    ----------
    df: pl.DataFrame
        df that contains all of the methylation points,
        labeled with the genes and group it belongs to

    Returns
    -------
    pl.DataFrame
        df that contains the gene name and std based on normalized position

    """
    # Calc normalized position
    lowest_pos = df.select("start").min().item()
    highest_pos = df.select("start").max().item()
    df = df.with_columns(
        ((pl.col("start") - lowest_pos) / (highest_pos - lowest_pos)).alias("norm pos"))

    # Calculate std for every gene
    std_methylations = df.group_by(["gene"]).agg(
        pl.col("norm pos").std().alias("Methylation variation")
    )

    return std_methylations


def sort_df(df: pl.DataFrame):
    """Sorts df

    Sorts df, based on std. Descending

    Parameters
    ----------
    df: pl.DataFrame
        df that contains the gene name and std based on normalized position

    Returns
    -------
    pl.DataFrame
        df that contains the gene name and std based on normalized position, sorted

    """
    return df.sort(by="Methylation variation", descending=True)


def write_to_csv(df: pl.DataFrame):
    """Writes sorted std_df to csv.

    Parameters
    ----------
    df: pl.DataFrame
        df that contains the gene name and std based on normalized position, sorted

    Returns
    -------
    None

    """
    
    df.write_csv("/homes/rreilman/jaar2/kwartaal2/dashboard/app_methylation/data/gene_variation.csv")


def main():
    """Main"""
    
    # Reads data
    config = be.parse_config()
    main = be.read_data(config)
    annotated = be.load_bed_file(config)

    # Gets all genes and methylation points
    counted = get_top_x_genes(annotated, main)

    # Calculate, sort and write
    std_df = calc_std(counted)
    std_df = sort_df(std_df)
    write_to_csv(std_df)


if __name__ == "__main__":
    main()
