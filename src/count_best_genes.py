import backend as be
import polars as pl


def get_top_x_genes(gene_df: pl.DataFrame, df: pl.DataFrame):
    df_count_groups_list = []
    dict_start_end = {}
    for gene_info in gene_df.iter_rows():

        chr, start, end, _ = gene_info
        if chr not in dict_start_end:
            dict_start_end[chr] = []

        if [start,end] in dict_start_end[chr]:
            continue
        print(start, end)
        gene_df_filtered = be.filter_chr([chr], gene_df)
        gene_df_filtered = be.filter_ranges(start, end, gene_df_filtered)
        gene_list = pl.Series(gene_df_filtered.select("gene_name")).to_list()


        new_df = be.filter_df_gene(chr, start, end, df)
        counted_df = be.count_methylation_data(new_df)
        for gene in gene_list:
            counted_df_annot = counted_df.with_columns([
                pl.lit(gene).alias("gene")
            ])
            df_count_groups_list.append(counted_df_annot)
        dict_start_end[chr].append([start,end])

    df_count_groups = pl.concat(df_count_groups_list)
    return df_count_groups


def calc_std(counted_df: pl.DataFrame):
    std_methylations = counted_df.group_by("gene").agg(
        pl.col("n methylations").std().alias("Methylation variation")
    )
    return std_methylations


def sort_df(df: pl.DataFrame):
    return df.sort(by = "Methylation variation", descending = True)


def write_to_csv(df: pl.DataFrame):
    df.write_csv("/homes/rreilman/jaar2/kwartaal2/dashboard/app_methylation/data/test2.csv")



def main():
    config = be.parse_config()
    main = be.read_data(config)
    annotated = be.load_bed_file(config)

    counted = get_top_x_genes(annotated, main)

    std_df = calc_std(counted)
    std_df = sort_df(std_df) 
    write_to_csv(std_df)

if __name__ == "__main__":
    main()
