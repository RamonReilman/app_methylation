import panel as pn
import backend as be


pn.extension(design="material", sizing_mode="stretch_width")

config = be.parse_config()
main_data = be.read_data(config=config)
annotated_bed = be.load_bed_file(config=config)


def create_settings():
    chr_select = pn.widgets.MultiChoice(options = main_data["chr"]
                                        .unique().to_list(), name = "chromosome:")
    group_select = pn.widgets.MultiChoice(options = main_data["group_name"]
                                          .unique().to_list(), name = "group:")


    min_range = pn.widgets.IntInput(name = "range start")
    max_range = pn.widgets.IntInput(name = "range end")


    #use_genes = pn.widgets.Checkbox(name = "use genes to search")
    all_genes = pn.widgets.MultiChoice(options = annotated_bed["gene_name"]
                                       .unique().to_list(), name = "genes:")

    return pn.layout.WidgetBox("# Settings", "### Configure settings for plotting",
                               chr_select, group_select,
                               min_range, max_range,all_genes)


def update_df(chr_select, group_select, min_range, max_range, all_genes):
    filtered_data = main_data.clone()
    if chr_select:
        print("Filtering chromosome")
        filtered_data = be.filter_chr(chr_select, filtered_data)

    if group_select:
        print("Filter groups")
        filtered_data = be.filter_group(group_select, filtered_data)

    if min_range and max_range:
        print(min_range, max_range)
        filtered_data = be.filter_ranges(min_range, max_range, filtered_data)

    if all_genes:
        print("Filtering genes")
        filtered_data = be.filter_genes(all_genes, filtered_data, annotated_bed)

    print("Filtered data!")
    return filtered_data


def load_page():
    settings_box = create_settings()
    final_df = pn.bind(update_df,
                       chr_select=settings_box[2],
                       group_select=settings_box[3],
                       min_range=settings_box[4],
                       max_range=settings_box[5],
                       all_genes=settings_box[6])

    plots = pn.bind(be.plot_plots, final_df)
    print("Plotted!")


    return pn.template.MaterialTemplate(
        site = "Methylatie tijd",
        title = "gangang",
        sidebar = [settings_box],
        main = [plots])


def main():
    load_page().servable()


main()
