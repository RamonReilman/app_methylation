import panel as pn
import backend as be


pn.extension(design="material", sizing_mode="stretch_width",
             nthreads=30, loading_spinner='dots', loading_color='#2196F3')

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

pn.cache(max_items=10, per_session=True)
def update_df(chr_select, group_select, min_range, max_range, gene_list):
    if not chr_select and not group_select and not min_range and not max_range and not gene_list:
        return main_data

    filtered_data = main_data.clone()
    if chr_select:
        print("Filtering chromosome")
        filtered_data = be.filter_chr(chr_select, filtered_data)

    if group_select:
        print("Filter groups")
        filtered_data = be.filter_group(group_select, filtered_data)

    if min_range or max_range:

        min_range = max(min_range, filtered_data["start"].min())
        max_range = min(max_range, filtered_data["end"].max())
        print(min_range, max_range)
        filtered_data = be.filter_ranges(min_range, max_range, filtered_data)

    if gene_list:
        print("Filtering genes")
        filtered_data = be.filter_genes(gene_list, filtered_data, annotated_bed)

    print("Filtered data!")
    return filtered_data


def load_page():
    settings_box = create_settings()
    final_df = pn.bind(update_df,
                       chr_select=settings_box[2],
                       group_select=settings_box[3],
                       min_range=settings_box[4],
                       max_range=settings_box[5],
                       gene_list=settings_box[6])


    return pn.template.MaterialTemplate(
        site = "Methylatie tijd",
        title = "Website",
        sidebar = [settings_box],
        main = [pn.param.ParamFunction(pn.bind(be.plot_plots, final_df), loading_indicator = True)])


def main():
    load_page().servable()
    pn.state.onload(load_page)

main()
