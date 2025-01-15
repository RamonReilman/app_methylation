import panel as pn
import backend as be
import asyncio
from concurrent.futures import ThreadPoolExecutor
import nest_asyncio

nest_asyncio.apply()



pn.extension(design="material", sizing_mode="stretch_width",
             nthreads=30, loading_spinner='dots', loading_color='#2196F3')

config = be.parse_config()
main_data = be.read_data(config=config)
annotated_bed = be.load_bed_file(config=config)
gene_variation = be.read_variation_genes(config)


def create_settings():
    chr_select = pn.widgets.MultiChoice(options = main_data["chr"]
                                        .unique().to_list(), name = "chromosome:")
    group_select = pn.widgets.MultiChoice(options = main_data["group_name"]
                                          .unique().to_list(), name = "group:")


    min_range = pn.widgets.IntInput(name = "range start")
    max_range = pn.widgets.IntInput(name = "range end")


    #use_genes = pn.widgets.Checkbox(name = "use genes to search")
    all_genes = pn.widgets.MultiChoice(options = annotated_bed["gene_name"]
                                       .unique().to_list(), name = "genes:",
                                       max_items = 5)

    amount_rows_variation_file = gene_variation.with_row_index().select("index").max().item()
    top_genes_count = pn.widgets.IntInput(name=f"top x genes with highest count variation\n (max: {amount_rows_variation_file})",
                                          value=20,
                                          end=amount_rows_variation_file, start = 1)
    return pn.layout.WidgetBox("# Settings", "### Configure settings for plotting",
                               chr_select, group_select,
                               min_range, max_range,all_genes,top_genes_count)


pn.cache(max_items=10, per_session=True)
def update_df(chr_select, group_select, min_range, max_range, gene_list):
    if not chr_select and not group_select and not min_range and not max_range and not gene_list:
        return main_data

    filtered_data = main_data.clone()
    
    if gene_list:
        print("Filtering genes")
        filtered_data = be.filter_genes(gene_list, filtered_data, annotated_bed)

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


    print("Filtered data!")
    return filtered_data


async def update_tabs(plots, *args):
    tabs = pn.Tabs()
    for title, plot in plots:
        tabs.append((title, plot))
        
    for content in args:
        title, item = content
        tabs.append((title, item))

    return tabs

async def append_tabs(tabs, title, item):
    return tabs.append((title, item))

async def load_page():
    settings_box = create_settings()
    final_df = pn.bind(update_df,
                       chr_select=settings_box[2],
                       group_select=settings_box[3],
                       min_range=settings_box[4],
                       max_range=settings_box[5],
                       gene_list=settings_box[6])
    headed_gene_variation = pn.bind(be.head_variation,gene_variation ,settings_box[7])

    plots_func = pn.bind(lambda df, scatter: asyncio.run(be.plot_plots(df, scatter)), final_df, settings_box[6])

    tabs = pn.bind(lambda plots, args, : asyncio.run(update_tabs(plots, args)), plots_func, ("Gene variation" , headed_gene_variation))
    #tabs = pn.bind(lambda tabs, title, df: asyncio.run(append_tabs(tabs, title, df)), tabs, "test", headed_gene_variation)
    
    return pn.template.MaterialTemplate(
        site = "Methylatie",
        title =     "Website",
        sidebar = [settings_box],
        main = [tabs])


def run_loadpage():
    with ThreadPoolExecutor() as executor:
        task = executor.submit(asyncio.run, load_page())
        result = task.result()
    return result


def main():
    main_page = run_loadpage()
    main_page.servable()
    pn.state.onload(run_loadpage)
    

main()
