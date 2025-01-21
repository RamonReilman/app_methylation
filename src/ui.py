"""
Ui.py
Author: Ramon Reilman
Version: 1.0
Year: BFV2

Usage:
This script serves as the main way to run the methylation webiste
This script creates the ui, and calls functions from backend.py

run:
panel serve src/ui.py --port=5100 --allow-websocket-origin='*' --num-procs 60 --reuse-sessions --global-loading-spinner
"""
import asyncio
from concurrent.futures import ThreadPoolExecutor
import panel as pn
import nest_asyncio
import backend as be

nest_asyncio.apply()
pn.extension(design="material", sizing_mode="stretch_width",
             nthreads=30, loading_spinner='dots', loading_color='#2196F3')
pn.param.ParamMethod.loading_indicator = True


config = be.parse_config()
main_data = be.read_data(config=config)
annotated_bed = be.load_bed_file(config=config)
gene_variation = be.read_variation_genes(config)


def create_settings():
    """Creates settings widgets

    This functions creates all of the settings for
    the wdigetbox on the site of the webpage

    Parameters
    ----------
    None

    Returns
    -------
    pn.layout.Widgetbox
                        A panel object that creates a box on the side.
                        This box contains all the settings made in this function.

    """
    chr_select = pn.widgets.MultiChoice(options=main_data["chr"]
                                        .unique().to_list(), name="chromosome:")
    group_select = pn.widgets.MultiChoice(options=main_data["group_name"]
                                          .unique().to_list(), name="group:")
    highest_start = main_data.select("start").max().item()
    min_range = pn.widgets.IntInput(name="range start",
                                    start=0,
                                    end=highest_start)
    max_range = pn.widgets.IntInput(name="range end"
                                    ,start=0
                                    ,end=highest_start)

    all_genes = pn.widgets.MultiChoice(options=annotated_bed["gene_name"]
                                       .unique().to_list(), name="genes:",
                                       max_items=5,
                                       option_limit=10)

    amount_rows_variation_file = gene_variation.with_row_index().select("index").max().item()
    name_top_genes_count = f"top x genes to display\n ({amount_rows_variation_file})"
    top_genes_count = pn.widgets.IntInput(name=name_top_genes_count,
                                          value=20,
                                          end=amount_rows_variation_file, start=1)

    submit = pn.widgets.Button(name='Filter data!', button_type='primary')
    return pn.layout.WidgetBox("# Settings", "### Configure settings for plotting",
                               chr_select, group_select,
                               min_range, max_range, all_genes, top_genes_count,
                               submit)


@pn.cache(max_items=10, per_session=True)
def update_df(chr_select, group_select, min_range, max_range, gene_list):
    """Filters df based on user input

    This function will take values set by the user in the frontend
    and will filter the main data based on those values.

    Parameters
    ----------
    chr_select : list
                The value of what chromosomes want to be seen
                comes from widgetbox
    group_select : list
                The value of what group the user wants to see
                comes from widgetbox
    min_range : int
                The lowest start value the user wants to see
                comes from widgetbox
    max_rane : int
                The highest start value the user wants to see
                comes from widgetbox
    gene_list : list
                The value of what genes the user wants to see
                comes from widgetbox

    Returns
    -------
    filtered_data : pl.Dataframe
                A polars dataframe filter on what the user wants to see.

    """
    filtered_data = main_data.clone()

    if gene_list:
        print("Filtering genes")
        filtered_data = be.filter_genes(
            gene_list, filtered_data, annotated_bed)

    if chr_select:
        print("Filtering chromosome")
        filtered_data = be.filter_chr(chr_select, filtered_data)

    if group_select:
        print("Filter groups")
        filtered_data = be.filter_group(group_select, filtered_data)

    if min_range or max_range:
        # Checks to see if user filters between the lowest and highest value
        min_range = max(min_range, filtered_data["start"].min())
        max_range = min(max_range, filtered_data["end"].max())
        print(min_range, max_range)
        filtered_data = be.filter_ranges(min_range, max_range, filtered_data)

    print("Filtered data!")
    return filtered_data


@pn.cache
def create_tabs(plots, *args):
    """Creates all of the tabs

    The webpage works by splitting plots and tables in tabs
    This funtions will create all tabs

    Parameters
    ----------
    Plots : list
            List that contains the plots and title of the tabs
    args : list
            Can be anything i want to add to the tabs

    Returns
    -------
    pn.Tabs
            Panel tabs object that contains the names of the tabs
            together with the page's content

    """
    print("Updating tabs")
    tabs = pn.Tabs()

    # Puts plots in the tabs
    if isinstance(plots, list):
        for title, plot in plots:
            tabs.append((title, plot))
    else:
        tabs.append(("Empty Dataframe", be.pn.pane.Markdown(
            "# No methylation found for these filters!")))
    # Puts args content in tabs, if args is used
    if args:
        for content in args:
            title, item = content
            tabs.append((title, item))

    return tabs


async def submit_button(button, settings_box):
    """Filters and plots when filter button is pressed

    This function is called when the filter button is pressed
    It wil filter the main data, plot the plots with the filtered data
    and return it in a tabs object

    Parameters
    ----------
    Plots : list
            List that contains the plots and title of the tabs
    args : list
            Can be anything i want to add to the tabs

    Returns
    -------
    pn.Tabs
            Panel tabs object that contains the names of the tabs
            together with the page's content

    """
    # Clones data
    temp_data = main_data.clone()
    if button:

        # Filter data
        temp_data = update_df(chr_select=settings_box[2].value,
                              group_select=settings_box[3].value,
                              min_range=settings_box[4].value,
                              max_range=settings_box[5].value,
                              gene_list=settings_box[6].value)
    headed_gene_variation = be.head_variation(
        gene_variation, settings_box[7].value)
    # Plot data
    plots = asyncio.run(be.plot_plots(temp_data, settings_box[6].value))
    if "gene_name" in temp_data.columns:
        return create_tabs(plots,
                            ("Gene Variation", headed_gene_variation),
                            ("Filtered data", pn.pane.DataFrame(temp_data.select(["chr", "start", "end", "group_name", "gene_name"]).to_pandas())))
    return create_tabs(plots,
                            ("Gene Variation", headed_gene_variation))



async def load_page():
    """Loads page

    This function will load all widgets for the page
    Will also ready the contents of the page

    Parameters
    ----------
    None

    Returns
    -------
    pn.template.MaterialTemplate
            A template that contains all of the content of the website in tabs
            has a sidebar that contains all of the settings

    """
    settings_box = create_settings()

    # Binds button to submit button function
    tabs = pn.bind(lambda button, settings_box: asyncio.run(submit_button(button,
                                                                          settings_box)
                                                            ), settings_box[8], settings_box)

    return pn.template.MaterialTemplate(
        site="Methylatie",
        title="Website",
        sidebar=[settings_box],
        main=[tabs])


def run_loadpage():
    """Creates a loop to run the async load_page function in

    This function will create a event loop to run the load_page funtion in
    This is necesarry because it will conflict with panel

    Parameters
    ----------
    None

    Returns
    -------
    result
            This will contain the pn.template object created in load_page()

    """
    with ThreadPoolExecutor() as executor:
        task = executor.submit(asyncio.run, load_page())
        result = task.result()

    return result


def main():
    """
    One main to rule them all
    """
    main_page = run_loadpage()
    main_page.servable()
    pn.state.onload(run_loadpage)


if __name__ == "__main__":
    pn.serve(main)
else:
    main()
