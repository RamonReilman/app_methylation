import panel as pn


pn.extension()

def create_main_settings_col():
    which_information = pn.widgets.Select(name="Plot choice", 
                                          options = ["Barchart", "Loliplot", "Sequence"],
                                          value = "Barchart")

    which_groups = pn.widgets.MultiChoice(name = "Select groups",
                                          options = ["All", "barcode 11", "barcode 12", "barcode 13"],)

    wants_groups = pn.widgets.Select(name = "Genes grouped by function?",
                                     options = ["Group genes", "Just genes"], value = "Just genes")

    which_genes = pn.widgets.MultiChoice(name = "Select genes",
                                         options = ["TP53", "CDK1", "PICC", "BRCA2"])
    return pn.WidgetBox("# Settings", which_information, 
                        which_groups,wants_groups ,which_genes)

def create_ui():
    create_main_settings_col().servable()

create_ui()