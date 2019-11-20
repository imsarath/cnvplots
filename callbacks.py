import pandas as pd
import dash_core_components as dcc
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate


from cnvplots import CNVPlot


def update_callback(dashapp, cnr, seg):
    @dashapp.callback(
        Output('datatable-interactivity-container', "children"),
        [Input('datatable-interactivity', "derived_virtual_data"),
         Input('datatable-interactivity', "derived_virtual_selected_rows")])
    def update_graphs(rows, derived_virtual_selected_rows):
        highlight = True

        cnr_df = cnr if rows is None else pd.DataFrame(rows)

        if derived_virtual_selected_rows is None:
            derived_virtual_selected_rows = []
            highlight = False

        if cnr_df.empty:
            raise PreventUpdate
        else:
            return dcc.Graph(id='example-graph',
                             figure=CNVPlot(cnr_df,
                                            seg,
                                            derived_virtual_selected_rows,
                                            highlight=highlight))
