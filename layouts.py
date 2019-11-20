import dash_table
import dash_html_components as html


def layout(cnr):
    table_cols = ['chromosome', 'start', 'end', 'gene', 'depth', 'log2']

    layout = html.Div(children=[
        html.H1(children='CNV Plots'),

        html.Div(children='''
            Dash: A web application framework for Python.
        '''),

        html.Div(id='datatable-interactivity-container'),

        html.Div([
            dash_table.DataTable(
                id='datatable-interactivity',
                columns=[{"name": i, "id": i, "selectable": True}
                         for i in cnr[table_cols].columns],
                data=cnr.to_dict("records"),
                editable=True,
                filter_action="native",
                sort_action="native",
                sort_mode="multi",
                row_selectable="multi",
                selected_rows=[],
                page_action="native",
                page_current=0,
                page_size=10,
            )
        ], style={'padding': 40,
                  'width': 1000,
                  'margin-left': 'auto',
                  'margin-right': 'auto'}),
    ])

    return layout
