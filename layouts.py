import dash_table
import dash_core_components as dcc
import dash_html_components as html


def layout(cnr):
    table_cols = ['chromosome', 'start', 'end', 'gene', 'depth', 'log2']

    layout = html.Div(
            id='main-page',
            children=[
                html.Div(
                    id='app-page-header',
                    children=[
                        html.H1(children='CNV Plots')
                    ],
                    style={
                        'text-align': 'center',
                    }
                ),
                html.Div(
                    id='app-page-content',
                    children=[
                        html.Div(dcc.Markdown('''
                            Developed an application for interactive Genome-wide copy number
                            plot from cnvkit. Dash is a python framework which is being used 
                            to build this application. This scatter plot has bin-level log2 
                            coverages and segmentation calls(lines) together.

                            Dash-datatable has been used to create interactive table which
                            is shown below. It has multiple options, we can select and
                            filter values from the table and the plot will be altered based on
                            filtered data.

                            ### How To Use Filters:
                             - Each column has filter text-box. we can use relational
                             operators to filter numeric values ( <, >, =). we have to
                             type values and press enter to trigger the filter action.

                             - Chromosome, start, end, depth and log2 columns have numeric
                             values. i.e: (for log2 column - >0 rows greater than o or <0
                             rows less than 0)

                             - Gene column values are strings, so simply we can use gene
                             names to filter values. i.e : PTEN, Plot and Table will be
                             showing only rows with PTEN.
                             
                             - If we want to highlight the regions which you are interested in,
                             we can select those rows and those regions will be highlighted on
                             the plot.
                             - The plot will be updated on the flow while filtering the table.
                        '''), style={'width': 750,
                                     'text-align': 'justify',
                                     'margin-right': 'auto',
                                     'margin-left': 'auto'}),

                        html.Div(id='datatable-interactivity-container',
                                 style={
                                    'width': 1000,
                                    'margin-right': 'auto',
                                    'margin-left': 'auto'
                                 }),

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
                                  'margin-right': 'auto'})
                    ]
                )
            ])

    return layout
