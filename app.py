#!/usr/env/env python

# import libraries
import pandas as pd
import dash

from callbacks import update_callback
from layouts import layout

cnr = pd.read_csv('data/sample_cnvkit.cnr', sep='\t')
seg = pd.read_csv('data/sample_cnvkit.cns', sep='\t')

app = dash.Dash(__name__)
app.layout = layout(cnr)

update_callback(app, cnr, seg)

if __name__ == '__main__':
    app.run_server(debug=True)
