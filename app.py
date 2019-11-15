#!/usr/env/env python

# import libraries
import numpy as np
import pandas as pd

import plotly.graph_objs as go
import dash
import dash_table
import dash_bio as dashbio
import dash_core_components as dcc
import dash_html_components as html


def _get_hover_text(df, genename=None, position=None, log2=None):
    """
    CNV plot - hover text
    """
    hover_text = ''

    if genename is not None and genename in df.columns:
        hover_text = hover_text + '<br>GENE: ' \
                                + df[genename].astype(str)

    if position is not None and position in df.columns:
        hover_text = hover_text + '<br>POSITION: ' \
                                + df[position].astype(str)

    if log2 is not None and log2 in df.columns:
        hover_text = hover_text + '<br>LOG2: ' \
                                + df[log2].astype(str)

    return hover_text


def by_chromosome(segments):

    for chrom, subtable in segments.groupby("chromosome", sort=False):
        yield chrom, pd.DataFrame(subtable)


def CNVPlot(cnr,
            seg,
            chrm="chromosome",
            log2="log2",
            gene="gene",
            start="start",
            end="end",
            title="Copy Number Aberation",
            showgrid=False,
            xlabel=None,
            ylabel="Copy Ratio(log2)",
            point_size=5,
            showlegend=True,
            col=None):

    cnv = _CNVPlot(cnr,
                   seg,
                   chrm=chrm,
                   log2=log2,
                   start=start,
                   end=end,
                   gene=gene)

    return cnv.figure(
        title=title,
        showgrid=showgrid,
        xlabel=None,
        ylabel="Copy Ratio(log2)",
        point_size=5,
        showlegend=True,
        col=None)


class _CNVPlot():
    """docstring for CNVPlot"""
    def __init__(self,
                 cnratio,
                 segements,
                 chrm="chromosome",
                 log2="log2",
                 gene="gene",
                 start="start",
                 end="end"):

        self.data = pd.DataFrame(cnratio)
        self.segments = pd.DataFrame(segements)
        self.mp = 'MAP_POS'
        self.xlabel = ""
        self.ticks = []
        self.ticklabels = []
        self.nChr = len(cnratio[chrm].unique())
        self.chrName = chrm
        self.geneName = gene
        self.log2 = log2
        self.start = start
        self.end = end
        self.index = 'INDEX'
        self.pos = 'POSITION'
        self.start_pos = 'START_POS'
        self.end_pos = 'END_POS'

        # Fixes the bug where one chromosome is missing by adding a sequential
        # index column.
        self.data[self.mp] = 0.5 * (self.data[start] + self.data[end])

        idx = 0
        for i in self.data[chrm].unique():
            idx = idx + 1
            self.data.loc[self.data[chrm] == i, self.index] = int(idx)
            self.segments.loc[self.segments[chrm] == i, self.index] = int(idx)

        # Set the type to be the same as provided for chrm column
        self.data[self.index] = \
            self.data[self.index].astype(self.data[chrm].dtype)

        self.segments[self.index] = \
            self.segments[self.index].astype(self.segments[chrm].dtype)

        if self.nChr == 1:
            # For a single chromosome
            self.data[self.pos] = self.data[self.mp]
            self.ticks.append(int(len(self.data[self.pos]) / 2.) + 1)
            self.xlabel = "Chromosome %s position" % (self.data[chrm].unique())
            self.ticksLabels = self.ticks
        else:
            # For multiple chromosomes
            lastbase = 0
            seg_lastbase = 0
            for i in self.data[self.index].unique():
                if i == 1:
                    self.data.loc[self.data[self.index] == i, self.pos] = \
                        self.data.loc[self.data[self.index] == i, self.mp].values

                    self.segments.loc[self.segments[self.index] == i, self.start_pos] = \
                        self.segments.loc[self.segments[self.index] == i, self.start].values
                    self.segments.loc[self.segments[self.index] == i, self.end_pos] = \
                        self.segments.loc[self.segments[self.index] == i, self.end].values
                else:
                    prevbp = self.data.loc[self.data[self.index] == i - 1, self.mp]
                    seg_prevbp = self.segments.loc[self.segments[self.index] == i - 1, self.end]
                    # Shift the basepair position by the largest bp of the
                    # current chromosome
                    lastbase = lastbase + prevbp.iat[-1]
                    seg_lastbase = seg_lastbase + seg_prevbp.iat[-1]

                    self.data.loc[self.data[self.index] == i, self.pos] = \
                        self.data.loc[self.data[self.index] == i, self.mp].values \
                        + lastbase

                    self.segments.loc[self.segments[self.index] == i, self.start_pos] = \
                        self.segments.loc[self.segments[self.index] == i, self.start].values \
                        + lastbase

                    self.segments.loc[self.segments[self.index] == i, self.end_pos] = \
                        self.segments.loc[self.segments[self.index] == i, self.end].values \
                        + lastbase

                tmin = min(self.data.loc[self.data[self.index] == i, self.pos])
                tmax = max(self.data.loc[self.data[self.index] == i, self.pos])
                self.ticks.append(int((tmin + tmax) / 2.) + 1)

            self.xlabel = 'Chromosome'
            self.data[self.pos] = self.data[self.pos].astype(
                self.data[self.mp].dtype)

            self.ticksLabels = self.data[chrm].unique()  # All the ticks

    def figure(self,
               title="Copy Number Aberation",
               showgrid=True,
               xlabel=None,
               ylabel="Copy Ratio(log2)",
               point_size=5,
               showlegend=True,
               col=None):

        xmin = min(self.data[self.pos].values)
        xmax = max(self.data[self.pos].values)

        data_to_plot = []

        if self.nChr == 1:

            if col is None:
                col = ['black']

            # If single chromosome, ticks and labels automatic.
            layout = go.Layout(
                title=title,
                xaxis={
                    'title': self.xlabel if xlabel is None else xlabel,
                    'showgrid': showgrid,
                    'range': [xmin, xmax],
                },
                yaxis={'title': ylabel,
                       'range': [-2, 2]},
                hovermode='closest'
            )

            hover_text = _get_hover_text(self.data,
                                         genename=self.geneName,
                                         position=self.pos,
                                         log2=self.log2)

            data_to_plot.append(
                go.Scattergl(
                    x=self.data[self.pos].values,
                    y=self.data[self.log2].values,
                    mode="markers",
                    showlegend=showlegend,
                    marker={
                        'color': col[0],
                        'size': point_size,
                        'name': "chr%s" % self.data[self.chrName].unique()
                    },
                    text=hover_text
                )
            )
        else:
            # if multiple chrms, use the ticks and labels you created above.
            layout = go.Layout(
                title=title,
                xaxis={
                    'title': self.xlabel if xlabel is None else xlabel,
                    'showgrid': showgrid,
                    'range': [xmin, xmax],
                    'tickmode': "array",
                    'tickvals': self.ticks,
                    'ticktext': self.ticksLabels,
                    'ticks': "outside"
                },
                yaxis={'title': ylabel,
                       'range': [-2, 2]},
                hovermode='closest'
            )

            chrom_segments = dict(by_chromosome(self.segments))

            icol = 0
            if col is None:
                col = [
                    'black' if np.mod(i, 2)
                    else 'grey' for i in range(self.nChr)
                ]

            for i in self.data[self.index].unique():

                tmp = self.data[self.data[self.index] == i]

                chromo = tmp[self.chrName].unique()  # Get chromosome name

                hover_text = _get_hover_text(tmp,
                                             genename=self.geneName,
                                             position=self.pos,
                                             log2=self.log2,)

                data_to_plot.append(
                    go.Scattergl(
                        x=tmp[self.pos].values,
                        y=tmp[self.log2].values,
                        mode="markers",
                        showlegend=showlegend,
                        name="Chr%s" % chromo[0],
                        marker={
                            'color': col[icol],
                            'size': point_size
                        },
                        text=hover_text
                    )
                )

                for chrom, seg in chrom_segments[chromo[0]].iterrows():
                    data_to_plot.append(
                        go.Scattergl(
                            x=[seg[self.start_pos], seg[self.end_pos]],
                            y=[seg.log2, seg.log2],
                            mode="lines",
                            showlegend=False,
                            name="Chr%s" % chromo[0],
                            line={
                                'color': 'firebrick',
                                'width': 3
                            }
                        )
                    )

                icol = icol + 1

        return go.FigureWidget(data=data_to_plot, layout=layout)


cnr = pd.read_csv('data/tumor.merged.cnr', sep='\t')
seg = pd.read_csv('data/tumor.merged.cns', sep='\t')

app = dash.Dash(__name__)

app.layout = html.Div(children=[
    html.H1(children='Hello Dash'),

    html.Div(children='''
        Dash: A web application framework for Python.
    '''),

    dcc.Graph(
        id='example-graph',
        figure=CNVPlot(cnr=cnr, seg=seg)
    ),

    html.Div(
        dash_table.DataTable(
            id='table',
            columns=[{"name": i, "id": i} for i in cnr.columns[:-1]],
            data=cnr.to_dict("records"),
            fixed_rows={"headers": True, 'data': 0},
            style_cell={"width": '150px'}
        )
    ),
    html.Br(),
    html.Br()
])



if __name__ == '__main__':
    app.run_server(debug=True)
