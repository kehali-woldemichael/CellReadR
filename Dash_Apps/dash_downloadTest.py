# -*- coding: utf-8 -*-

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
# For working with sequence objects
from Bio.Seq import Seq
import pandas as pd

import sys
# Importing module of personal functions
sys.path.append(
    '/home/user1/Dropbox/Research/Neurobiology_PhD/Rotations/Huang/Projects/CellReadR/Code/')
from kCellReadR import *

app = Dash(__name__)

app.layout = html.Div([
                html.Label('Species: '),
                dcc.Input(id='input-1-submit',
                          type='text',
                          placeholder='Enter species'),
                html.Br(),

                html.Label('Gene Name: '),
                dcc.Input(id='input-2-submit',
                          type='text',
                          placeholder='Enter gene name'),
                html.Br(), html.Br(),

                html.Label('Length sesRNA: '),
                dcc.Input(id='input-3-submit',
                          type='number',
                          placeholder='Between 200-300 bp'),
                html.Br(),

                html.Label('Min number of TGGs: '),
                dcc.Input(id='input-4-submit',
                          type='number',
                          placeholder='Min number of TGG'),
                html.Br(),

                html.Label('Max number of stop codons: '),
                dcc.Input(id='input-5-submit',
                          type='number',
                          placeholder='Max number of stops'),

                html.Br(),html.Br(),
                html.Button('Submit', id='btn-submit'),
                html.Button('Download', id='btn-download'),
                html.Br(),
                html.Hr(),
                html.Label('Output'), html.Br(),html.Br(),
                html.Div(id='output-submit'),
                html.Br(), html.Hr(),

                html.A(id='download-link', children='Download File'),
                dcc.Dropdown(
                    id='dropdown',
                    options=[{'label': i, 'value': i} for i in ['NYC', 'LA' 'SF']],
                    value='NYC',
                    clearable=False
                )
            ])

def generate_table(df, max_rows=10):
    return html.Table(className="responsive-table",
                      children=[
                          html.Thead(
                              html.Tr(
                                  children=[
                                      html.Th(col.title()) for col in df.columns.values],
                                  style={'color':app_colors['text']}
                                  )
                              ),
                          html.Tbody(
                              [
                              html.Tr(
                                  children=[
                                      html.Td(data) for data in d
                                      ], style={'color':app_colors['text'],
                                                'background-color':quick_color(d[2])}
                                  )
                               for d in df.values.tolist()])
                          ]
    )


app_colors = {
    'background': '#FFFFFF',
    'text': '#000000',
    'sentiment-plot':'#41EAD4',
    'volume-bar':'#FBFC74',
    'someothercolor':'#FF206E',
}


def quick_color(s):
    # except return bg as app_colors['background']
    if s >= POS_NEG_NEUT:
        # positive
        return "#FFFFFF"
    elif s <= -POS_NEG_NEUT:
        # negative:
        return "#FFFFFF"

    else:
        return app_colors['background']


def pos_neg_neutral(col):
    if col >= POS_NEG_NEUT:
        # positive
        return 1
    elif col <= -POS_NEG_NEUT:
        # negative:
        return -1

    else:
        return 0


POS_NEG_NEUT = 0.1

MAX_DF_LENGTH = 100


@app.callback(Output('output-submit', 'children'),
              Output('download-link', 'href'),
                [Input('btn-submit', 'n_clicks')],
                [State('input-1-submit', 'value'),
                 State('input-2-submit', 'value'),
                 State('input-3-submit', 'value'),
                 State('input-4-submit', 'value'),
                 State('input-5-submit', 'value')])
def update_output(clicked, species, geneName, lenSes, numTGG, numStop):
    if clicked:
        # Loading reference sequences
        rC_exon_records, C_exon_records, CDS, cDNA = \
            load_referenceSequences(geneName, species)

        rC_parameters = parameters_sesRNA('Reverse', 0, lenSes, numTGG, numStop, 'All upstream')
        rC_multiExon_sesRNAs, rC_sequenceMetrics, rC_sesRNA_objs = \
            generate_sesRNAs_multiExon(rC_exon_records, CDS, rC_parameters)

        return generate_table(rC_sequenceMetrics, max_rows=10)

@app.callback(Output('download-link', 'href'),
              [Input('btn-download', 'n_clicks')])
def download_data(clicked):
    if clicked:
        csv_string = rC_sequenceMetrics.to_csv(index = False, endcoding = 'utf-8')
        return csv_string

if __name__ == '__main__':
    app.run_server(debug=True)
