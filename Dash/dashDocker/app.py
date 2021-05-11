

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
# For working with sequence objects
from Bio.Seq import Seq
# import pandas as pd

import sys
# Importing module of personal functions
from kCellReadR import *

app = Dash(__name__)

app.layout = html.Div([
                html.H1('CellReadR: sesRNA selection'),
                html.H2('Sequence Selection'),
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

                html.Button('Submit', id='btn-geneInfo'),
                html.Br(), html.Br(),
                html.Label('Number of variants (CDS, cDNA): '),
                html.Div(id='output-geneInfo'),
                html.Div(id='output-cDNA'),
                html.Br(),

                html.Hr(),
                html.H2('Sequence search'),
                html.Label('Isoform: '),
                dcc.Input(id='input-3-submit',
                          type='number',
                          placeholder='Isoform #'),
                html.Br(),

                html.Label('Search seq: '),
                dcc.Input(id='input-4-submit',
                          type='text',
                          placeholder='CDS, cDNA, ...'),
                html.Br(),

                html.H2('sesRNA Features'),
                html.Label('Length sesRNA: '),
                dcc.Input(id='input-5-submit',
                          type='number',
                          placeholder='Between 200-300 bp'),
                html.Br(),

                html.Label('Min number of TGGs: '),
                dcc.Input(id='input-6-submit',
                          type='number',
                          placeholder='Min number of TGG'),
                html.Br(),

                html.Label('Max number of stop codons: '),
                dcc.Input(id='input-7-submit',
                          type='number',
                          placeholder='Max number of stops'),
                html.Br(),

                html.Label('Min GC content:'),
                dcc.Input(id='input-8-submit',
                          type='number',
                          placeholder='Min GC'),
                html.Br(),

                html.Label('Max GC content:'),
                dcc.Input(id='input-9-submit',
                          type='number',
                          placeholder='Max GC'),
                html.Br(),

                html.Label('TGG distance from center:'),
                dcc.Input(id='input-10-submit',
                          type='number',
                          placeholder='TGG bp from center'),
                html.Br(),

                html.Label('TGG distance from stop:'),
                dcc.Input(id='input-11-submit',
                          type='number',
                          placeholder='TGG bp from stop'),
                html.Br(),

                html.Label('ATG choice: '),
                dcc.Dropdown(id='input-12-submit',
                    options=[
                        {'label': 'None', 'value': 'None'},
                        {'label': 'All upstream', 'value': 'All upstream'},
                        {'label': 'Upstream central', 'value': 'Upstream central'}
                    ],
                    value=['None'],
                    placeholder='None',
                    style=dict(
                        width='40%',
                        display='inline-block',
                        verticalAlign='middle')
                ),
                html.Br(),

                html.Br(),html.Br(),
                html.Button('Submit', id='btn-submit'),
                html.Br(),html.Br(),

                html.Br(),
                html.Hr(),
                html.Label('Output'), html.Br(),html.Br(),
                html.Div(id='output-submit'),
                html.Br(), html.Hr(),

                html.Label('Sequence choice: '),
                dcc.Input(id='input-13-submit',
                          type='number',
                          placeholder='Chosen sequence'),
                html.Br(), html.Br(),

                html.Button('View sequence', id='seq-num'),
                html.Button('Download', id='btn-download'),
                html.Br(), html.Br(),
                html.Hr(),
                html.Label('Chosen sequence'), html.Br(),html.Br(),
                dcc.Textarea(
                    id='output-strSeq',
                    value='DNA sesRNA',
                    style={'width': '80%', 'height': 100},
                ),
                html.Hr(),

                html.Label('Number of TGGs to convert: '),
                dcc.Input(id='input-14-submit',
                          type='number',
                          placeholder='Number TGG'),
                html.Br(), html.Br(),
                html.Button('Convert sequence', id='btn-numTGG'),
                html.Br(), html.Br(),
                html.Label('Converted RNA'), html.Br(),html.Br(),
                dcc.Textarea(
                    id='outputRNA-strSeq',
                    value='RNA sesRNA',
                    style={'width': '80%', 'height': 100},
                ),
                dcc.Download(id='download-seq'),
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

@app.callback([Output('output-geneInfo', 'children'),
               Output('output-cDNA', 'children')],
                [Input('btn-geneInfo', 'n_clicks')],
                [State('input-1-submit', 'value'),
                 State('input-2-submit', 'value'),
                 State('input-3-submit', 'value'),
                 State('input-4-submit', 'value'),])
def generate_geneInfo(clicked, species, geneName, isoForm, searchSeq):
    if clicked:
        # Loading reference sequences
        rC_exon_records, C_exon_records, CDS, cDNA = \
            load_referenceSequences(geneName, species)

        return len(CDS), len(cDNA)
    else:
        return None, None

@app.callback(Output('output-submit', 'children'),
                [Input('btn-submit', 'n_clicks')],
                [State('input-1-submit', 'value'),
                 State('input-2-submit', 'value'),
                 State('input-3-submit', 'value'),
                 State('input-4-submit', 'value'),
                 State('input-5-submit', 'value'),
                 State('input-6-submit', 'value'),
                 State('input-7-submit', 'value'),
                 State('input-8-submit', 'value'),
                 State('input-9-submit', 'value'),
                 State('input-10-submit', 'value'),
                 State('input-11-submit', 'value'),
                 State('input-12-submit', 'value')])
def update_output(clicked, species, geneName, isoForm, searchSeq, lenSes, numTGG, numStop,
                  minGC, maxGC, nearCenter, fromStop, atgChoice):
    if clicked:
        # Loading reference sequences
        rC_exon_records, C_exon_records, CDS, cDNA = \
            load_referenceSequences(geneName, species)

        rC_parameters = parameters_sesRNA('Reverse', (isoForm-1), lenSes, numTGG,
                                          numStop, atgChoice,
                                          minGC, maxGC,
                                          nearCenter, fromStop)

        if searchSeq == 'CDS': targetSeq = CDS
        elif searchSeq == 'cDNA': targetSeq = cDNA
        rC_multiExon_sesRNAs, rC_sequenceMetrics, rC_sesRNA_objs = \
            generate_sesRNAs_multiExon(rC_exon_records, targetSeq, rC_parameters)

        return generate_table(rC_sequenceMetrics, max_rows=10)

@app.callback(Output("output-strSeq", "value"),
                [Input("seq-num", "n_clicks")],
                [State('input-1-submit', 'value'),
                 State('input-2-submit', 'value'),
                 State('input-3-submit', 'value'),
                 State('input-4-submit', 'value'),
                 State('input-5-submit', 'value'),
                 State('input-6-submit', 'value'),
                 State('input-7-submit', 'value'),
                 State('input-8-submit', 'value'),
                 State('input-9-submit', 'value'),
                 State('input-10-submit', 'value'),
                 State('input-11-submit', 'value'),
                 State('input-12-submit', 'value'),
                 State('input-13-submit', 'value')],
                 prevent_initial_call=True,
)
def display_chosenSeq(clicked, species, geneName, isoForm, searchSeq, lenSes, numTGG, numStop,
                  minGC, maxGC, nearCenter, fromStop, atgChoice, seqNum):
    if clicked:
        filename = "sesRNA.txt"
        # Alternatively:
        # filename = f"{uuid.uuid1()}.txt"

        # Loading reference sequences
        rC_exon_records, C_exon_records, CDS, cDNA = \
            load_referenceSequences(geneName, species)

        parameters = parameters_sesRNA('Reverse', (isoForm-1), lenSes, numTGG,
                                          numStop, atgChoice,
                                          minGC, maxGC,
                                          nearCenter, fromStop)

        if searchSeq == 'CDS': targetSeq = CDS
        elif searchSeq == 'cDNA': targetSeq = cDNA
        multiExon_sesRNAs, sequenceMetrics, sesRNA_objs = \
            generate_sesRNAs_multiExon(rC_exon_records, targetSeq, parameters)

        output_strSeq = str(multiExon_sesRNAs[seqNum-1])

        return output_strSeq

@app.callback(Output("outputRNA-strSeq", "value"),
                [Input("btn-numTGG", "n_clicks")],
                [State('input-1-submit', 'value'),
                 State('input-2-submit', 'value'),
                 State('input-3-submit', 'value'),
                 State('input-4-submit', 'value'),
                 State('input-5-submit', 'value'),
                 State('input-6-submit', 'value'),
                 State('input-7-submit', 'value'),
                 State('input-8-submit', 'value'),
                 State('input-9-submit', 'value'),
                 State('input-10-submit', 'value'),
                 State('input-11-submit', 'value'),
                 State('input-12-submit', 'value'),
                 State('input-13-submit', 'value'),
                 State('input-14-submit', 'value')],
                 prevent_initial_call=True,
)
def display_chosenRNA(clicked, species, geneName, isoForm, searchSeq, lenSes, numTGG, numStop,
                  minGC, maxGC, nearCenter, fromStop, atgChoice, seqNum, numTGG_convert):
    if clicked:
        filename = "sesRNA.txt"
        # Alternatively:
        # filename = f"{uuid.uuid1()}.txt"

        # Loading reference sequences
        rC_exon_records, C_exon_records, CDS, cDNA = \
            load_referenceSequences(geneName, species)

        parameters = parameters_sesRNA('Reverse', (isoForm-1), lenSes, numTGG,
                                          numStop, atgChoice,
                                          minGC, maxGC,
                                          nearCenter, fromStop)

        if searchSeq == 'CDS': targetSeq = CDS
        elif searchSeq == 'cDNA': targetSeq = cDNA
        multiExon_sesRNAs, sequenceMetrics, sesRNA_objs = \
            generate_sesRNAs_multiExon(rC_exon_records, targetSeq, parameters)

        outputRNA_strSeq = str(convert_DNA(multiExon_sesRNAs[seqNum-1], numTGG_convert))

        return outputRNA_strSeq

@app.callback(Output("download-seq", "data"),
                [Input("btn-download", "n_clicks")],
                [State('input-1-submit', 'value'),
                 State('input-2-submit', 'value'),
                 State('input-3-submit', 'value'),
                 State('input-4-submit', 'value'),
                 State('input-5-submit', 'value'),
                 State('input-6-submit', 'value'),
                 State('input-7-submit', 'value'),
                 State('input-8-submit', 'value'),
                 State('input-9-submit', 'value'),
                 State('input-10-submit', 'value'),
                 State('input-11-submit', 'value'),
                 State('input-12-submit', 'value'),
                 State('input-13-submit', 'value')],
                 prevent_initial_call=True,
)
def create_download_file(clicked, species, geneName, isoForm, searchSeq, lenSes, numTGG, numStop,
                  minGC, maxGC, nearCenter, fromStop, atgChoice, seqNum):
    if clicked:
        filename = "sesRNA.txt"
        # Alternatively:
        # filename = f"{uuid.uuid1()}.txt"

        # Loading reference sequences
        rC_exon_records, C_exon_records, CDS, cDNA = \
            load_referenceSequences(geneName, species)

        parameters = parameters_sesRNA('Reverse', (isoForm-1), lenSes, numTGG,
                                          numStop, atgChoice,
                                          minGC, maxGC,
                                          nearCenter, fromStop)

        if searchSeq == 'CDS': targetSeq = CDS
        elif searchSeq == 'cDNA': targetSeq = cDNA
        multiExon_sesRNAs, sequenceMetrics, sesRNA_objs = \
            generate_sesRNAs_multiExon(rC_exon_records, targetSeq, parameters)

        output_strSeq = str(multiExon_sesRNAs[seqNum-1])

        return dict(content=output_strSeq, filename=filename)

@app.server.route('/dash/urlToDownload')
def download_csv():
    return send_file('output/downloadFile.csv',
                     mimetype='text/csv',
                     attachment_filename='downloadFile.csv',
                     as_attachment=True)

if __name__ == "__main__":
    app.run_server(host="0.0.0.0", port=8080, debug=True)

