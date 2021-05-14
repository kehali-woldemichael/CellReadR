# -*- coding: utf-8 -*-

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

app = Dash(__name__)

app.layout = html.Div([
                html.Label('Text 1: '),
                dcc.Input(id='input-1-submit', type='text', placeholder='Enter name'),
                html.Br(),
                html.Label('Text 2: '),
                dcc.Input(id='input-2-submit', type='text', placeholder='Enter surname'),
                html.Br(),html.Br(),
                html.Button('Submit', id='btn-submit'),
                html.Br(),
                html.Hr(),
                html.Label('Output'), html.Br(),html.Br(),
                html.Div(id='output-submit'),
                html.Br(), html.Hr()
            ])

@app.callback(Output('output-submit', 'children'),
                [Input('btn-submit', 'n_clicks')],
                [State('input-1-submit', 'value'),State('input-2-submit', 'value')])
def update_output(clicked, input1, input2):
    if clicked:
        x = input1 + " " + input2
        return 'Name is => ' + x

if __name__ == '__main__':
    app.run_server(debug=True)
