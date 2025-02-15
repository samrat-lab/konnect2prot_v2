from utils.preprocess_pipeline import preprocess_and_find_DEGs
from dash import dcc, html, dash_table, callback
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash import html, dcc, dash_table
from flask_session import Session
from flask import session
from flask import request 
from pymongo import MongoClient
from skimage.io import imread
#from app import app
from functools import reduce
from uuid import uuid4 
import plotly.graph_objects as go
import dash_bio as dashbio
import dash_uploader as du
import plotly.express as px
import pandas as pd
import networkx as nx
import numpy as np
import collections
import io
import json
import dash
import requests
import base64



def get_tab1_sidebar():
    my_div_tab1=html.Div([
        html.H4("Instructions",className='panel-header') ,
        html.Div([html.B("1.Data Tab"),html.P("To start the analysis and for uploading your data")],className='input-label'),
        html.Div([html.B("2. EDA Tab:"),html.P("To find DEGs and Quantitative Exploration of the data ")],className='input-label'),
        html.Div([html.B("3. Visualization  Tab:"),html.P("Qualitative Exploration such as Volcano plot , PCA and Pathway Analysis")],className='input-label'),
        html.Div([html.B("4. Network Tab:"),html.P("Previous version of konnect2prot")],className='input-label'),
        html.Div([html.B("*For more information go to tutorial page")],className='input-label')],className='panel-container')
    return my_div_tab1


def get_data_sidebar(sess_id):
    my_div_data = html.Div([
        html.Div([
                #html.H1(f'Menu',id = 'menu_status'),
                
                html.H2('Upload Your File (.csv , .tsv , .xls)'),
                html.P('Preferable file size < 100 Mb'),
                 #-----------------------------Review Code ---------------------------------------------
                html.Div([du.Upload(
                    text='Click To Import',
                    text_completed='Completed: ',
                    pause_button=False,
                    cancel_button=True,
                    max_file_size=50000,  # 1800 Mb
                    filetypes=['csv', 'tsv' , 'xls'],
                    id='upload-data',
                    upload_id = sess_id
                    #upload_id = session.get('static_uuid')
                )]),
                html.Div(id='callback-output', style={'display': 'none'}),
                html.Button("Download Sample File", id="download-btn", n_clicks=0, style={"padding": "10px", "font-size": "16px"},className="highlight-anchor"),
                html.Div(['Or prepare your data with row as feature/gene/protein and sample in column']),
                dcc.Download(id="file-download")
                ##html.Div(id='upload-data-status')  # Placeholder for upload status messages
            ], className="p-4 max-w-lg mx-auto bg-white shadow-lg rounded-lg"),
        html.Div(
            [
                html.Div([
                    html.H4('Outlier Removal', className='panel-header'),
                    html.Label("Lower Percentile for Flooring:  "),
                    dcc.Input(id="floor-percentile", type="number", min=0, max=100, step=0.1, value=0,className='input-box'),
                    html.Br(),
                    html.Label("Upper Percentile for Capping:  "),
                    dcc.Input(id="cap-percentile", type="number", min=0, max=100, step=0.1, value=100,className='input-box'),
                    html.Br(),
                    html.Button("Apply Outlier Removal", id="apply-outlier", n_clicks=0, style={"padding": "10px", "font-size": "16px"},className="submit-button"),
                ], className="p-4 max-w-lg mx-auto my-4 bg-white shadow-lg rounded-lg"),
            ],
            className='row'
        ),])
    return my_div_data
         #-----------------------------Review Code ---------------------------------------------
def get_eda_sidebar(uploaded_data_store):

    if uploaded_data_store is None:
        return html.Div("No data uploaded.")

    df = pd.DataFrame.from_records(uploaded_data_store)
    columns = df.columns[1:]  # Exclude the first column


    my_div_eda = html.Div([


    # Row for Preprocessing Parameters
    html.Div([
        html.Div([
            html.H4('Preprocessing Parameters', className='panel-header'),

            html.H5('Select Normalization Method', className='input-label'),
            dcc.Dropdown(
                id='norm-method-dropdown',
                options=[
                    {'label': 'log2', 'value': 'log2'},
                    {'label': 'log10', 'value': 'log10'},
                    {'label': 'none', 'value': 'none'}
                ],
                value='none',
                placeholder='Select Normalization Method',
                className='dropdown'
            ),

            html.H5('Select Scaling Method', className='input-label'),
            dcc.Dropdown(
                id='scaling-method-dropdown',
                options=[
                    {'label': 'Standard Scaling', 'value': 'standard'},
                    {'label': 'Min-Max Scaling', 'value': 'minmax'},
                    {'label': 'none', 'value': 'none'}
                ],
                value='none',
                placeholder='Select Scaling Method',
                className='dropdown'
            ),

            html.H5('P Value (Significance Threshold)', className='input-label'),
            dcc.Input(id='alpha-input', type='number', value=0.05, step=0.01, className='input-box'),
                ], className='panel-container')
            ], className='row'),
    # Row for Group Naming
    html.Div([
        html.Div([
            html.H4('Define Groups', className='panel-header'),

            html.H5('Provide a name for Set 1 (e.g., Disease)', className='input-label'),
            dcc.Input(id='set1-name', type='text', placeholder='Enter name for Set 1', className='input-box'),

            html.H5('Provide a name for Set 2 (e.g., nonDisease)', className='input-label'),
            dcc.Input(id='set2-name', type='text', placeholder='Enter name for Set 2', className='input-box'),

            html.H5('Select columns for Set 1', className='input-label'),
            dcc.Dropdown(
                id='set1-dropdown',
                options=[{'label': col, 'value': col} for col in columns],
                multi=True,
                placeholder='Select columns for Set 1',
                className='dropdown'
            ),

            html.H5('Columns for Set 2 (auto-populated)', className='input-label'),
            dcc.Dropdown(
                        id='set2-dropdown',
                        options=[],
                        multi=True,
                        disabled=True,
                        className='dropdown'
                    ),
                ], className='panel-container')
            ], className='row'),

    

    # Row for Submit Button
    html.Div([
        html.Button('Submit', id='submit-columns', n_clicks=0, className='submit-button')
    ], className='row'),

    # # Row for Output
    # html.Div([
    #     html.Div(id='column-output')
    # ], className='row'),

    dcc.Download(id='download-degs')  # Download component
])


    return my_div_eda

def get_vis_sidebar():
    my_div_visualize=html.Div([html.P("log2(Fold Change) Threshold",className='panel-header'),
                                        dcc.Slider(
                                        id='fold-change-threshold',
                                        min=0,
                                        max=5,
                                        step=0.1,
                                        value=0.0,
                                        marks={i: f'{i}' for i in range(6)},
                                        tooltip={"placement": "bottom", "always_visible": True}
                                    ),
                                html.P("Enrichment Type",className='panel-header'),
                                dcc.Dropdown(id='enrich_cat',clearable=False,
                                options=[{'label': 'Process','value': 'GO_Biological_Process_2018'},{'label':'Pathway','value':'KEGG_2019_Human'},{'label':'Disease','value':'ClinVar_2019'}],
                                value = 'KEGG_2019_Human'
                             )],className='panel-container')
    return my_div_visualize


def get_net_sidebar():
    my_div_network = html.Div([
        html.P([
            html.P('Enter Protein Name(s) or'),
            html.A(id='example_button', children='Try Example', className="ml-2 text-blue-500 underline")
        ], className="items-center space-x-2"),
        
        dcc.RadioItems(
            id='query_type', 
            options=[
                {'label': 'Gene Symbol', 'value': 'gene_symbol'}, 
                {'label': 'Uniprot ID', 'value': 'uniprot_id'}
            ],
            value='gene_symbol', 
            labelStyle={'display': 'inline-block'},
            className="flex justify-center space-x-4 mt-2"
        ),
        #dcc.Input(className="custom-input"),
        #html.Button(className="custom-button"),
        dcc.Textarea(
            id='Textarea Input',
            placeholder='Type here or click the Try Example',
            value='',
            #className="border rounded-md w-full p-2 mt-2 resize-y",
            className = "custom-input",
            style={"height": "100px", "min-height": "50px"}  # Allows for resizing
        ),
        
        html.Div('OR', className="text-center font-bold mt-4"),
        
        dcc.Upload(
            id='upload-prot_list',
            children=html.Div(['Drag and Drop a list or ', html.A('Select File', className="text-blue-500 underline")], 
                              className="form-control"),
            className="w-full h-14 border border-dashed border-gray-400 rounded-lg text-center flex items-center justify-center my-4"
        ),
        
        dcc.Dropdown(
            id='dropdown-update-layout',
            value='grid',
            clearable=False,
            className="dcc_control border border-gray-300 rounded-md w-full p-2",
            options=[{'label': name.capitalize(), 'value': name} for name in ['klay', 'dagre', 'spread', 'cola', 'euler', 'breadthfirst', 'cose-bilkent', 'preset', 'grid', 'random', 'circle', 'cose', 'concentric']]
        ),
        
        html.Button(children='Search', id='save button', type='submit', n_clicks=0, className='bg-blue-600 text-white py-2 px-4 rounded-lg mt-4 w-full'),
        #html.Button(children='Search', id='save button', type='submit', n_clicks=0, className='custom-button'),
        
        dcc.Download(id="unmatch_prot"),
        
        html.Div('Choose filter(s)', className='control_label font-semibold mt-4'),
        
        html.Label(["Location", dcc.Dropdown(id="location", multi=False, className="border border-gray-300 rounded-md w-full mt-1")]),
        html.Label(["Molecular Function", dcc.Dropdown(id="function", multi=False, className="border border-gray-300 rounded-md w-full mt-1")]),
        html.Label(["Biological Process", dcc.Dropdown(id="process", multi=False, className="border border-gray-300 rounded-md w-full mt-1")]),
        html.Label(["Pathways", dcc.Dropdown(id="pathway", multi=False, className="border border-gray-300 rounded-md w-full mt-1")]),
        html.Label(["Disease", dcc.Dropdown(id="disease", multi=False, className="border border-gray-300 rounded-md w-full mt-1")]),
        html.Label(["Tissue expression", dcc.Dropdown(id="tissue", multi=False, className="border border-gray-300 rounded-md w-full mt-1")]),
        
        html.Div([
            html.Button(children='Reset Filters', id='reset', type='submit', n_clicks=0, className='bg-gray-300 text-gray-700 py-1 px-4 rounded-lg text-sm'),
            html.Button(children='Apply Filters', id='apply_fil', type='submit', n_clicks=0, className='bg-blue-500 text-white py-1 px-4 rounded-lg text-sm')
        ], className='flex px-8 space-x-4 mt-4'),
        
        html.Br(),
        
        html.Button(children='Analyse', id='cp_calculate', type='submit', n_clicks=0, className='bg-green-500 text-white py-2 px-4 rounded-lg mx-auto my-4 flex justify-center items-center'),
        
        dcc.Loading(html.Div(id='analyse_status'))
    ], className="p-4 max-w-lg mx-auto bg-white shadow-lg rounded-lg")

    return my_div_network


