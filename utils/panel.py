from dash import dcc, html, dash_table
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from io import BytesIO
import base64
import scipy.stats as stats

def create_panel1(title, image_src, loading_id):
    """Helper function to create a panel layout for plots."""
    return html.Div([
        html.H4(title),
        dcc.Loading(
            id=loading_id,
            type="circle",
            children=html.Div([
                html.Img(src=image_src, style={'width': '100%', 'margin': 'auto'})
            ], style={'textAlign': 'center'})
        )
    ], style={
        'width': '48%', 'display': 'inline-block', 'border': '1px solid #ddd', 
        'padding': '20px', 'margin': '20px 10px', 'border-radius': '5px', 
        'box-shadow': '2px 2px 10px rgba(0,0,0,0.1)'
    })

   
# Reusable panel component
panel_style = {
    'padding': '10px',
    'margin': '10px',
    'boxShadow': '0 3px 5px rgba(57, 63, 72, 0.3)',
    'borderRadius': '8px',
    'backgroundColor': '#f8f9fa'
}
def create_panel(title, content):
    return html.Div([
        html.H5(title, style={'text-align': 'center', 'margin-bottom': '10px'}),
        content
    ], style=panel_style, className='pretty_container row')

