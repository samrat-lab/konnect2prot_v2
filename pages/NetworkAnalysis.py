import pandas as pd
from dash import dcc, html, dash_table
import dash_cytoscape as cyto
from utils.panel import *


#cyto.load_extra_layouts()
singal_tooltip='Black nodes are signaling pathways\nRed nodes are proteins\n\nBCR\t Signaling by B cell receptor\n HH\t Signaling by Hedgehog\n HIPPO\t Signaling by Hippo\n IIP\t Innate immune pathways\n JAK/STAT\t Janus Activating Kinase/ Signal Transducer and Activator of Transcription\n Notch\t Signaling by Notch\n TCR\t Signaling by T cell receptor\n TGF\t Signaling by Transforming Growth Factor beta\n TLR\t  Toll-Like receptor signaling\n TNF\t Signaling by Tumor necrosis factor\n NHR\t Signaling by Nuclear hormone receptor\n RTK\t Signaling by Receptor tyrosine kinase\n WNT/Wingless\t Signaling by WNT\n'
example_prots='RB1\nPPP1CA\nTP53\nCSNK2A1\nCDK1\nCHEK1\nEEF2K\nEGFR\nERBB2\nCDC7\nAR\nBRCA1\nMAPK6\nSIRT1\nNME1\nEIF2AK2'
colors = [['red','tee'], ['blue','triangle'], ['grey','none']]
TABLE_COLUMNS=['Name','Degree','InDegree','OutDegree','Betweenness','Clustering','Closeness']
result_columns=['Name', 'Degree', 'InDegree', 'OutDegree', 'Betweenness', 'Clustering',
       'Closeness' ,'Protein Class', 'Location', 'Inhibitors',
       'PDB complex(Count)']

tabs_styles = {
    'height': '44px'
}
tab_style = {
    'borderBottom': '1px solid #d6d6d6',
    'padding': '6px',
    'fontWeight': 'bold'
}

tab_selected_style = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': '#119DFF',
    'color': 'white',
    'padding': '6px'
}

sub_head_css = {
    'textAlign': 'center',
    'font':'Courier New',
    'font-weight':'900',
    'color':'white',
    'font-size':'x-large'
}


#client = MongoClient("mongodb://127.0.0.1:27017")
#clinet = MongoClient("localhost:27017")
#db = client.ppi
#################function starts#############################

def generate_network(cytoscape_id='cytoscape-update-layout', message='Generating Network...'):
    """
    Generate a Cytoscape network visualization with dummy data and custom styling.

    Parameters:
    cytoscape_id (str): The ID to assign to the Cytoscape component.
    message (str): A message to print when the function is called.

    Returns:
    cyto.Cytoscape: A Dash Cytoscape component with predefined nodes, edges, and styling.
    """
    print(f'====================================== {message} ==============================')

    # Generate dummy nodes
    node_data = [
        ('la', 'PPI', 34.03, -118.25),
        ('nyc', 'Pathways', 40.71, -74),
        ('to', 'Processes', 43.65, -79.38),
        ('mtl', 'Function', 45.50, -73.57),
        ('van', 'Disease', 49.28, -123.12),
        ('chi', 'Ligand', 41.88, -87.63),
        ('bos', 'Localization', 42.36, -71.06),
        ('hou', 'Mutation', 29.76, -95.37)
    ]
    nodes = [
        {
            'data': {'id': short, 'label': label},
            'position': {'x': 20 * lat, 'y': -20 * long}
        }
        for short, label, long, lat in node_data
    ]

    # Generate dummy edges
    edge_data = [
        ('van', 'la'), ('la', 'chi'), ('hou', 'chi'), 
        ('to', 'mtl'), ('mtl', 'bos'), ('nyc', 'bos'), 
        ('to', 'hou'), ('to', 'nyc'), ('la', 'nyc'), 
        ('nyc', 'bos')
    ]
    edges = [
        {'data': {'source': source, 'target': target}}
        for source, target in edge_data
    ]

    # Combine nodes and edges into elements
    elements = nodes + edges

    # Define styles
    default_node_style = {
        "selector": 'node',
        'style': {
            "opacity": 1,
            'height': 15,
            'width': 15,
            'background-color': '#222222',
            "label": "data(label)"
        }
    }
    default_edge_style = {
        "selector": 'edge',
        'style': {
            'target-arrow-color': 'black',
            "curve-style": "bezier",
            "opacity": 0.3,
            'width': 2
        }
    }
    custom_styles = [
        {
            'selector': '.search',
            'style': {'background-color': 'red'}
        }
    ]
    color_styles = [
        {
            "selector": f'.{color}',
            'style': {'line-color': color, "target-arrow-shape": arrow}
        }
        for color, arrow in colors
    ]

    # Combine all styles into a stylesheet
    default_stylesheet = [default_node_style, default_edge_style] + custom_styles + color_styles

    # Return Cytoscape component
    return cyto.Cytoscape(
        id=cytoscape_id,
        zoom=0.5,
        stylesheet=default_stylesheet,
        layout={'name': 'circle'},
        style={'width': '100%', 'height': '700px'},
        elements=elements
    )



my_layout = html.Div([
# Panel for Network Visualization
html.Div([
    html.Div([
            dcc.Loading(id='cyto_loading', children=generate_network('cytoscape-update-layout', 'Main network')),
            html.Div(id='net_info', className='four columns'),
            html.Div(id='intermediate-value', style={'display': 'none'}),
            html.Div(id='unmatch-value', style={'display': 'none'}),
            dcc.Store(id='memory_net', storage_type='session'),
        ],className="grid-item full-width")
    ],className="dashboard-grid"),

# Panel for Local Properties
html.Div([
        html.Div([
            html.P("Local Properties (Click Node/Edge to retrieve)",className="panel-title"),
            html.Div(id='func_info', style={'display': 'none'})
            ],className="grid-item full-width")  # Remove?
    ],className="dashboard-grid"),

    # Panel for Tabs (Node info, PDB Complex, etc.)
html.Div([
        html.Div([
            dcc.Tabs(id="tabs-styled-with-inline", value='tab-0', children=[
                dcc.Tab(label='Node info', value='tab-0', style=tab_style, selected_style=tab_selected_style),
                dcc.Tab(label='PDB Complex', value='tab-1', style=tab_style, selected_style=tab_selected_style),
                dcc.Tab(label='Ligands', value='tab-2', style=tab_style, selected_style=tab_selected_style),
                dcc.Tab(label='Mutation (disease)', value='tab-3', style=tab_style, selected_style=tab_selected_style),
                dcc.Tab(label='Residue information', value='tab-4', style=tab_style, selected_style=tab_selected_style),
                dcc.Tab(label='Interaction trivia', value='tab-5', style=tab_style, selected_style=tab_selected_style),
            ], style=tabs_styles),
            dcc.Loading(html.Div(id='tabs-content-inline'))
        ],className="grid-item full-width")
    ],className="dashboard-grid"),


    # Section 1: Global Properties
html.Div([
        html.Div([
            html.P("Global Properties", className="panel-title")
            ],className="grid-item full-width"),
    ], className="dashboard-grid"),

    # Section 2: Enrichment Panel
html.Div("Enrichment panel", className="banner"),
html.Div([
        html.Div([
            html.Div("Pathway Graph", className="panel-title"),
            dcc.Loading(dcc.Graph(id="pathway_graph")),
        ], className="grid-item"),
        html.Div([
            html.Div("Process Graph", className="panel-title"),
            dcc.Loading(dcc.Graph(id="process_graph")),
        ], className="grid-item"),
    ], className="dashboard-grid"),

    # Section 3: Data panel
html.Div([
        html.Div([
            html.Div("Protein classes", className="panel-title"),
            dcc.Loading(dcc.Graph(id="pclass_graph")),
        ], className="grid-item"),
        html.Div([
            html.Div("Disease enrichment", className="panel-title"),
            dcc.Loading(dcc.Graph(id="disease_graph")),
        ], className="grid-item"),
    ], className="dashboard-grid"),

    # Section 4: Topological Panel
html.Div("Topological panel", className="banner"),
html.Div([
        html.Div([
            html.Div("Complex Parameters", className="panel-title"),
            dash_table.DataTable(
                id="central_para",
                columns=[{"name": i, "id": i} for i in TABLE_COLUMNS],
                sort_action="native",
                filter_action="native",
                page_action="native",
                page_current=0,
                page_size=12,
                style_cell={"fontSize": 12, "font-family": "Arial"},
                style_table={"overflowX": "auto"},
                style_header={
                    "backgroundColor": "rgb(230, 230, 230)",
                    "fontWeight": "bold",
                },
            ),
        ], className="grid-item"),
        html.Div([
            html.Div("Degree vs Betweenness", className="panel-title"),
            dcc.Loading(dcc.Graph(id="network_para_graph")),
        ], className="grid-item"),
    ], className="dashboard-grid"),

# Section 5.1: Results for Top Spreaders
html.Div("Results for top spreaders", className="banner"),
html.Div([
        html.Div([
            html.Div("Top Stats", className="panel-title"),
            dash_table.DataTable(
                id="top_stats",
                columns=[{"name": i, "id": i} for i in result_columns],
                sort_action="native",
                page_action="native",
                page_current=0,
                page_size=12,
                style_data={"whiteSpace": "normal", "height": "auto", "lineHeight": "15px"},
                style_cell={"fontSize": 12, "font-family": "Arial"},
                style_table={"overflowX": "auto", "table-layout": "fixed"},
                style_header={"backgroundColor": "rgb(230, 230, 230)", "fontWeight": "bold"},
                editable=True,
                export_format="xlsx",
            ),
        ], className="grid-item full-width"),
    ], className="dashboard-grid"),

# Section 5.2: Pathway cluster
html.Div([
        html.Div([
            html.Div("Pathway clustergram", className="panel-title"),
            dcc.Graph(id="top_pathway"),
        ], className="grid-item full-width")
    ], className="dashboard-grid"),

# Section 6: Additional Graphs
html.Div([
        html.Div([
            html.Div("Tissue specificity", className="panel-title"),
            dcc.Graph(id="location_graph"),
        ], className="grid-item"),
        html.Div([
            html.Div("Location of top spreaders", className="panel-title"),
            dcc.Graph(id="sub_network_para_graph"),
        ], className="grid-item"),
    ], className="dashboard-grid"),

# Section 7: Hallmarks and Signaling
html.Div([
        # Hallmark Network
        html.Div([
            html.Div("Hallmark Network", className="panel-title"),
            generate_network('hallmark_network', 'hallmark network'),
        ], className="grid-item"),
        # Signalling Network
        html.Div([
            html.Div("Signalling Network", className="panel-title"),
            generate_network('signa_network', 'signalling network'),
        ], className="grid-item"),
    ], className="dashboard-grid"),

    # Section 8: Footer
    # html.Div([
    #     html.Div(html.P(["Copyright ©2024 Dr. Samrat Chatterjee's LAB (THSTI)"])),
    # ], id="mainContainer", style={"display": "flex", "flex-direction": "column"})
])





