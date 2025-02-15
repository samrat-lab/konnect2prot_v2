from dash import dcc, html, get_asset_url
import dash_cytoscape as cyto
cyto.load_extra_layouts()


def index_layout(app):
    # Define the navbar
    navbar = html.Div([
        #html.Img(src=app.get_asset_url("k2p_logo4.png"), className="navbar-logo"),  # Add logo image
        html.A(href="/",children=[html.Img(src=app.get_asset_url("New_logo_k2p_caps.png"), className="navbar-logo")])
        
    ], className="navbar")
    tab_selected_style = {
        'backgroundColor': '#dec209',
        'color': 'black',
        'padding': '6px',
        'text-align': 'center',
        'text-overflow': 'ellipsis'
    }
    # Define the tab bar
    tab_bar = dcc.Tabs(
        id='tabs',
        children=[
            dcc.Tab(label='Data', value='data-tab', className='tab-button',  selected_style=tab_selected_style),
            dcc.Tab(label='EDA', value='eda-tab', className='tab-button',  selected_style=tab_selected_style),
            dcc.Tab(label='Visualization', value='visualization-tab', className='tab-button',  selected_style=tab_selected_style),
            dcc.Tab(label='Network Analysis (k2p 1.0)', value='network-tab', className='tab-button',  selected_style=tab_selected_style),
        ],
        className='tab-bar',
        #persistence=True,
        #persistence_type="session"
    )

    # Define the sidebar container
    sidebar = html.Div(id='my_sidebar', className='sidebar')
    network_tab_content_hidden = html.Div([
        html.Div(id='func_info'),
        html.Div(id='click_info'),
        html.Div(id='docking'),
        html.Div(id='ligand'),
        html.Div(id='pocket_info'),
        html.Div(id='disease_mutation'),
        ],className='row',style={'display': 'none'})

# Define the full layout

    layout = html.Div([
        navbar,
        tab_bar,
        sidebar,
        network_tab_content_hidden,
        html.Div(id = 'page-content', className='content'),
        html.Div(id = 'callback-output', style={'display': 'none'}),  # Hidden div to store file path
        dcc.Store(id='n_clicks_store', data=None),
        dcc.Store(id = 'degs-store'),# Storing the deg
        dcc.Store(id = 'data_label'),# Selected disease label by user
        dcc.Store(id = 'callback-counter',data = 0),# Not required in program, debug feature
        dcc.Store(id='central_para'),# Consider removing [with caution], not required
        dcc.Store(id = 'upload_status')# Store the file upload status

    ])
    return layout
