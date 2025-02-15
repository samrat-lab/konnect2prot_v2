from utils.preprocess_pipeline import preprocess_and_find_DEGs
from statsmodels.stats.multitest import multipletests
from dash import dcc, html, dash_table, callback
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from plotly.graph_objects import Figure
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from scipy.stats import pearsonr
from pymongo import MongoClient
import plotly.express as px
import pandas as pd
import numpy as np
import dash
import json



client = MongoClient("mongodb://127.0.0.1:27017")
#clinet = MongoClient("localhost:27017")
db = client.ppi

def layout():


    my_eda = html.Div([
                    html.Div([
                        html.Div([
                            html.P("Quantitative exploration", className="panel-title")
                            ],className="grid-item full-width"),
                        ], className="dashboard-grid"),
                    html.Div([
                            html.Div([
                                html.Div("Sample Summary", className="panel-title"),
                                html.Div(id='sample_info')
                            ], className="grid-item"),
                            html.Div([
                                html.Div("DEG Summary", className="panel-title"),
                                html.Div(id='deg_info')
                            ], className="grid-item"),
                        ], className="dashboard-grid"),
                    
                    html.Div([
                        html.Div([
                            html.Div("Fold Change Table", className="panel-title"),
                            dcc.Loading(html.Div(id = "my_deg_table")),
                            html.Div(html.Button(
                            children=[
                                html.I(className="fas fa-download", style={"marginRight": "8px"}),  # Font Awesome download icon
                                html.Span("Download DEGs"),  # Button text
                            ],
                            id="download-degs-button",
                            n_clicks=0,
                            style={
                                "backgroundColor": "white",  # White background
                                "color": "#333",  # Dark grey text color
                                "border": "2px solid #333",  # Dark grey border
                                "borderRadius": "8px",  # Rounded corners
                                "padding": "10px 20px",  # Padding
                                "fontSize": "16px",  # Font size
                                "boxShadow": "0px 4px 6px rgba(0, 0, 0, 0.1)",  # Subtle shadow
                                "cursor": "pointer",  # Pointer cursor
                                "display": "inline-flex",  # Align text and icon
                                "alignItems": "center",  # Center-align content
                            },
                        ),style={"textAlign": "left", "marginTop": "10px"}),
                        ], className="grid-item full-width"),
                    ], className="dashboard-grid"),
                     
                    
                    html.Div([
                            html.Div([
                                html.Div("Pathway exploration", className="panel-title"),
                                dcc.Dropdown(id="kegg_pathways", multi=False),
                                #html.Div(id='pathway_deg')
                                dcc.Loading(dcc.Graph(id='pathway_deg'))
                            ], className="grid-item"),
                            html.Div([
                                html.Div("Complex exploration", className="panel-title"),
                                dcc.Dropdown(id="protein_complex", multi=False),
                                #html.Div(id='complex_deg')
                                dcc.Loading(dcc.Graph(id='complex_deg'))
                            ], className="grid-item"),
                        ], className="dashboard-grid"),
                    html.Div([
                            html.Div([
                                html.Div("Similarity exploration", className="panel-title"),
                                dcc.Dropdown(id="gene_corr_query", multi=False),
                                # html.Div([
                                #     dcc.Input(id="corr_thres",value=0 ,type="number", placeholder="Pearson value (abs)"),
                                #     dcc.Input(id="corr_cat", type="text", placeholder="Enter extact set name"),
                                #     ],className='row'),
                                html.Div([
                                    html.Label("Pearson value (abs):", style={"margin-right": "10px", "font-weight": "bold"}),
                                    dcc.Input(
                                        id="corr_thres", 
                                        value=0, 
                                        type="number", 
                                        placeholder="Enter value",
                                        style={"margin-right": "20px", "padding": "5px", "width": "150px"}
                                    ),
                                    html.Label("Exact set name:", style={"margin-right": "10px", "font-weight": "bold"}),
                                    dcc.Input(
                                        id="corr_cat", 
                                        type="text", 
                                        placeholder="Enter name",
                                        style={"padding": "5px", "width": "150px"}
                                    )
                                ], style={"display": "flex", "align-items": "center", "gap": "20px", "margin-bottom": "15px"}),

                                
                                dcc.Loading(dcc.Graph(id='gene_corr_mat'))
                            ], className="grid-item"),
                            html.Div([
                                html.Div("Expression distribution", className="panel-title"),
                                dcc.Dropdown(id="gene_box_query", multi=True),
                                #html.Div(id='complex_deg')
                                dcc.Loading(dcc.Graph(id='gene_box')),
                                html.Div(id = 'gene_box_message')

                            ], className="grid-item"),
                        ], className="dashboard-grid"),

                    

                ])





    
    return my_eda


# # Helper Functions
# def generate_set_table(set1_name, set1_columns, set2_name, set2_columns):
#     """Generate an HTML table summarizing Set 1 and Set 2."""
#     return html.Table(
#         [
#             html.Thead(
#                 html.Tr([
#                     html.Th('Set', style={'background-color': '#d1e7dd', 'font-weight': 'bold'}),
#                     html.Th('Columns', style={'background-color': '#d1e7dd', 'font-weight': 'bold'})
#                 ])
#             ),
#             html.Tbody([
#                 html.Tr([
#                     html.Td(set1_name, style={'font-weight': 'bold', 'background-color': '#d4edda'}),
#                     html.Td(', '.join(set1_columns), style={'background-color': '#f1f9f5'})
#                 ]),
#                 html.Tr([
#                     html.Td(set2_name, style={'font-weight': 'bold', 'background-color': '#f8d7da'}),
#                     html.Td(', '.join(set2_columns), style={'background-color': '#f9ecee'})
#                 ])
#             ])
#         ],
#         style={'width': '100%', 'border': '1px solid black', 'padding': '8px', 'text-align': 'left', 'margin': '20px 0'}
#     )

def generate_set_table(meta_df):

    set_names = meta_df['Label'].unique()
    set1_name,  set2_name = set_names[0], set_names[1]
   

    set1_columns = meta_df[meta_df['Label'] == set1_name]['ID'].tolist()
    set2_columns = meta_df[meta_df['Label'] == set2_name]['ID'].tolist()

    return html.Table(
        [
            html.Thead(
                html.Tr([
                    html.Th('Set', style={'background-color': '#d1e7dd', 'font-weight': 'bold'}),
                    html.Th('Columns', style={'background-color': '#d1e7dd', 'font-weight': 'bold'})
                ])
            ),
            html.Tbody([
                html.Tr([
                    html.Td(set1_name, style={'font-weight': 'bold', 'background-color': '#d4edda'}),
                    html.Td(', '.join(set1_columns), style={'background-color': '#f1f9f5'})
                ]),
                html.Tr([
                    html.Td(set2_name, style={'font-weight': 'bold', 'background-color': '#f8d7da'}),
                    html.Td(', '.join(set2_columns), style={'background-color': '#f9ecee'})
                ])
            ])
        ],
        style={'width': '100%', 'border': '1px solid black', 'padding': '8px', 'text-align': 'left', 'margin': '20px 0'}
    )


def generate_significant_summary_table(significant_count, non_significant_count, total_genes):
    """Generate an HTML table summarizing the counts of significant and non-significant results."""
    return html.Table(
        [
            html.Thead(
                html.Tr([
                    html.Th('Count Type', style={'font-weight': 'bold'}),
                    html.Th('Count', style={'font-weight': 'bold'})
                ])
            ),
            html.Tbody([
                html.Tr([html.Td('Significant p-values'), html.Td(significant_count)]),
                html.Tr([html.Td('Non Significant p-values'), html.Td(non_significant_count)]),
                html.Tr([html.Td('Total Genes/Proteins'), html.Td(total_genes)])
            ])
        ],
        style={'width': '50%', 'border': '1px solid black', 'padding': '8px', 'text-align': 'left', 'margin': '20px 0'}
    )






def generate_deg_table(degs):
    """Generate a Dash DataTable for displaying significant DEGs."""

    #print('here is the table', degs.shape)
    return dash_table.DataTable(
        filter_action='native',
        data=degs.to_dict('records'),
        columns=[{"name": i, "id": i} for i in degs.columns],
        style_table={'height': '300px', 'overflowY': 'auto'},
        style_cell={'textAlign': 'center', 'padding': '8px'},
        style_header={'backgroundColor': '#f7f7f7', 'fontWeight': 'bold'},
        style_data_conditional=[
            {
                'if': {'column_id': 'Significant', 'filter_query': '{Significant} = True'},
                'backgroundColor': '#d4edda', 'color': 'black'
            }
        ]
    )

# def get_eda_knowledge():

#     pathway_db=db.kegg_pathway
#     path_res=my_db.find({},{'external_id':1,'pathway':1,'_id':0})
#     path_options = [{'label': name.capitalize(), 'value': path_id} for name in location if not pd.isna(name)]

#     complex_db=db.protein_complex
#     complex_res=my_db.find({},{'complex_id':1,'complex_name':1,'_id':0})
#     complex_options = [{'label': name.capitalize(), 'value': com_id} for name in location if not pd.isna(name)]

#     return path_options, complex_options

# import pandas as pd

def get_eda_knowledge():

    try:
        # Query the KEGG pathway database and convert to DataFrame
        pathway_db = db.kegg_pathway
        path_res = pd.DataFrame(list(pathway_db.find({}, {'external_id': 1, 'pathway': 1, '_id': 0})))
        #print('pathway_fetch',path_res.shape)

        # Process pathways into dropdown options
        if not path_res.empty:
            path_res.dropna(subset=['external_id', 'pathway'], inplace=True)
            path_options = path_res.apply(
                lambda row: {'label': row['pathway'].capitalize(), 'value': row['external_id']}, axis=1
            ).tolist()
        else:
            path_options = []

        # Query the protein complex database and convert to DataFrame
        complex_db = db.protein_complex
        complex_res = pd.DataFrame(list(complex_db.find({}, {'complex_id': 1, 'complex_name': 1, '_id': 0})))
        #print('--------complex--------')
        #print(complex_res)

        # Process complexes into dropdown options
        if not complex_res.empty:
            complex_res.dropna(subset=['complex_id', 'complex_name'], inplace=True)
            complex_options = complex_res.apply(
                lambda row: {'label': row['complex_name'].capitalize(), 'value': row['complex_id']}, axis=1
            ).tolist()
        else:
            complex_options = []

        #print(f'EDA_drop {len(path_options)} and {len(complex_options)}')

        return path_options, complex_options
        #a = [{'label': 'Folate biosynthesis - homo sapiens (human)', 'value': 'path:hsa00790'}, {'label': 'Complement and coagulation cascades - homo sapiens (human)', 'value': 'path:hsa04610'}, {'label': 'Citrate cycle (tca cycle) - homo sapiens (human)', 'value': 'path:hsa00020'}, {'label': 'Renin-angiotensin system - homo sapiens (human)', 'value': 'path:hsa04614'}]
        #return a, a

    except Exception as e:
        print(f"Error fetching or processing data from MongoDB: {e}")
        return [], []



################################for drop down########################################
# Helper function for generating bar plots

def fig_msg(my_title="No genes match your query."):

    fig = Figure()
    fig.update_layout(
        title=my_title,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        annotations=[
            dict(
                text="No Data Available",
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=16, color="red")
            )
        ]
    )
    return fig

def generate_bar_plot(sub_df, title):
    if sub_df.empty:
        # Return an empty figure with a "No Data" message
        fig = Figure()
        fig.update_layout(
            title="No genes match your query.",
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            annotations=[
                dict(
                    text="No Data Available",
                    xref="paper",
                    yref="paper",
                    showarrow=False,
                    font=dict(size=16, color="red")
                )
            ]
        )
        return fig

    # Create a bar plot for DEGs
    fig = px.bar(
        sub_df,
        x='Gene_id',
        y='log2_Fold_Change',
        color='p-value',
        color_continuous_scale='Viridis',
        labels={
            'Gene_id': 'Gene ID',
            'log2_Fold_Change': 'Log2 Fold Change',
            'p-value': 'P-value'
        },
        #title=title
    )

    # Customize figure layout
    fig.update_layout(
        xaxis_title="Gene ID",
        yaxis_title="Log2 Fold Change",
        coloraxis_colorbar_title="P-value",
        plot_bgcolor='white',
        margin=dict(l=40, r=40, t=40, b=40),
        paper_bgcolor='white',
        font=dict(size=12)
    )
    fig.update_xaxes(tickangle=45, showgrid=True, gridcolor='lightgrey')
    fig.update_yaxes(showgrid=True, gridcolor='lightgrey')

    # Add a horizontal line at y=0
    fig.add_shape(
        type="line",
        x0=-0.5,  # Extend slightly before the first bar
        x1=len(sub_df) - 0.5,  # Extend slightly after the last bar
        y0=0,
        y1=0,
        line=dict(color="black", width=1),
        xref="x",
        yref="y"
    )

    return fig

# Generic function for processing DEGs data
def process_deg_query(query_value, degs_data, db_collection, gene_field, query_field):
    if not degs_data:
        raise PreventUpdate

    df = pd.DataFrame.from_records(degs_data)

    # Query database and extract relevant gene names
    query_result = pd.DataFrame(list(db_collection.find({query_field: query_value}, {gene_field: 1, '_id': 0})))
    print('eda_query',query_result.shape)

    if query_result.empty or gene_field not in query_result.columns:
        return pd.DataFrame()  # Return empty DataFrame for no matches

    # Flatten gene list if comma-separated
    if query_field=='external_id': 
        sep = ','
    else:
        sep=';'

    query_result[gene_field] = query_result[gene_field].apply(
            lambda x: x.split(sep) if isinstance(x, str) else []
        )
    gene_list = query_result[gene_field].explode().str.strip().dropna().unique()

    # Filter DEGs data for matching genes
    return df[df['Gene_id'].isin(gene_list)]

def get_gene_options(my_list):

    options = [{'label': name.capitalize(), 'value': name} for name in my_list if not pd.isna(name)]

    return options

# Callback to update Set 2 dropdown dynamically
@callback(
    Output('set2-dropdown', 'value'),
    Output('set2-dropdown', 'options'),
    Input('set1-dropdown', 'value'),
    State('upload_status', 'data'), 
    prevent_initial_call=False
)
def update_set2_columns(set1_columns, status_data):
    if not status_data:
        return [], []

    file_path = status_data['uploaded_files'][0]
    if file_path:
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        if file_path.endswith('.tsv'):
            df = pd.read_csv(file_path,sep="\t")
        if file_path.endswith('.xls'):
            df = pd.read_excel(file_path)
   
    total_columns = df.columns[1:]  # All columns except the first one

    if set1_columns:
        set2_columns = [col for col in total_columns if col not in set1_columns]
    else:
        set2_columns = total_columns

    return set2_columns, [{'label': col, 'value': col} for col in set2_columns]

# Callback for Submit button to process data and store DEGs in dcc.Store
# Main Callback

@callback(
    Output('my_deg_table', 'children'),
    Output('sample_info', 'children'),
    Output('deg_info', 'children'),
    Output('degs-store', 'data'),
    Output('data_label', 'data'),
    Output('kegg_pathways', 'options'),
    Output('protein_complex', 'options'),
    Output('gene_corr_query', 'options'),
    Output('gene_box_query', 'options'),
    Output('n_clicks_store', 'data'),  # To track the previous value of n_clicks
    Input('submit-columns', 'n_clicks'),
    State('n_clicks_store', 'data'),  # Store for previous n_clicks
    State('set1-dropdown', 'value'),
    State('set1-name', 'value'),
    State('set2-name', 'value'),
    State('norm-method-dropdown', 'value'),
    State('scaling-method-dropdown', 'value'),
    State('alpha-input', 'value'),
    State('upload_status', 'data'),
    State('degs-store', 'data'),
    State('data_label', 'data'),
    prevent_initial_call=True
)
def process_column_selection(n_clicks, prev_n_clicks, set1_columns, set1_name, set2_name, norm_method, scaling_method, alpha, status_data, degs_store, data_label):
    # Check if the callback is triggered by a change in n_clicks
    
    if not status_data:
        raise PreventUpdate

    ctx = dash.callback_context
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # Fetch knowledge base for dropdowns
    path_list, complex_list = get_eda_knowledge()

    # Handle cases where existing data is reused
    if (prev_n_clicks is not None and n_clicks == prev_n_clicks and degs_store and data_label) or ( n_clicks==0 and degs_store and data_label):
        print("Reusing existing data")
        degs = pd.DataFrame(degs_store)
        meta_data = pd.DataFrame(data_label)
        significant_degs = degs[degs['Significant'] == True]
        gene_option = get_gene_options(degs['Gene_id'].to_list())

        return (
            generate_deg_table(significant_degs),
            generate_set_table(meta_data),
            generate_significant_summary_table(len(significant_degs), len(degs) - len(significant_degs), len(degs)),
            dash.no_update,
            dash.no_update,
            path_list,
            complex_list,
            gene_option,
            gene_option,
            n_clicks
        )

    # File and metadata processing
    file_path = status_data['uploaded_files'][0]
    if file_path:
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        if file_path.endswith('.tsv'):
            df = pd.read_csv(file_path,sep="\t")
        if file_path.endswith('.xls'):
            df = pd.read_excel(file_path)

    if not set1_name or not set2_name:
        return html.Div("Please provide names for both sets.", style={'color': 'red'}), None, None, None, None, [], [], [], [], n_clicks
    if not set1_columns:
        return html.Div("Please select columns for the first set.", style={'color': 'red'}), None, None, None, None, [], [], [], [], n_clicks
    print(n_clicks)
    # Process columns for grouping
    total_columns = df.columns[1:]
    set2_columns = [col for col in total_columns if col not in set1_columns]
    meta_data = pd.DataFrame({'ID': set1_columns + set2_columns, 'Label': [set1_name] * len(set1_columns) + [set2_name] * len(set2_columns)})

    # Perform DEG analysis
    degs = preprocess_and_find_DEGs(df, meta_data, norm_method=norm_method, scaling_method=scaling_method, alpha=alpha)
    significant_degs = degs[degs['Significant'] == True]
    degs_store = degs.to_dict('records')
    gene_option = get_gene_options(degs['Gene_id'].to_list())

    # Fetch EDA knowledge
    path_list, complex_list = get_eda_knowledge()

    return (
        generate_deg_table(significant_degs),
        generate_set_table(meta_data),
        generate_significant_summary_table(len(significant_degs), len(degs) - len(significant_degs), len(degs)),
        degs_store,
        meta_data.to_dict('records'),
        path_list,
        complex_list,
        gene_option,
        gene_option,
        n_clicks  # Update the n_clicks_store with the new n_clicks
    )



# Download DEGs data from dcc.Store
@callback(
    Output('download-degs', 'data'),
    Input('download-degs-button', 'n_clicks'),
    State('degs-store', 'data'),
    prevent_initial_call=True
)
def download_degs(n_clicks, degs_data):

    if n_clicks > 0 and degs_data:
        df = pd.DataFrame.from_records(degs_data)
        return dcc.send_data_frame(df.to_csv, "DEGs.csv")



# Callback for pathway DEGs
@callback(
    Output('pathway_deg', 'figure'),
    Input('kegg_pathways', 'value'),
    State('degs-store', 'data'),
    prevent_initial_call=True
)
def download_pathway_degs(path_query, degs_data):
    sub_df = process_deg_query(
        query_value=path_query,
        degs_data=degs_data,
        db_collection=db.kegg_pathway,
        gene_field='hgnc_symbol_ids',
        query_field='external_id'
    )
    return generate_bar_plot(sub_df, "Differentially Expressed Genes (DEGs) for Selected Pathway")

# Callback for complex DEGs
@callback(
    Output('complex_deg', 'figure'),
    Input('protein_complex', 'value'),
    State('degs-store', 'data'),
    prevent_initial_call=True
)
def download_complex_degs(complex_query, degs_data):
    sub_df = process_deg_query(
        query_value=complex_query,
        degs_data=degs_data,
        db_collection=db.protein_complex,
        gene_field='subunits_gene_name',
        query_field='complex_id'
    )
    return generate_bar_plot(sub_df, "Differentially Expressed Genes (DEGs) for Selected Protein Complex")


# Callback for complex DEGs
@callback(
    Output('gene_corr_mat', 'figure'),
    Input('gene_corr_query', 'value'),
    Input('corr_thres','value'),
    Input('corr_cat','value'),
    State('data_label', 'data'),
    State('upload_status', 'data'),
    prevent_initial_call=True
)
def get_similar_genes(query_gene, thres_corr, data_cat,data_label, status_data):
    thres_pval = 0.01
    #thres_corr = 0

    if not data_label and not status_data:
        raise PreventUpdate

    if not query_gene or not data_cat:

        return fig_msg('Select gene from dropdown and fill data category') 

    print('data query',data_cat)
    # if query_gene not in expression_data.index:
    #     return generate_bar_plot(pd.DataFrame(), "No genes match your query.")  # Return empty plot
    
    # Calculate correlation and p-value for each gene
    # Process file and generate DEGs
    meta_df = pd.DataFrame(data_label)
    selected_columns = list(meta_df.loc[meta_df['Label']==data_cat,'ID'])
    print('new check columns selected ', len(selected_columns))
    file_path = status_data['uploaded_files'][0]
    
    if file_path:
        if file_path.endswith('.csv'):
            expression_data= pd.read_csv(file_path,index_col=0)
        if file_path.endswith('.tsv'):
            expression_data= pd.read_csv(file_path,sep="\t" ,index_col=0)
        if file_path.endswith('.xls'):
            expression_data=pd.read_excel(file_path,index_col=0)
     
    correlations = []
    expression_sub_df = expression_data[selected_columns]
    query_values = expression_data.loc[query_gene,selected_columns].values
    correlations = []
    #query_values = expression_data.loc[query_gene,group_labels['Label']==data_cat].values
    
    print(query_values.shape)
    for gene in expression_data.index:
        if gene != query_gene:
            gene_values = expression_sub_df.loc[gene].values
            print(gene)
            r, p = pearsonr(query_values, gene_values)
            correlations.append({'Gene_id': gene, 'Correlation': r, 'p-value': p})
    
    # Create a DataFrame for results
    corr_df = pd.DataFrame(correlations)

    # Filter for significant correlations
     # Apply FDR correction to p-values
    _, fdr_pvals, _, _ = multipletests(corr_df['p-value'], method='fdr_bh')
    corr_df['FDR_p-value'] = fdr_pvals

    # Filter for significant correlations based on FDR threshold
    significant_corr_df = corr_df[(corr_df['FDR_p-value'] < thres_pval)&(corr_df['Correlation'].abs() > thres_corr)]
    #significant_corr_df = corr_df[corr_df['p-value'] < threshold]

    # Sort by correlation strength
    significant_corr_df = significant_corr_df.sort_values(by='Correlation', ascending=False)

    # Generate a bar plot for significant correlations
    return generate_bar_plot(
        significant_corr_df.rename(columns={'Correlation': 'log2_Fold_Change'}),
        "Genes Significantly Correlated with {}".format(query_gene)
    )


# Callback for complex DEGs
@callback(
    Output('gene_box', 'figure'),
    Output('gene_box_message','children'),
    Input('gene_box_query', 'value'),
    State('degs-store', 'data'),
    State('data_label','data'),
    State('upload_status', 'data'),
    prevent_initial_call=True
)
def get_feature_analysis(selected_features, deg_data, data_labels,status_data):
    """
    Generates a box plot for selected features and evaluates their predictive capability.

    Parameters:
        selected_features: List of feature names to analyze.
        expression_data: DataFrame with features in rows and samples in columns.
        meta_data: DataFrame with 'ID' (sample index) and 'Label' (group label).

    Returns:
        - Plotly box plot figure if valid features are selected, else an empty figure with a message.
        - Text containing the predictive capability metrics (accuracy and AUROC).
    """
    if not deg_data and not status_data:
        raise PreventUpdate
    # Validate feature selection
    if len(selected_features) > 5:
        return px.scatter(), "Please select 5 or fewer features for analysis."

    if len(selected_features) < 1:
        return px.scatter(), "Please select at least 1 feature for analysis."

    # Ensure meta_data has the required columns.
    meta_data = pd.DataFrame(data_labels)
    if not {'ID', 'Label'}.issubset(meta_data.columns):
        raise ValueError("Metadata is missing required columns: 'ID' and 'Label'.")

    #loading expression data
    file_path = status_data['uploaded_files'][0]
    if file_path:
        if file_path.endswith('.csv'):
            expression_data= pd.read_csv(file_path,index_col=0)
        if file_path.endswith('.tsv'):
            expression_data= pd.read_csv(file_path,sep="\t" ,index_col=0)
        if file_path.endswith('.xls'):
            expression_data=pd.read_excel(file_path,index_col=0)

    # Align expression_data with metadata
    
    meta_data = meta_data.set_index('ID')
    common_samples = expression_data.columns.intersection(meta_data.index)
    expression_data = expression_data.loc[selected_features, common_samples]
    meta_data = meta_data.loc[common_samples]

    if expression_data.empty:
        return px.scatter(), "No matching features or samples found for analysis."

    # Box plot
    box_data = expression_data.T  # Transpose for box plot
    box_data['Group'] = meta_data['Label']

    fig = px.box(
        box_data,
        #x='Group',
        #y=selected_features[0] if len(selected_features) == 1 else None,
        points="all",
        color='Group',
        title="Box Plot of Selected Features",
        labels={'value': 'Expression', 'Group': 'Sample Group'},
    )
    fig.update_layout(

        plot_bgcolor="white",
        paper_bgcolor="white"
        )

    fig.update_xaxes(showgrid=True, gridcolor='lightgrey',tickangle=45)
    fig.update_yaxes(showgrid=True, gridcolor='lightgrey')

    # 5-Fold Cross-Validation for Predictive Capability
    X = expression_data.T  # Features in columns, samples in rows
    y = meta_data['Label']
    classifier = LogisticRegression(max_iter=1000)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    accuracy_scores = []
    auroc_scores = []

    for train_idx, test_idx in cv.split(X, y):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        classifier.fit(X_train, y_train)
        y_pred = classifier.predict(X_test)
        y_pred_prob = classifier.predict_proba(X_test)[:, 1]

        accuracy_scores.append(accuracy_score(y_test, y_pred))
        auroc_scores.append(roc_auc_score(y_test, y_pred_prob))

    mean_accuracy = np.mean(accuracy_scores)
    mean_auroc = np.mean(auroc_scores)

    metrics_text = (
        f"Predictive Capability of Selected Features:\n"
        f" - Mean Accuracy: {mean_accuracy:.2f}\n"
        f" - Mean AUROC: {mean_auroc:.2f}"
    )

    return fig, metrics_text







# def download_degs(path_query, degs_data):

#     if not degs_data:
#         raise PreventUpdate
#     df = pd.DataFrame.from_records(degs_data)


#     pathway_db = db.kegg_pathway
#     path_res = pd.DataFrame(list(pathway_db.find({'external_id':path_query}, {'hgnc_symbol_ids': 1, '_id': 0})))

#     path_genes = path_res['hgnc_symbol_ids']

#     sub_df = df.loc[df['Gene_id'].isin(path_genes)]# improve that there might be gene not prsent in df
#     fig = #create a bar plot where Gene_id will be in x axis log_2_Fold_Change will be on Y axis and p-value will be reflected as color of bar

#     return fig
