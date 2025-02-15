import pandas as pd
from dash import dcc, html, callback, Input, Output, State
import plotly.express as px
import plotly.figure_factory as ff
import numpy as np
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
from io import BytesIO
import base64
from utils.get_plots import get_enrichment_data, plot_enrichment_bubble
from utils.panel import *








# Layout function with PCA and Hclust Heatmap panels
def layout():
    my_viz = html.Div([
                    html.Div([
                        html.Div([
                            html.P("Qualitative exploration", className="panel-title")
                            ],className="grid-item full-width"),
                        ], className="dashboard-grid"),

                    html.Div([
                            html.Div([
                                html.Div("Volcano plot", className="panel-title"),
                                # dcc.Slider(
                                #         id='fold-change-threshold',
                                #         min=0,
                                #         max=5,
                                #         step=0.1,
                                #         value=0.0,
                                #         marks={i: f'{i}' for i in range(6)},
                                #         tooltip={"placement": "bottom", "always_visible": True}
                                #     ),
                                dcc.Loading(dcc.Graph(id="volcano-plot")),
                                html.Div(id='error-volcano', style={'color': 'red'})
                            ], className="grid-item"),
                            html.Div([
                                html.Div("Top 10 enrichment", className="panel-title"),
                             #    dcc.Dropdown(id='enrich_cat',clearable=False,
                             #    options=[{'label': 'Process','value': 'GO_Biological_Process_2018'},{'label':'Pathway','value':'KEGG_2019_Human'},{'label':'Disease','value':'ClinVar_2019'}],
                             #    value = 'KEGG_2019_Human'
                             # ),
                                dcc.Loading(dcc.Graph(id="pathway-plot")),
                                html.Div(id='error-pathway', style={'color': 'red'})
                            ], className="grid-item"),
                        ], className="dashboard-grid"),

                    html.Div([
                            html.Div([
                                html.Div("PCA (All Features)", className="panel-title"),
                                dcc.Loading(dcc.Graph(id="pca-plot1"))
                                
                            ], className="grid-item"),
                            html.Div([
                                html.Div("PCA (Selected Features)", className="panel-title"),
                                dcc.Loading(dcc.Graph(id="pca-plot2")),
                            ], className="grid-item"),
                        ], className="dashboard-grid")
                    ])

        
    
    
    return my_viz

# Callback for generating the Volcano Plot
@callback(
    [Output('volcano-plot', 'figure'),
     Output('error-volcano', 'children')],
    [Input('fold-change-threshold', 'value'),
     Input('degs-store', 'data')]
)
def generate_volcano_plot(fold_change_threshold, degs_store):
    try:
        if degs_store is None:
            return px.scatter(), "No DEGs data found."

        # Convert DEGs data to DataFrame
        degs = pd.DataFrame.from_records(degs_store)

        # Validate required columns
        required_columns = ['Fold Change', 'p-value', 'Significant']
        if not all(col in degs.columns for col in required_columns):
            return px.scatter(), "DEGs data is missing required columns."

        # Compute log2 fold change and -log10(p-value)
        degs['log2_Fold_Change'] = np.log2(degs['Fold Change'].replace(0, np.nan))
        degs['-log10_p-value'] = -np.log10(degs['p-value'].replace(0, np.nan))

        # Add a color column based on thresholds
        degs['Highlight'] = degs.apply(
            lambda row: 'Upregulated' if (row['log2_Fold_Change'] >= fold_change_threshold) & (row['Significant'] == True) else
                ('Downregulated' if (row['log2_Fold_Change'] <= -fold_change_threshold) & (row['Significant'] == True) else
                 'Non-Significant'),
            axis=1
        )

        # Generate the volcano plot
        volcano_fig = px.scatter(
            degs,
            x='log2_Fold_Change',
            y='-log10_p-value',
            color='Highlight',
            color_discrete_map={'Upregulated': 'blue', 'Downregulated': 'red', 'Non-Significant': 'grey'},
            #title=f"(Threshold: |Log2 Fold Change| ≥ {fold_change_threshold})",
            labels={'log2_Fold_Change': 'Log2 Fold Change', '-log10_p-value': '-log10(p-value)'},
            hover_data={
                'Fold Change': ':.2f',
                'p-value': ':.2e',
                'log2_Fold_Change': ':.2f',
                '-log10_p-value': ':.2f'
            }
        )

        # Improve layout
        volcano_fig.update_traces(marker=dict(size=8, opacity=0.7))
        volcano_fig.update_layout(
            template='plotly_white',
            xaxis=dict(title='Log2 Fold Change'),
            yaxis=dict(title='-log10(p-value)'),
            legend_title=dict(text='Significance'),
        )

        return volcano_fig, ""

    except Exception as e:
        return px.scatter(), f"An error occurred: {str(e)}"

@callback(
    [Output('pathway-plot', 'figure'),
     Output('error-pathway', 'children')],
    [Input('fold-change-threshold', 'value'),
     Input('degs-store', 'data'),
     Input('enrich_cat','value')]
     #State('enrich_cat','label')]
)
def generate_enrichment_plot(fold_change_threshold, degs_store,data_name):
    #data_name = 'KEGG_2019_Human'
    #my_prots = [ 'RB1', 'PPP1CA', 'TP53', 'CSNK2A1', 'CDK1', 'CHEK1', 'EEF2K', 'EGFR', 'ERBB2', 'CDC7', 'AR', 'BRCA1', 'MAPK6', 'SIRT1', 'NME1', 'EIF2AK2' ]
    # ctx = dash.callback_context
    # ctx_msg = json.dumps({
            
    #         'triggered': ctx.triggered,
    #         'inputs': ctx.inputs,
    #         'state' : ctx.states
    #     }, indent=2)

    # print('enriching',ctx_msg)

    try:
        if degs_store is None:
            return px.scatter(), "No DEGs data found."

        degs = pd.DataFrame.from_records(degs_store)

        if 'Fold Change' not in degs.columns or 'p-value' not in degs.columns or 'Significant' not in degs.columns:
            return px.scatter(), "DEGs data is missing required columns."

        degs['log2_Fold_Change'] = np.log2(degs['Fold Change'].replace(0, np.nan)).dropna()
        filtered_degs = degs.loc[(degs['log2_Fold_Change'].abs() >= fold_change_threshold) & (degs['Significant']==True)]

        


        gene_list = filtered_degs['Gene_id'].to_list()

        # Fetch enrichment data
        len_deg = len(gene_list)
        # if len_deg <100:
        #     print('here is deg')
        #     print(gene_list)
        # else:
        #     print('more than 100 deg')
        #     print(len(gene_list))
        my_prots = [gene.split('/')[0].strip() for gene in gene_list if isinstance(gene, str) and gene.strip()]
        enrich_table, err = get_enrichment_data(my_prots, data_name)

        # Ensure there was no error before plotting
        if not err:
            #print("here is the enrichment",enrich_label)
            if 'Clin' in data_name:
                enrich_label = 'Disease'
            elif 'Bio' in data_name:
                enrich_label = 'Process'
            else:
                enrich_label = 'Pathway'
            enrichment_fig = plot_enrichment_bubble(enrich_table[:10],enrich_label)
            #fig.show()
        else:
            print(f"Error fetching enrichment data: {err}")
        

        return enrichment_fig, ""

    except Exception as e:
        return px.scatter(), f"An error occurred: {str(e)}"





# Callback for PCA Plot
@callback(
    Output('pca-plot1', 'figure'),
    Input('degs-store', 'data'),
    State('upload_status','data'),
    State('data_label','data')
)
def generate_pca_plot(degs_store,status_data,data_labels):
    if degs_store is None:
        return px.scatter()


    file_path = status_data['uploaded_files'][0]
    if file_path:
        if file_path.endswith('.csv'):
            my_df= pd.read_csv(file_path)
        if file_path.endswith('.tsv'):
            my_df= pd.read_csv(file_path,sep="\t")
        if file_path.endswith('.xls'):
            my_df=pd.read_excel(file_path)
    meta_data = pd.DataFrame(data_labels)
    #df = temp_df.drop(temp_df.columns[0],axis=1)

    sample_ids = list(my_df.columns[1:])  # Exclude the index column
    group_labels = meta_data.set_index('ID').loc[sample_ids, 'Label'].values

    # Perform PCA
    df = my_df.drop(my_df.columns[0],axis=1)
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(df.T)  # Transpose to make samples rows and features columns

    #group1 = data[:, np.array(group_labels) == unique_groups[0]]
    #group2 = data[:, np.array(group_labels) == unique_groups[1]]
    
    # Create a DataFrame for PCA results
    pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
    pca_df['Group'] = group_labels
    unique_groups = np.unique(group_labels)
    set1_name, set2_name = unique_groups[0], unique_groups[1]
    
    # Create scatter plot
    fig = px.scatter(
        pca_df, 
        x='PC1', 
        y='PC2', 
        color='Group', 
        #title="PCA Plot of Gene Expressions",
        labels={'PC1': 'Principal Component 1', 'PC2': 'Principal Component 2'},
        color_discrete_map={set1_name: 'blue', set2_name: 'orange'}
    )
    fig.update_traces(marker=dict(size=10, line=dict(width=1, color='DarkSlateGrey')))
    fig.update_layout(
        template='plotly_white',
        xaxis=dict(title='PC1'),
        yaxis=dict(title='PC2'),
        legend_title=dict(text='Groups')
    )

    return fig


# Callback for DEG PCA Plot
@callback(
    Output('pca-plot2', 'figure'),
    Input('degs-store', 'data'),
    Input('fold-change-threshold', 'value'),
    State('upload_status','data'),
    State('data_label','data')
)
def generate_pca_plot(degs_store, fold_change_threshold, status_data, data_labels):
    if degs_store is None:
        return px.scatter()

    #file_path = status_data['uploaded_files'][0]
    #df = pd.read_csv(file_path,index_col=0)

    file_path = status_data['uploaded_files'][0]
    if file_path:
        if file_path.endswith('.csv'):
            my_df= pd.read_csv(file_path,index_col=0)
        if file_path.endswith('.tsv'):
            my_df= pd.read_csv(file_path,sep="\t",index_col=0)
        if file_path.endswith('.xls',index_col=0):
            my_df=pd.read_excel(file_path)
    
    meta_data = pd.DataFrame(data_labels)
    #df = temp_df.drop(temp_df.columns[0],axis=1)

    sample_ids = list(my_df.columns)  # Exclude the index column
    group_labels = meta_data.set_index('ID').loc[sample_ids, 'Label'].values

    

    # DEG data
    degs = pd.DataFrame.from_records(degs_store)

    # Validate required columns
    required_columns = ['Fold Change', 'p-value', 'Significant']
    if not all(col in degs.columns for col in required_columns):
        return px.scatter(), "DEGs data is missing required columns."

    # Compute log2 fold change and -log10(p-value)
    degs['log2_Fold_Change'] = np.log2(degs['Fold Change'].replace(0, np.nan))
    degs['-log10_p-value'] = -np.log10(degs['p-value'].replace(0, np.nan))

    # Filter only significant features based on the fold change threshold and significance status
    significant_degs = degs[(abs(degs['log2_Fold_Change']) >= fold_change_threshold) & (degs['Significant'] == True)]

    # Get the significant gene names
    significant_genes = significant_degs['Gene_id'].tolist()
    significant_genes = [gene for gene in significant_genes if gene is not None]

    # Filter the df to keep only the rows with significant genes
    #df = my_df.drop(my_df.columns[0],axis=1)
    filtered_df = my_df.loc[significant_genes, :]  # Keep only rows (features) corresponding to significant genes
    filtered_df = filtered_df.reset_index(drop=True)

    if filtered_df.empty:
        return px.scatter(), "No significant genes found after applying the thresholds."

    # Perform PCA on the filtered data (samples are in columns)
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(filtered_df.T)  # Transpose to make samples rows and features columns

    # Create a DataFrame for PCA results
    pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
    pca_df['Group'] = group_labels
    unique_groups = np.unique(group_labels)
    set1_name, set2_name = unique_groups[0], unique_groups[1]

    # Create scatter plot
    fig = px.scatter(
        pca_df, 
        x='PC1', 
        y='PC2', 
        color='Group', 
        #title="PCA Plot of Significant Gene Expressions",
        labels={'PC1': 'Principal Component 1', 'PC2': 'Principal Component 2'},
        color_discrete_map={set1_name: 'blue', set2_name: 'orange'}
    )
    fig.update_traces(marker=dict(size=10, line=dict(width=1, color='DarkSlateGrey')))
    fig.update_layout(
        template='plotly_white',
        xaxis=dict(title='PC1'),
        yaxis=dict(title='PC2'),
        legend_title=dict(text='Groups')
    )

    return fig




