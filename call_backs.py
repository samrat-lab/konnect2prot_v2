from pages import  Data, Visualization, NetworkAnalysis, EDA
from pages.my_sidebars import get_data_sidebar, get_net_sidebar, get_eda_sidebar, get_vis_sidebar, get_tab1_sidebar
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash import html, dcc, dash_table, get_asset_url
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
import os
import json
import dash
import requests
import base64



example_prots='RB1\nPPP1CA\nTP53\nCSNK2A1\nCDK1\nCHEK1\nEEF2K\nEGFR\nERBB2\nCDC7\nAR\nBRCA1\nMAPK6\nSIRT1\nNME1\nEIF2AK2'

client = MongoClient("mongodb://127.0.0.1:27017")
#clinet = MongoClient("localhost:27017")
db = client.ppi
def get_func_info(gene,trim_str='Function [CC]\nFUNCTION: '):
    
    #db=client.ppi_test
    my_db=db.database1
    res=my_db.find({'name':gene.upper()},{'id':1,'name':1,'_id':0})
    #print('tap test', res[0])
    u_id=res[0]['id']
    g_id=res[0]['name']
    
    url=f'https://rest.uniprot.org/uniprotkb/{u_id}?format=json&fields=cc_function'
    my_hnd=requests.get(url)
    url_json=my_hnd.json()
    #return url_txt
    page_txt=url_json['comments'][0]['texts'][0]['value']
    #print('uniprot connection check', page_txt)
    if page_txt:
        func_info=my_hnd.text.split(' (')[0][len(trim_str):]+'.'
        flines=page_txt.split('.')
        func_info=flines[0]+'.'+flines[-1]+'.'
        #func_info=page_txt
    else:
        func_info='Not Found'
    return g_id,u_id,func_info


def query_db_net(prots):
    query=prots.copy()
    if prots:
    #print(prots)
    #my_res=data.find({'name':{'$in':prots}});
    #db = client.ppi
        my_data = db.edge_lists2
        my_res = my_data.find({'$or':[{ "ENTITYA": {'$in':prots }},{ "ENTITYB": {'$in':prots }}]},{'ENTITYA':1,'ENTITYB':1})
        cnt=0
        prots=[]
        for res in my_res:
            cnt+=1
            prots.append(res['ENTITYA'])
            prots.append(res['ENTITYB'])
            #print(res)
    prots=list(set(prots))
    my_res = my_data.find({'$and':[{ "ENTITYA": {'$in':prots }},{ "ENTITYB": {'$in':prots }}]},{'ENTITYA':1,'ENTITYB':1,'EFFECT':1})
    
    g = nx.DiGraph()
    colors = ['red', 'blue', 'green', 'yellow', 'pink']
    cy_nodes=[]
    cy_edges=[]
    #disease=[]
    nodes=set()
    for res in my_res:

        #disease=disease+res['Disease']
        source = res['ENTITYA']
        target = res['ENTITYB']
        
        if 'up' in res['EFFECT']:
                color='blue'
        elif 'down' in res['EFFECT']:
            color ='red'
        else: 
            color ='grey'
        
        cy_edges.append({  # Add the Edge Node
                'data': {'source': source, 'target': target},
                'classes': color
            })
        if source in query:
            cat='search'
        else:
            cat='not_search'
        
        if source not in nodes:  # Add the source node
            nodes.add(source)
            cy_nodes.append({"data": {"id": source,'label':source}, "classes": cat})

        if target in query:
        #print(target)
            cat='search'
        else:
            cat='not_search'

        if target not in nodes:  # Add the target node
            nodes.add(target)
            cy_nodes.append({"data": {"id": target,'label':target }, "classes": cat})
    return cy_edges,cy_nodes,nodes


def fetch_dropdown(nodes):
    location = []
    function = []
    process = []
    pathway = []
    disease = []
    tissue_exp = []

    db = client.ppi
    my_data = db.database1
    my_res = my_data.find({ "name": {'$in': list(nodes) }}, 
                          {'location': 1, 'biological_proc': 1, 'molecular_func': 1, 
                           'pathway_class': 1, 'diseases': 1, 'tissue_exp': 1, '_id': 0})

    for res in my_res:
        location.extend(res.get('location', []))
        function.extend(res.get('molecular_func', []))
        process.extend(res.get('biological_proc', []))
        pathway.extend(res.get('pathway_class', []))
        tissue_exp.extend(res.get('tissue_exp', []))
        disease.extend(res.get('diseases', []))

    location = list(set(location))
    function = list(set(function))
    process = list(set(process))
    pathway = list(set(pathway))
    disease = list(set(disease))
    tissue_exp = list(set(tissue_exp))

    options = [{'label': name.capitalize(), 'value': name} for name in location if not pd.isna(name)]
    options1 = [{'label': name.capitalize(), 'value': name} for name in function if not pd.isna(name)]
    options2 = [{'label': name.capitalize(), 'value': name} for name in process if not pd.isna(name)]
    options3 = [{'label': name.capitalize(), 'value': name} for name in pathway if not pd.isna(name)]
    options4 = [{'label': name.capitalize(), 'value': name} for name in disease if not pd.isna(name)]
    options5 = [{'label': str(name), 'value': name} for name in tissue_exp if not pd.isna(name)]

    return options, options1, options2, options3, options4, options5

# Helper function: Fetch data from the database
def fetch_data_from_db(collection, filter_, projection):
    my_db = db[collection]
    return list(my_db.find(filter_, projection))

# Helper function: Generate a DataTable
def generate_data_table(data, columns, width="150px", page_size=5):
    return html.Div(
        [
            dash_table.DataTable(
                data=data,
                columns=[{'name': col, 'id': col} for col in columns],
                fixed_rows={'headers': True, 'data': 0},
                style_data={'border': '1px solid black'},
                style_cell={'width': width},
                page_size=page_size,
                export_format='xlsx',
            )
        ],
        style={'height': '100%', 'width': '100%'}
    )

def uploadstatus_to_dict(status: du.UploadStatus) -> dict:
    return {
        'uploaded_files': [str(file) for file in status.uploaded_files],  # convert Path to str
        'n_total': status.n_total,
        'uploaded_size_mb': status.uploaded_size_mb,
        'total_size_mb': status.total_size_mb,
        'upload_id': status.upload_id,
        'is_completed': status.is_completed,
        'n_uploaded': status.n_uploaded
    }

def dict_to_uploadstatus(data: dict) -> du.UploadStatus:
    return du.UploadStatus(
        uploaded_files=data['uploaded_files'],
        n_total=data['n_total'],
        uploaded_size_mb=data['uploaded_size_mb'],
        total_size_mb=data['total_size_mb'],
        upload_id=data['upload_id']
    )


# New multiquery function
def multi_filter_query(loc,fun,proc,path,dis,exp):
    #db=client.ppi
    my_db = db.database1
    #Updating location query
    if loc:
        locate=loc
        print('query',locate)
    else:
        locate={'$exists': True}
    #Updating function query
    if fun:
        mol_fun=fun
        print('query',mol_fun)
    else:
        mol_fun={'$exists': True}
    #Updating process query
    if proc:
        bio_proc=proc
        print('query',bio_proc)
    else:
        bio_proc={'$exists': True}
    #Updating pathway query
    if path:
        pathway=path
        print('query',pathway)
    else:
        pathway={'$exists': True}
    #Updating disease query
    if dis:
        disease=dis
        print('query',disease)
    else:
        disease={'$exists': True}
    #Updating expression query
    if exp:
        tiss_exp=exp
        print('query',tiss_exp)
    else:
        tiss_exp={'$exists': True}
    
        
    
    
    my_res=my_db.find(
        {
            '$and': [ 
                { 'location':locate },
                { 'molecular_func':mol_fun },
                { 'biological_proc':bio_proc },
                { 'pathway_class':pathway },
                { 'diseases':disease },
                { 'tissue_exp': tiss_exp }
            ]
        },
        {'_id':0,'name':1}
    )
    
    tmp_df=pd.DataFrame(list(my_res))
    #print('Cant you see',tmp_df)
    if tmp_df.empty:
        return []
    else:
        return tmp_df['name'].to_list()


def get_current_valid_edges(current_nodes, all_edges):
    """Returns edges that are present in Cytoscape:
    its source and target nodes are still present in the graph.
    """
    valid_edges = []
    node_ids = {n['data']['id'] for n in current_nodes}

    for e in all_edges:
        if e['data']['source'] in node_ids and e['data']['target'] in node_ids:
            valid_edges.append(e)
    return valid_edges

def update_nodes(el_tmp,ul):

    el_copy=el_tmp.copy()
    #valid_edges = []
    current_nodes=[]
    all_edges=[]
    cnt=0
    for itm in el_copy:
        # if the element is a node
        if 'source' not in itm['data']:
            if itm['data']['label'] in ul:
                current_nodes.append(itm)
        else:
            all_edges.append(itm)
    val_edges = get_current_valid_edges(current_nodes, all_edges)
    return val_edges+current_nodes

def cyto_2_graph(con_elements):
    g=nx.DiGraph()
    my_sub_node=get_nodes(con_elements)
    if len(con_elements)>18:
        cnt=0
        for el in con_elements:

            try:
                g.add_edge(el['data']['source'],el['data']['target'])
            except Exception as e:
                g.add_node(el['data']['id'])
                continue

    g_sub=g.subgraph(my_sub_node)
    return g_sub

def get_nodes(el):
    nodes=[]
    for itm in el:
        try:
            nodes.append(itm['data']['label'])
        except Exception as e:
            continue
    return nodes

def get_edges(el):
    edges=[]
    for itm in el:
        try:
            tmp=itm['data']['source']
            edges.append(itm['data'])
        except Exception as e:
            continue
    return edges

def query_signa_net(top_prots):
    #query=top_prots.copy()
    my_data = db.signaling
    my_res=my_data.find({ 'source_name':{'$in':top_prots} },{'source_name':1,'source_pathways':1,'_id':0})
    temp=list(my_res)
    network_df=pd.DataFrame(temp)
    sig_paths=list(network_df['source_pathways'].unique())
    my_sig_prot=network_df.groupby(['source_pathways']).size()
    my_res=my_data.find({ 'source_pathways':{'$in':sig_paths} },{'source_name':1,'source_pathways':1,'_id':0})
    temp=list(my_res)
    network_df=pd.DataFrame(temp)
    total_sig_prot=network_df.groupby(['source_pathways']).size()
    path_score=my_sig_prot.divide(total_sig_prot)
    #print("This is pathway score",path_score)
    #path_score
    #########score calculated above#########################
 
    my_data = db.signaling
    
    my_res=my_data.find({ 'source_name':{'$in':top_prots} },{'source_name':1,'source_pathways':1,'_id':0})
    cnt=0
    g = nx.Graph()
    colors = ['red', 'blue', 'green', 'yellow', 'pink']
    cy_nodes=[]
    cy_edges=[]
    #disease=[]
    nodes=set()
    for res in my_res:

        #disease=disease+res['Disease']
        source = res['source_name']
        target = res['source_pathways']
        
        
        cy_edges.append({  # Add the Edge Node
                'data': {'source': source, 'target': target}
                #,'classes': color
            })
        

        if source not in nodes:  # Add the source node
            nodes.add(source)
            cy_nodes.append({'data': {'id': source,'label':source,'size':10}, 'classes': 'search'})



        if target not in nodes:  # Add the target node
            nodes.add(target)
            #cy_nodes.append({'data': {'id': target,'label':target,'size':path_score[target]*1000 }, 'classes': 'not_search'})
            cy_nodes.append({'data': {'id': target,'label':target,'size':10 }, 'classes': 'not_search'})


    return cy_edges,cy_nodes,nodes



def generate_dummy_network():
    """Generate a default dummy network with placeholder nodes and edges."""
    dummy_nodes = [
        {'data': {'id': 'Dummy1', 'label': 'Placeholder Node 1', 'size': 20}, 'classes': 'dummy'},
        {'data': {'id': 'Dummy2', 'label': 'Placeholder Node 2', 'size': 20}, 'classes': 'dummy'}
    ]
    dummy_edges = [
        {'data': {'source': 'Dummy1', 'target': 'Dummy2'}}
    ]
    return dummy_edges, dummy_nodes, set(['Dummy1', 'Dummy2'])


def query_hallmark_net(top_prots):
    # Validate input
    if not top_prots:
        return generate_dummy_network()

    try:
        my_data = db.hallmarks

        # Query hallmark data
        my_res = my_data.find(
            {'Gene Symbol': {'$in': top_prots}},
            {'Gene Symbol': 1, 'Hallmark': 1, 'References': 1, '_id': 0}
        )
        temp = list(my_res)

        # If no results, return dummy data
        if not temp:
            return generate_dummy_network()

        # Build Cytoscape graph
        g = nx.Graph()
        cy_nodes = []
        cy_edges = []
        nodes = set()

        for res in temp:
            source = res.get('Gene Symbol')
            target = res.get('Hallmark')
            ref = res.get('References')

            # Skip if essential data is missing
            if not source or not target:
                continue

            # Add edge
            cy_edges.append({'data': {'source': source, 'target': target, 'Ref': ref}})

            # Add source node
            if source not in nodes:
                nodes.add(source)
                cy_nodes.append({'data': {'id': source, 'label': source}, 'classes': 'search'})

            # Add target node
            if target not in nodes:
                nodes.add(target)
                cy_nodes.append({'data': {'id': target, 'label': target}, 'classes': 'not_search'})

        return cy_edges, cy_nodes, nodes

    except Exception as e:
        # Log exception for debugging purposes (replace with proper logging in production)
        print(f"Error in query_hallmark_net: {e}")
        return generate_dummy_network()



#def get_top_tissue(top_prots):
def query_locations(top_prots):
    my_data = db.database1
    my_res=my_data.find({ 'name':{'$in':top_prots}},{'tissue_exp':1,'_id':0})
    all_tissues=[]
    for exp in my_res:
        #high_exp = list(filter(lambda x: '(High)' in x, exp['tissue_exp'])
        try:
            high_exp= [i for i in exp['tissue_exp'] if '(High)' in i]
            temp=[s.strip(' (High)') for s in high_exp]
            all_tissues.extend(temp)
        except:
            continue            
        #break
    all_tissues=list(set(all_tissues))
    my_exp_df=pd.DataFrame(index=top_prots,columns=all_tissues)
    my_exp_df=my_exp_df.fillna(0)
    my_res=my_data.find({ 'name':{'$in':top_prots}},{'tissue_exp':1,'name':1,'_id':0})
    for exp in my_res:
        #high_exp = list(filter(lambda x: '(High)' in x, exp['tissue_exp'])
        try:
            high_exp= [i for i in exp['tissue_exp'] if '(High)' in i]
            temp=[s.strip(' (High)') for s in high_exp]
            #all_tissues.extend(temp)
            my_exp_df.loc[exp['name'],temp]=1 
        except:
            continue
    #tiss_ = px.imshow(my_exp_df) 
    return my_exp_df


def get_path_cluster(top_prots):
    my_data = db.database1
    my_res=my_data.find({ 'name':{'$in':top_prots}},{'name':1,'pathways':1,'_id':0})

    all_tissues=[]
    for exp in my_res:
        #high_exp = list(filter(lambda x: '(High)' in x, exp['tissue_exp'])
        #path_way= [i for i in exp['pathways']]
        try:
            #path_ls=exp['pathways']
            #path_way= [i for i in path_ls]
            path_way= [i for i in exp['pathways']]
            #temp=[s.strip(' (High)') for s in high_exp]
            all_tissues.extend(path_way)
        except:
            continue            
        #break
    all_tissues=list(set(all_tissues))

    my_exp_df=pd.DataFrame(index=top_prots,columns=all_tissues)
    my_exp_df=my_exp_df.fillna(0)

    my_res=my_data.find({ 'name':{'$in':top_prots}},{'pathways':1,'name':1,'_id':0})
    for exp in my_res:
        #high_exp = list(filter(lambda x: '(High)' in x, exp['tissue_exp'])
        #high_exp= [i for i in exp['tissue_exp'] if '(High)' in i]
        try:
            temp=[i for i in exp['pathways']]
        #all_tissues.extend(temp)
            my_exp_df.loc[exp['name'],temp]=1 
        except:
            continue   
    #break
#all_tissues=list(set(all_tissues))
    df_to_return=my_exp_df.loc[:,my_exp_df.sum(axis=0)>3]
    return df_to_return



def get_dummy_fig(message="Load data to process"):

    try:
        fig = px.bar(x=["Step 1", "Step 2", "Step 3"], y=[0, 0, 0])
        fig.update_traces(marker_color="lightgray")
        fig.update_layout(
            xaxis=dict(showticklabels=False, title=None),
            yaxis=dict(showticklabels=False, title=None),
            annotations=[dict(
                text=message, xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=20, color="gray"), align="center"
            )]
        )
        return fig
    except Exception as e:
        print(f"Error: {e}")
        return go.Figure().add_annotation(
            text="Error creating plot", 
            xref="paper", yref="paper", showarrow=False,
            font=dict(size=16, color="red")
        )


def no_click_response():
    """Response when no analysis is triggered."""
    fig = get_dummy_fig()
    return [], fig, fig, fig, [], [], [], fig, [], []

def small_network_response():
    """Response for small or invalid networks."""
    fig = get_error_figure("Network size too small to analyze.")
    return [], fig, fig, fig, [], [], [], fig, [], ["Unable to analyze, network size too small."]


def get_error_figure(message):
    """Generate an error figure with a message."""
    fig = px.imshow(imread("assets/error_message.png"))
    fig.update_layout(title=message)
    return fig


def calculate_basic_stats(graph):

    # Node-level statistics
    degree = dict(graph.degree())
    in_degree = dict(graph.in_degree()) if graph.is_directed() else None
    out_degree = dict(graph.out_degree()) if graph.is_directed() else None
    betweenness = nx.betweenness_centrality(graph)
    clustering = nx.clustering(graph)
    closeness = nx.closeness_centrality(graph)

    # Combine node-level statistics into a DataFrame
    stats = {
        'Degree': degree,
        'InDegree': in_degree,
        'OutDegree': out_degree,
        'Betweenness': betweenness,
        'Clustering': clustering,
        'Closeness': closeness
    }
    node_stats_df = pd.DataFrame(stats).fillna(0)  # Fill NaN for undirected graphs
    node_stats_df.reset_index(level=0, inplace=True)
    node_stats_df = node_stats_df.rename(columns = {'index':'Name'})

    # Graph-level average statistics
    degree_values = np.array(list(degree.values()))
    non_zero_degrees = degree_values[degree_values != 0]
    avg_stats = {
        'num_nodes': graph.number_of_nodes(),
        'num_edges': graph.number_of_edges(),
        'average_degree': non_zero_degrees.mean() if non_zero_degrees.size > 0 else 0,
        'clustering_coefficient': nx.average_clustering(graph),
    }

    return node_stats_df, avg_stats


def get_network_elements(top_proteins):
    """Fetch network elements for hallmark and signaling networks."""
    hallmark_edges, hallmark_nodes, _ = query_hallmark_net(top_proteins)
    signa_edges, signa_nodes, _ = query_signa_net(top_proteins)
    return hallmark_edges + hallmark_nodes, signa_edges + signa_nodes



def get_top_stats(complex_para, top_proteins):

    # Fetch data from database
    top_complex_para = complex_para.loc[complex_para['Name'].isin(top_proteins),:]
    my_data = db.database1
    my_res = my_data.find(
        {'name': {'$in': top_proteins}},
        {'name': 1, 'location': 1, 'pclass': 1, 'pathways': 1, 'ligands': 1, 'pdb_id': 1, '_id': 0}
    )
    data = list(my_res)
    df = pd.DataFrame(data)

    # Build location DataFrame
    location_df = df.explode('location').dropna(subset=['location'])[['name', 'location']].rename(columns={'name': 'Name', 'location': 'Location'})

    # Build protein class DataFrame
    df_class = df[['name', 'pclass', 'location']].rename(
        columns={'name': 'Name', 'pclass': 'Protein Class', 'location': 'Location'}
    ).dropna()

    my_res=my_data.find({ 'name':{'$in':top_proteins}},{'name':1,'pdb_id':1,'ligands':1,'_id':0})
    temp=list(my_res)
    df=pd.DataFrame(temp)
    temp_struct_df=[]
    for id_,item in df.iterrows():
        prot_name=item['name']
        if item['ligands']!='NA':
            temp=pd.DataFrame(item['ligands'])
            lig_cnt=len(temp[temp == True].index)
            temp = temp.apply(lambda x : True if x['interaction_types'] == 'inhibitor' else False, axis = 1)
        else:
            lig_cnt=0
        if item['pdb_id']!='NA':
            pdb_cnt=len(item['pdb_id'])
        else:
            pdb_cnt=0
        temp_struct_df.append({'Name':prot_name,'Inhibitors':lig_cnt,'PDB complex(Count)':pdb_cnt}) 
    my_struct_df = pd.DataFrame(temp_struct_df)
    # Merge DataFrames
    top_stats = reduce(lambda left, right: pd.merge(left, right, on='Name', how='outer'), [top_complex_para, df_class, my_struct_df])

    top_stats[['Protein Class','Location']] = top_stats[['Protein Class','Location']].astype('string')
    # Generate location graph
    top_stats_graph_loc = px.bar(
        location_df,
        x="Location",
        y="Name",
        color="Name",
        title="Location of Top Spreaders",
        labels={"Location": "Protein Location", "Name": "Proteins"}
    )
    top_stats_graph_loc.update_yaxes(visible=False, showticklabels=False)

    return top_stats, top_stats_graph_loc


def get_location_graph(top_proteins):
    """Generate the location graph for top proteins."""
    # Assuming location data is fetched from MongoDB
    df = query_locations(top_proteins)
    #print('---data check----',df.shape)
    #print(df.head())
    #return px.bar(df, x="Location", y="Name", color="Name", title="Location of Top Proteins")

    my_fig=px.imshow(df,title="High tissue specificity of the spreaders",
                                    color_continuous_scale=[
                                    (0.00, "blue"),
                                    (0.5, "blue"),
                                    (0.5, "yellow"),
                                    (1.00, "yellow")]
                                                            )

    my_fig.update_layout(
                    xaxis=dict(tickmode='linear'),
                    yaxis=dict(tickmode='linear'),
                    coloraxis_colorbar=dict(
                    title="<b>Specificity</b>",
                    tickvals=[0,1],
                    ticktext=["<i>Not High</i>","<i>High</i>"]))
    return my_fig

def get_pathway_graph(top_proteins):
    """Generate pathway clustergram or dummy figure."""
    path_clus = get_path_cluster(top_proteins)
    if path_clus.empty:
        return get_dummy_fig()

    columns = list(path_clus.columns.values)
    rows = list(path_clus.index)
    path_graph=dashbio.Clustergram(
    data=path_clus.loc[rows].values,
    #title="Pathway clustergram of the triggers",
    column_labels=columns,
    row_labels=rows,

    color_map= [(0.00, "#EF553B"),   (0.5, "#EF553B"),
                (0.5, "#636EFA"),  (1.00, "#636EFA")],
    line_width=2,
    #hidden_labels='row',
    #cluster='column',
    height=800,
    width=1200
    )
    #path_graph.update_layout(showlegend=True)
    path_graph.data[-1].colorbar.tickvals=[-0.2,0.3]
    path_graph.data[-1].colorbar.ticktext=['Cluster 1','Cluster 2']

    return path_graph


def display_network_info(stats):
    """Create an HTML block for network statistics."""
    return html.Div(
        [
            html.H5("Network Statistics"),
            html.P(f"No. of Nodes: {stats['num_nodes']}"),
            html.P(f"No. of Edges: {stats['num_edges']}"),
            html.P(f"Average Degree: {stats['average_degree']:.3f}"),
            html.P(f"Clustering Coefficient: {stats['clustering_coefficient']:.3f}"),
        ],
        style={'padding': '10px', 'box-shadow': '0 3px 5px rgba(57, 63, 72, 0.3)'}
    )

def get_enrichment_data(my_prots, data_name):
    """
    Helper function to fetch enrichment data from the Enrichr API.
    """
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(my_prots)
    description = 'Example gene list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        return None, "connection_err"

    data = json.loads(response.text)
    user_list_id = data['userListId']

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = f'?userListId={user_list_id}&backgroundType={data_name}'
    response = requests.get(ENRICHR_URL + query_string)
    if not response.ok:
        return None, "connection_err"

    data = json.loads(response.text)
    return data[data_name], None


def create_bar_plot(data, title, x_label, y_label, top_n=10):
    """
    Helper function to create a bar plot from enrichment data.
    """
    rank, name, pval = zip(*[(d[0], d[1], d[2]) for d in data[:top_n]])
    fig = px.bar(
        x=-np.log10(pval),
        y=rank,
        #title=title,
        orientation='h',
        text=name,
        labels={ "x": x_label, "y": y_label }
    )
    fig.update_yaxes(autorange="reversed", showticklabels=False)
    fig.update_traces(textposition="inside")
    return fig

def register_callbacks(app):

    @app.callback(
    Output("file-download", "data"),
    Input("download-btn", "n_clicks"),
    #Input("file-dropdown", "value"),
    prevent_initial_call=True
    )
    def download_file(n_clicks):
        if n_clicks>0:
            print('directory')
            print(os.getcwd())
            file_path = 'assets/TCGA_GBM_small.csv'
            return dcc.send_file(file_path)
        return None

    ################################Data Callbacks begin####################################################
    @du.callback(
        #Output("callback-output", "children"),
        output = Output("upload_status","data"),
        id="upload-data",
    )
    def callback_on_completion(status: du.UploadStatus):
        #-----------------------------Review Code ---------------------------------------------
        if 'static_uuid' not in session:
            session['static_uuid'] = str(uuid4())

        session_upload_id = request.form.get("upload_id")

    # Use the static UUID for upload
        static_uuid = session['static_uuid']
         #-----------------------------Review Code ---------------------------------------------
        print(f'Static UUID: {static_uuid}')
        ctx = dash.callback_context
        ctx_msg = json.dumps({
            
            'triggered': ctx.triggered,
            'inputs': ctx.inputs,
            'state' : ctx.states
        }, indent=2)
        print('Triger 3: only after file upload',status)
        #print('Triger 3:', ctx_msg)
        print(f'Flask session is {session_upload_id}')
        
        # Display uploaded file paths
        if status.is_completed:
            csv_file_path = status.uploaded_files[0]
            print('debug 1', str(csv_file_path))
            #print('debug 2', type(status))
            return uploadstatus_to_dict(status)
            #return str(csv_file_path), uploadstatus_to_dict(status)


    @app.callback(
        Output('callback-counter','data'),
        Input('callback-counter','data')
        )
    def callback_init(init_cnt):

        init_cnt += 1
        print(f'======================== Callback init count {init_cnt}==========================')

        return init_cnt


    # Update page and selected tab
    @app.callback(
        Output('page-content','children'),
        Input('tabs','value'),# This will show the page content
        State('upload_status','data'), prevent_initial_call=True)
        #Input("callback-output", "children"), # This will create the page content
        #State('callback-output', 'children') # This is redundent
        
    def load_page(selected_tab,status_data):
        # Check condition that file uploader has data
        
        # if status_data is None:
        #     if 'data-tab' in selected_tab:
        #         return html.Div(['Please upload data'])#Data.layout()
        #     elif 'tab-1' in selected_tab:
        #         return html.Div(['Welcome to k2p'])

        #else:
        ctx = dash.callback_context
        ctx_msg = json.dumps({
            
            'triggered': ctx.triggered,
            'inputs': ctx.inputs,
            'state' : ctx.states
        }, indent=2)

        template_div = html.Div([

        html.Div([

            html.Div(["Please follow this flow chart for best results"],className="grid-item full-width")


            ],className="dashboard-grid"),



        html.Div(
            [ html.Div([html.Img(src=app.get_asset_url("landing_page.jpg"),  # Replace with your image in the `assets` folder
                                style={
                                        "width": "100%",  # Ensure image width fills parent
                                        "height": "100%",  # Ensure image height fills parent
                                        "object-fit": "cover"  # Maintains aspect ratio while filling
                                        })],className='grid-item full-width')
                      # Applies your provided CSS class
            ],className="dashboard-grid"
        )
        ])


        # print('New_check',status_data)
        # #file_path = status_data.uploaded_files[0]
        # file_path = status_data['uploaded_files'][0]

        # print('Here',file_path)
        #print('Trigger 1: page-content',ctx_msg)
        #if file_path:
          
            # user_upload = pd.read_csv(file_path)
            # uploaded_data_store = user_upload.to_dict('records')

        if 'data-tab' in selected_tab:
            return Data.layout()
        elif 'eda-tab' in selected_tab:
            return EDA.layout() #html.Div(['Content goes here']) #
        elif 'visualization-tab' in selected_tab:
            return Visualization.layout()
        elif 'network-tab' in selected_tab:
            return NetworkAnalysis.my_layout
        else:
            return template_div#'Welcome to K2P'
        #else:
            #raise PreventUpdate
            #return 'Please upload your data' 



    
    @app.callback(
        Output('my_sidebar', 'children'),
        Input('tabs', 'value'),  # Track active tab
        State('upload_status','data'),
        #State('callback-output', 'children'),
        prevent_initial_call=True
    )
    def display_sidebar(selected_tab,status_data):
        # Preventing update of side bar without value

        ctx = dash.callback_context
        ctx_msg = json.dumps({
            
            'triggered': ctx.triggered,
            'inputs': ctx.inputs,
            'state' : ctx.states
        }, indent=2)
        print('new_side bar', ctx_msg)

        # Condition when page lands
        if selected_tab == 'tab-1':
            my_div_tab1=get_tab1_sidebar()
            return my_div_tab1
            


        # Switching conditions when tab gets selected.
        if selected_tab in ['data-tab'] :
            upload_id = session.get('static_uuid')
            my_div_data = get_data_sidebar(upload_id)
            return my_div_data
        elif selected_tab == 'eda-tab':
            if status_data is None:
                raise PreventUpdate
            else:
                file_path = status_data['uploaded_files'][0]
                print('Inside sidebar trigger path:', file_path)

                # Determine file format and load accordingly
                if file_path.endswith('.csv'):
                    user_upload = pd.read_csv(file_path)
                elif file_path.endswith('.tsv'):
                    user_upload = pd.read_csv(file_path, sep="\t")
                elif file_path.endswith('.xls'):
                    user_upload = pd.read_excel(file_path)
            
                else:
                    raise ValueError("Unsupported file format. Please upload a .csv, .tsv, or .xls/.xlsx file.")

        # Convert uploaded data to a dictionary for use in Dash components
                uploaded_data_store = user_upload.to_dict('records')
                my_div_eda = get_eda_sidebar(uploaded_data_store)

            return my_div_eda
        elif selected_tab == 'visualization-tab':
            my_div_visualization = get_vis_sidebar()
            return my_div_visualization 
        elif selected_tab == 'network-tab':
            my_div_network = get_net_sidebar()
            return my_div_network
        else:
            return html.Div('No sidebar available for this tab.')


    ################################Network Callbacks begin####################################################

    @app.callback(Output('Textarea Input', 'value'),
                  Input('example_button', 'n_clicks'),prevent_initial_call=True)
    def update_layout(n_clicks):
        
        trigger_source = dash.callback_context.triggered
        
        
        #if curr_zoom:
        if trigger_source[0]['prop_id']=='example_button.n_clicks':
            return example_prots

    @app.callback(Output('cytoscape-update-layout', 'layout'),
              [Input('dropdown-update-layout', 'value')])
    def update_layout(layout):
        #print('change layout')
        return {
            'name': layout,
            'animate': True
        }

    @app.callback(
                  Output('intermediate-value','children'),
                  Output('unmatch-value','children'),
                  Input('save button', 'n_clicks'),
                  #State('upload-data', 'contents'),
                  #State('upload-data', 'filename'),
                  #State('upload-data', 'last_modified'),
                  State('Textarea Input', 'value'),
                  State('query_type','value')
              #[State('Filename Input', 'value')]
              )
    def storing_prots(n_clicks,inp_prot,qtype):
    #def storing_prots(n_clicks,list_of_contents, list_of_names, list_of_dates,inp_prot,qtype):


        ctx = dash.callback_context
        trigger_source=ctx.triggered
        ctx_msg = json.dumps({
            
            'triggered': ctx.triggered,
            'inputs': ctx.inputs,
            'state' : ctx.states
        }, indent=2)
        
        #print('The entry click with id type',ctx_msg)
        #print('states', ctx.states.keys())
        #trigger_source=ctx.triggered
        qtype_chk=ctx.states['query_type.value']

        if n_clicks == 0:
            print('network trigger callbacks', ctx_msg)
            raise PreventUpdate
        else:
            my_data = db.database1
            # if list_of_contents is not None:
            #         print('using file to fetch')
            #         div_children = parse_contents(list_of_contents, list_of_names) 
            if inp_prot is not None:
                print('using text area')
                div_children = [x.upper() for x in inp_prot.split("\n") if x]
            if qtype_chk=='gene_symbol':
                my_res=my_data.find({},{'name':1,'_id':0})
                db_prots=pd.DataFrame(list(my_res))
                nf=list(set(div_children)-set(db_prots['name'].to_list()))


                return div_children,nf
            elif qtype_chk=='uniprot_id':
                my_res=my_data.find({},{'id':1,'_id':0})
                db_prots=pd.DataFrame(list(my_res))
                nf=list(set(div_children)-set(db_prots['id'].to_list()))
                new_div=id_conversion(div_children)
                return new_div,nf

    @app.callback(
       Output("location", "value"),
       Output("function", "value"),
       Output("process", "value"),
       Output("pathway", "value"),
       Output("disease", "value"),
       Output("tissue", "value"),
       Output('apply_fil', 'n_clicks'),
       Input('reset','n_clicks')
        )
    def reset_filer(n_clicks):
        if n_clicks==0:
            raise PreventUpdate
        else:
            return None,None,None,None,None,None, 1
    
    @app.callback(
                Output(component_id='cytoscape-update-layout', component_property ='elements'),
                Input('intermediate-value','children'),
                Input('apply_fil', 'n_clicks'),
                State('location', 'value'),
                State('function', 'value'),
                State('process', 'value'),
                State('pathway', 'value'),
                State('disease', 'value'),
                State('tissue', 'value'),
                State('memory_net', 'data')
          #[State('Filename Input', 'value')]
          )

    #def storing_prots(n_clicks,sel_loc,sel_fun,sel_proc,sel_path,sel_dis,sel_tiss,list_of_contents, list_of_names, list_of_dates,session_net,inp_prot,curr_net):
    def storing_prots(div_children,filters,sel_loc,sel_fun,sel_proc,sel_path,sel_dis,sel_tiss,sess_net):

        #part1
        ctx = dash.callback_context
        #print('states', ctx.states.keys())
        trigger_source=ctx.triggered
        ctx_msg = json.dumps({
            
            'triggered': ctx.triggered,
            'inputs': ctx.inputs
        }, indent=2)
        #print('Trigger of cb for network generation',ctx_msg)

        if div_children is None and filters == 0:
            raise PreventUpdate
        else:
            

            query=div_children
            #query=prots
            #print('user input',prots)
            #data = db.database1#check that after making entry this line shud be re-runned or not ?
            if trigger_source[0]['prop_id']=='intermediate-value.children':
                cy_edges,cy_nodes,nodes=query_db_net(query)
                #print('Today we found them',nodes)
                print(f'Network info: Number of nodes {len(cy_nodes)}and Number of edges {len(cy_edges)}')
                #options,options1,options2,options3,options4,options5=fetch_dropdown(nodes)
                return  cy_edges+cy_nodes


            elif trigger_source[0]['prop_id']=='apply_fil.n_clicks':
                ctx = dash.callback_context
                call_state=ctx.states


                loc_fil=call_state["location.value"]
                fun_fil=call_state["function.value"]
                proc_fil=call_state["process.value"]
                path_fil=call_state["pathway.value"]
                dis_fil=call_state["disease.value"]
                exp_fil=call_state["tissue.value"]

                filter_ls=multi_filter_query(loc_fil,fun_fil,proc_fil,path_fil,dis_fil,exp_fil)
                print(f'Network Query {len(filter_ls)}')
                new_el=update_nodes(sess_net,filter_ls)
                return new_el


    @app.callback(
       Output("location", "options"),
       Output("function", "options"),
       Output("process", "options"),
       Output("pathway", "options"),
       Output("disease", "options"),
       Output("tissue", "options"),
       Output('memory_net', 'data'),
       Input('intermediate-value','children')    
        )
    def drop_down_feed(div_children):

        ctx = dash.callback_context
        #print('states', ctx.states.keys())
        trigger_source=ctx.triggered
        ctx_msg = json.dumps({
            
            'triggered': ctx.triggered,
            'inputs': ctx.inputs,
            'state' : ctx.states
        }, indent=2)
        #print('Trigger of cb for new drop down',ctx_msg)


        if div_children is None:
            print('No trigger for down')
            raise PreventUpdate
        else:
            query=div_children
            if trigger_source[0]['prop_id']=='intermediate-value.children':
                cy_edges,cy_nodes,nodes=query_db_net(query)
                #print('Today we found them',nodes)
                
                options,options1,options2,options3,options4,options5=fetch_dropdown(nodes)
                print(f'Number of nodes {len(cy_nodes)}and Number of edges {len(cy_edges)} and options{len(options)}')

                return  options, options1, options2, options3,options4,options5,cy_edges+cy_nodes

    
    @app.callback(Output('ligand', 'children'), Input('cytoscape-update-layout', 'tapNodeData'))
    def display_ligand_data(data):
        if not data:
            raise PreventUpdate
        
        query_result = fetch_data_from_db('database1', {'name': data['id']}, {'_id': 0, 'ligands': 1})
        ligands = query_result[0].get('ligands', 'NA')

        if ligands != 'NA':
            columns = ['ligand_name', 'drug_concept_id', 'interaction_types', 'PMIDs']
            df = pd.DataFrame(ligands)[columns]
            return generate_data_table(df.to_dict('records'), columns)
        else:
            return html.Div(["No data found"], style={'height': '100%', 'width': '100%'})

    # Callback for 'docking' data
    @app.callback(Output('docking', 'children'), Input('cytoscape-update-layout', 'tapNodeData'))
    def display_docking_data(data):
        if not data:
            raise PreventUpdate

        query_result = fetch_data_from_db('database1', {'name': data['id']}, {'_id': 0, 'structure_info': 1})
        structure_info = query_result[0].get('structure_info', [])

        if structure_info:
            df = pd.DataFrame(structure_info).drop(columns=['Gene_ID'], errors='ignore')
            return generate_data_table(df.to_dict('records'), df.columns)
        else:
            return html.Div(["No data found"], style={'height': '100%', 'width': '100%'})

    # Callback for 'disease_mutation' data
    @app.callback(Output('disease_mutation', 'children'), Input('cytoscape-update-layout', 'tapNodeData'))
    def display_disease_mutation_data(data):
        if not data:
            raise PreventUpdate

        query_result = fetch_data_from_db('database1', {'name': data['id']}, {'_id': 0, 'disease_mut': 1})
        disease_mut = query_result[0].get('disease_mut', 'NA')

        if disease_mut != 'NA':
            df = pd.DataFrame(disease_mut).drop(columns=['gene', 'uniprot'], errors='ignore')
            return generate_data_table(df.to_dict('records'), df.columns)
        else:
            return html.Div(["No data found"], style={'height': '100%', 'width': '100%'})

    # Callback for 'func_info' data
    @app.callback(Output('func_info', 'children'), Input('cytoscape-update-layout', 'tapNodeData'))
    def display_func_info(data):
        if not data:
            raise PreventUpdate

        gene, protein, func_txt = get_func_info(data['id'])
        label = f'Uniprot: {protein}, Gene: {gene}'
        return html.Div(
            [
                html.H5(label),
                html.P(func_txt),
            ],
            style={'padding': '10px', 'margin': '10px', 'border-color': '#dce8fc'}
        )

    @app.callback(Output('pocket_info', 'children'), Input('cytoscape-update-layout', 'tapNodeData'))
    def display_residue_data(data):
        if not data:
            raise PreventUpdate

        db = client.ppi_test
        info_db = db.database1
        #query_result = fetch_data_from_db('database1', {'name': data['id']}, {'_id': 0, 'disease_mut': 1})
        #disease_mut = query_result[0].get('disease_mut', 'NA')
        #pocket_df = pd.DataFrame(list(info_db.find({'name':data['id']},{'id':1,'pocket_data':1,'_id':0}))[0]['pocket_data'])
        entry = list(info_db.find({'name': data['id']}, {'id': 1, 'pocket_data': 1, '_id': 0}))

        # Check if the entry exists and has the 'pocket_data' attribute
        if entry and 'pocket_data' in entry[0]:
            pocket_df = pd.DataFrame(entry[0]['pocket_data'])
        else:
            pocket_df = pd.DataFrame()

        print('checking residue',pocket_df.shape)
        if not pocket_df.empty:
            df = pocket_df.drop(columns=[' pocket'],axis=1, errors='ignore')
            return generate_data_table(df.to_dict('records'), df.columns)
        else:
            return html.Div(["No data found"], style={'height': '100%', 'width': '100%'})

    # Callback for 'click_info' (edge data)
    @app.callback(Output('click_info', 'children'), Input('cytoscape-update-layout', 'tapEdgeData'))
    def display_click_info(data):
        if not data:
            raise PreventUpdate

        query_result = fetch_data_from_db(
            'edge_list_ref',
            {'$and': [{'ENTITYA': data['source']}, {'ENTITYB': data['target']}]},
            {'ENTITYA': 1, 'IDA': 1, 'ENTITYB': 1, 'IDB': 1, 'EFFECT': 1, 'MECHANISM': 1, 'PMID': 1, '_id': 0}
        )
        df = pd.DataFrame(query_result)

        if not df.empty:
            return generate_data_table(df.to_dict('records'), df.columns)
        else:
            return html.Div(["No data found"], style={'height': '100%', 'width': '100%'})

    # Callback for tabs content
    @app.callback(
        Output('tabs-content-inline', 'children'),
        Input('tabs-styled-with-inline', 'value'),
        State('func_info', 'children'),
        State('click_info', 'children'),
        State('docking', 'children'),
        State('ligand', 'children'),
        State('disease_mutation', 'children'),
        State('pocket_info','children')
    )
    def render_content(tab, func_info, click_info, docking, ligand, disease_mutation,pocket_info):
        if tab == 'tab-0':
            return func_info
        elif tab == 'tab-1':
            return docking
        elif tab == 'tab-2':
            return ligand
        elif tab == 'tab-3':
            return disease_mutation
        elif tab == 'tab-4':
            return pocket_info
        elif tab == 'tab-5':
            return click_info

##############################################################################
    @app.callback(
    [
        Output('central_para', 'data'),
        Output('network_para_graph', 'figure'),
        Output('location_graph', 'figure'),
        Output('sub_network_para_graph', 'figure'),
        Output('hallmark_network', 'elements'),
        Output('signa_network', 'elements'),
        Output('top_stats', 'data'),
        Output('top_pathway', 'figure'),
        Output('net_info', 'children'),
        Output('analyse_status', 'children')
    ],
    [Input('cp_calculate', 'n_clicks')],
    [State('cytoscape-update-layout', 'elements')]
    )
    def calculate_network_para(n_clicks, curr_elements):

        ctx = dash.callback_context
        ctx_msg = json.dumps({
            
            'triggered': ctx.triggered,
            'inputs': ctx.inputs,
            'state' : ctx.states
        }, indent=2)

        
        """Analyze the network parameters and update the dashboard components."""
        if n_clicks <= 0:
            #print('Analysis CBs',ctx_msg)
            return no_click_response()

        print('Starting analysis...')
        graph = cyto_2_graph(curr_elements)

        if graph.number_of_edges() == 0:
            return small_network_response()

        # Topological analysis

        # Basic network statistics
        network_stats, average_stats = calculate_basic_stats(graph)

        
        # Top proteins and related queries
        top_proteins = nx.voterank(graph, 15)
        hallmark_elements, signa_elements = get_network_elements(top_proteins)

        # Visualization and statistics
        graph_deg_vs_bet = px.scatter(network_stats, x="Degree", y="Betweenness",title="Degree vs Betweenness", hover_data=['Name'],color="Degree")
        top_stats, location_graph = get_top_stats(network_stats, top_proteins)
        print('here is final table')
        print(top_stats.dtypes)
        tissue_graph = get_location_graph(top_proteins)
        pathway_graph = get_pathway_graph(top_proteins)
        net_info = display_network_info(average_stats)

        print('Analysis complete.')
        return (
            network_stats.to_dict('records'),  # Central Parameters Data
            graph_deg_vs_bet,  # Network Parameters Graph
            tissue_graph,  # Location Graph
            location_graph,  # Pathway Graph
            hallmark_elements,  # Hallmark Network
            signa_elements,  # Signaling Network
            top_stats.to_dict('records'),  # Top Stats Data
            pathway_graph,  # Pathway Figure
            net_info,  # Network Info
            ['Scroll down to global properties panel for results']  # Analysis Status
        )

    @app.callback(
    Output('process_graph', 'figure'),
    Input('cp_calculate', 'n_clicks'),
    State('cytoscape-update-layout', 'elements')
    )
    def make_process_graph(clicked, sess_net):
        if clicked > 0:
            my_prots = get_nodes(sess_net)
            data_name = 'GO_Biological_Process_2018'
            data, err = get_enrichment_data(my_prots, data_name)

            if err == "connection_err":
                #disp_img = imread("../assets/error_message.png")
                return get_error_figure('Connection Error')

            return create_bar_plot(data, "Enriched processes of interactome", "-log(P value)", "Processes")
        return get_dummy_fig()


    @app.callback(
        Output('pclass_graph', 'figure'),
        Input('cp_calculate', 'n_clicks'),
        State('cytoscape-update-layout', 'elements')
    )
    def make_protein_class_graph(n_clicks, sess_net):
        if n_clicks > 0:
            total_prot = get_nodes(sess_net)
            #db = client.ppi
            my_data = db.database1
            my_res = my_data.find({'name': {'$in': total_prot}}, {'pclass': 1})

            p_class = []
            for res in my_res:
                try:
                    p_class.extend(res['pclass'])
                except Exception:
                    continue

            frequency = collections.Counter(p_class)
            fig = px.bar(
                x=list(frequency.keys()),
                y=list(frequency.values()),
                title="Protein class abundance of interactome",
                labels={"x": "Protein classes", "y": "Abundance"}
            )
            return fig
        return get_dummy_fig()


    @app.callback(
        Output('pathway_graph', 'figure'),
        Input('cp_calculate', 'n_clicks'),
        State('cytoscape-update-layout', 'elements')
    )
    def make_pathway_graph(clicked, sess_net):
        if clicked > 0:
            my_prots = get_nodes(sess_net)
            data_name = 'KEGG_2019_Human'
            data, err = get_enrichment_data(my_prots, data_name)

            if err == "connection_err":
                #disp_img = imread("../assets/connection_err.png")
                return get_error_figure('Connection Error')

            return create_bar_plot(data, "Enriched pathways of interactome", "-log(P value)", "Pathways")
        return get_dummy_fig()


    @app.callback(
        Output('disease_graph', 'figure'),
        Input('cp_calculate', 'n_clicks'),
        State('cytoscape-update-layout', 'elements')
    )
    def make_disease_graph(n_clicks, sess_net):
        if n_clicks > 0:
            my_prots = get_nodes(sess_net)
            data_name = 'ClinVar_2019'
            data, err = get_enrichment_data(my_prots, data_name)

            if err == "connection_err":
                #disp_img = imread("G:/PPI/website/website(konnect2prot) v4.5.10/assets/connection_err.png")
                return get_error_figure('Connection Error')

            if not data:
                disp_img = imread("G:/PPI/website/website(konnect2prot) v4.5.10/assets/err_msg_size.png")
                return px.imshow(disp_img).update_xaxes(showticklabels=False).update_yaxes(showticklabels=False)

            return create_bar_plot(data, "Enriched diseases of interactome", "-log(P value)", "Diseases")
        return get_dummy_fig()

    # -----------------------------------
    # Helper Functions
    # -----------------------------------

    
