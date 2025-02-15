import plotly.express as px
import pandas as pd
import requests
import json



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


    enrich_table = json.loads(response.text)[data_name]
    '''
    enrich_table is list of list, in each entry there will be items in below order:
    
    Rank, Term name, P-value, Odds ratio, Combined score, Overlapping genes, Adjusted p-value
    '''
    
    
    return enrich_table, None


# Function to plot enrichment results
def plot_enrichment_bubble(enrich_table, Y_label = "Pathway"):

    # Convert enrich_table into a DataFrame
    columns = ['Rank', 'Pathway', 'P-value', 'Odds ratio', 'Combined score', 'Overlapping genes', 'Adjusted P-value','dummy1','dummy2']
    enrich_df = pd.DataFrame(enrich_table, columns=columns)
    #print("here is the enrich_table",enrich_df.shape)

    # Process Overlapping genes column (it's a string, we count the genes)
    enrich_df['Num overlapping genes'] = enrich_df['Overlapping genes'].apply(lambda x: len(x))

    # Create bubble plot
    fig = px.scatter(
        enrich_df,
        x='Odds ratio',
        y='Pathway',
        size='Num overlapping genes',
        color='P-value',
        color_continuous_scale='Viridis_r',  # Colors for p-values (inverted for low p-value emphasis)
        hover_data={
            'P-value': ':.3e',  # Scientific notation for p-values
            'Adjusted P-value': ':.3e',
            'Num overlapping genes': True,
            'Odds ratio': ':.2f',
            'Pathway': True,
        },
        #title='Pathway Enrichment Analysis',
        labels={
            'Odds ratio': 'Odds Ratio',
            'Pathway': 'Enriched Pathway',
            'P-value': 'P-value',
            'Num overlapping genes': 'Overlapping Genes',
        },
    )
    fig.update_traces(marker=dict(line=dict(width=0.5, color='DarkSlateGrey')))
    fig.update_layout(
        xaxis=dict(title='Odds Ratio'),
        yaxis=dict(title=Y_label),
        template='plotly_white',
    )
    return fig