import pandas as pd
import numpy as np
from dash import dcc, html, dash_table, callback
from dash.dependencies import Input, Output, State
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from scipy.stats import ttest_ind, shapiro
import flask


from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind, shapiro
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import numpy as np
import pandas as pd

def preprocess_and_find_DEGs(my_df, meta_df, norm_method, scaling_method, alpha=0.05, normality_test=False):

    #print('meta_table')
    #print(meta_df)
    sample_ids = list(my_df.columns[1:])  # Exclude the index column
    group_labels = meta_df.set_index('ID').loc[sample_ids, 'Label'].values

    # Extract data
    data = my_df.iloc[:, 1:].values  # Drop the first column (gene names)

    # Normalization
    if norm_method == 'log2':
        data = np.log2(data + 1)
    elif norm_method == 'log10':
        data = np.log10(data + 1)

    # Scaling
    if scaling_method == 'standard':
        scaler = StandardScaler()
        data = scaler.fit_transform(data.T).T
    elif scaling_method == 'minmax':
        scaler = MinMaxScaler()
        data = scaler.fit_transform(data.T).T

    # Differential Expression Analysis
    unique_groups = np.unique(group_labels)
    if len(unique_groups) != 2:
        raise ValueError("Only two-group comparisons are supported in this function.")

    group1 = data[:, np.array(group_labels) == unique_groups[0]]
    group2 = data[:, np.array(group_labels) == unique_groups[1]]

    # Parametric test
    t_stats, p_values = ttest_ind(group1, group2, axis=1)

    # Normality test (optional)
    if normality_test:
        p_values_normality = [shapiro(group)[1] for group in [group1, group2]]
        print(f"Shapiro-Wilk normality test p-values: {p_values_normality}")

    # Fold change and log2 fold change
    mean_group1 = np.mean(group1, axis=1)
    mean_group2 = np.mean(group2, axis=1)
    fold_change = mean_group2 / mean_group1
    log2_fold_change = np.log2(fold_change)

    # FDR correction
    _, fdr_p_values, _, _ = multipletests(p_values, method='fdr_bh')

    # Extract gene names from the index column
    gene_names = my_df.iloc[:, 0].tolist()

    # Construct the DEGs DataFrame
    degs = pd.DataFrame({
        'Gene_id': gene_names,
        'p-value': p_values,
        'FDR': fdr_p_values,
        'Fold Change': fold_change,
        'log2_Fold_Change': log2_fold_change,
        'Significant': p_values < alpha
    })

    return degs.sort_values('p-value').reset_index(drop=True)

# def preprocess_and_find_DEGs(my_df, group_labels, norm_method, scaling_method, alpha=0.05, normality_test=False):

#     data = my_df.drop(my_df.columns[0], axis=1).values

#     #print("check gene names",data.index[:10])
#     # 1. Normalization
#     if norm_method == 'log2':
#         data = np.log2(data + 1)
#     elif norm_method == 'log10':
#         data = np.log10(data + 1)

#     # 2. Scaling
#     if scaling_method == 'standard':
#         scaler = StandardScaler()
#         data = scaler.fit_transform(data.T).T
#     elif scaling_method == 'minmax':
#         scaler = MinMaxScaler()
#         data = scaler.fit_transform(data.T).T

#     # 3. Differential Expression Analysis
#     unique_groups = np.unique(group_labels)
#     if len(unique_groups) != 2:
#         raise ValueError("Only two-group comparisons are supported in this function.")

#     group1 = data[:, np.array(group_labels) == unique_groups[0]]
#     group2 = data[:, np.array(group_labels) == unique_groups[1]]
    
#     # Parametric test
#     t_stats, p_values = ttest_ind(group1, group2, axis=1)

#     # Normality test (optional)
#     if normality_test:
#         p_values_normality = [shapiro(group)[1] for group in [group1, group2]]
#         print(f"Shapiro-Wilk normality test p-values: {p_values_normality}")
    
#     mean_group1 = np.mean(group1, axis=1)
#     mean_group2 = np.mean(group2, axis=1)
#     fold_change = mean_group2 / mean_group1
#     log2_fold_change = np.log2(fold_change)

#     #gene_names = [f'Gene_{i+1}' for i in range(data.shape[0])]  # Placeholder gene names
#     gene_names = list(my_df.iloc[:,0])

#     degs = pd.DataFrame({
#         'Gene_id': gene_names,
#         'p-value': p_values,
#         'Fold Change':fold_change,
#         'log2_Fold_Change': log2_fold_change,
#         'Significant': p_values < alpha
#     })
#     degs=degs.reset_index(drop=True)
#     return degs.sort_values('p-value')