from dash import dcc, html, callback, dash_table, ctx
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from io import BytesIO
import base64
import scipy.stats as stats

# Set the backend to Agg to avoid GUI errors in the Dash environment
plt.switch_backend('Agg')


def layout():
    my_data = html.Div([
        html.H3('Uploaded Data Preview', className="panel-title"),
        
        # Data Table Section
        html.Div(
            [
                html.Div([
                    html.H4('Data Table', className="panel-title"),
                    dcc.Loading(
                        dash_table.DataTable(
                            id="data-table",
                            style_table={'overflowX': 'scroll'},
                            style_cell={"textAlign": "center", "padding": "10px"},
                            style_header={
                                "backgroundColor": "#f8f9fa",
                                "fontWeight": "bold",
                                "border": "1px solid #ddd",
                            },
                        ),
                    )
                ], className="grid-item full-width"),
            ],
            className="dashboard-grid"
        ),

        # Outlier Detection Section
        

        # Plots Section
        html.Div(
            [
                # Box Plot
                html.Div(
                    [
                        html.Div("Box Plot of Gene Expression Across Samples", className="panel-title"),
                        dcc.Loading(html.Img(id="box-plot", style={"width": "100%"})),
                    ],
                    className="grid-item",
                ),

                # Mean-Variance Trend
                html.Div(
                    [
                        html.Div("Mean-Variance Trend (log2(σ) vs Average Log-Expression)", className="panel-title"),
                        dcc.Loading(html.Img(id="mean-variance-plot", style={"width": "100%"})),
                    ],
                    className="grid-item",
                ),
            ],
            className="dashboard-grid",
        ),

        html.Div(
            [
                # Density Plot
                html.Div(
                    [
                        html.Div("Density Plot of Gene Expression for Each Sample", className="panel-title"),
                        dcc.Loading(html.Img(id="density-plot", style={"width": "100%"})),
                    ],
                    className="grid-item",
                ),

                # QQ Plot
                html.Div(
                    [
                        html.Div("QQ Plot for Each Sample", className="panel-title"),
                        dcc.Loading(html.Img(id="qq-plot", style={"width": "100%"})),
                    ],
                    className="grid-item",
                ),
            ],
            className="dashboard-grid",
        )
    ])
    return my_data


def fig_to_base64(fig):
    """Convert a matplotlib figure to a base64 image for rendering in Dash."""
    buffer = BytesIO()
    fig.savefig(buffer, format='png', bbox_inches='tight')
    plt.close(fig)
    buffer.seek(0)
    img_str = base64.b64encode(buffer.read()).decode('utf-8')
    return 'data:image/png;base64,' + img_str


def load_data(file_path):
    """Load data from a CSV or TSV file."""
    if file_path.endswith('.csv'):
        return pd.read_csv(file_path)
    elif file_path.endswith('.tsv'):
        return pd.read_csv(file_path , sep="\t")
    elif file_path.endswith(".xls"):
        return pd.read_excel(file_path)
    else:
        raise ValueError("Unsupported file format. Only CSV and TSV files are allowed.")


@callback(
    Output("data-table", "data" ,allow_duplicate=True),
    Output("data-table", "columns" ,allow_duplicate=True),
    Output("upload_status", "data" , allow_duplicate=True),
    Input("apply-outlier", "n_clicks"),
    State("upload_status", "data"),
    State("floor-percentile", "value"),
    State("cap-percentile", "value"),
    prevent_initial_call=True
)
def apply_outlier_detection_and_update_table(n_clicks, status_data, lower_percentile, upper_percentile):
    """Apply flooring and capping for outlier removal and update the data table."""
    if not status_data:
        raise PreventUpdate

    file_path = status_data['uploaded_files'][0]
    df = load_data(file_path)

    # Apply flooring and capping
    if n_clicks>0:
        sample_columns = df.columns[1:]  # Exclude the ID column
        for col in sample_columns:
            lower_bound = np.percentile(df[col], lower_percentile)
            upper_bound = np.percentile(df[col], upper_percentile)
            df[col] = np.clip(df[col], lower_bound, upper_bound)

        # Save updated dataset
        updated_file_path = file_path.replace(".csv", "_updated.csv").replace(".tsv", "_updated.tsv")
        df.to_csv(updated_file_path, sep='\t' if file_path.endswith('.tsv') else ',', index=False)
        status_data['uploaded_files'][0] = updated_file_path

    # Update table
    columns = [{"name": col, "id": col} for col in df.columns]
    data = df.head(5).to_dict("records")

    return data, columns, status_data


@callback(
    Output("data-table", "data" , allow_duplicate=True),
    Output("data-table", "columns" , allow_duplicate=True),
    Input("upload_status", "data"),
    prevent_initial_call=True
)
def display_uploaded_data(status_data):
    """Display the uploaded data when the file is first uploaded."""
    if not status_data:
        raise PreventUpdate

    file_path = status_data['uploaded_files'][0]
    df = load_data(file_path)
    columns = [{"name": col, "id": col} for col in df.columns]
    data = df.head(5).to_dict("records")

    return data, columns


@callback(
    Output("box-plot", "src"),
    Input("upload_status", "data"),
)
def generate_box_plot(status_data):
    """Generate box plot as a static image."""
    if not status_data:
        raise PreventUpdate

    file_path = status_data['uploaded_files'][0]
    df = load_data(file_path)
    melted_df = df.melt(id_vars=df.columns[0], var_name="Sample", value_name="Expression")

    fig, ax = plt.subplots(figsize=(12, 10))
    sns.boxplot(data=melted_df, x="Expression", y="Sample", ax=ax, color="#4CAF50", linewidth=1.2)
    ax.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.set_xlabel("Gene Expression", fontsize=24)
    ax.set_ylabel("Sample", fontsize=24)
    sns.despine(offset=10, trim=True)

    return fig_to_base64(fig)


@callback(
    Output("mean-variance-plot", "src"),
    Input("upload_status", "data"),
)
def generate_mean_variance_plot(status_data):
    """Generate mean-variance trend as a static image."""
    if not status_data:
        raise PreventUpdate

    file_path = status_data['uploaded_files'][0]
    df = load_data(file_path)
    sample_columns = df.columns[1:]

    log_mean = np.log2(df[sample_columns].mean(axis=1) + 1e-5)
    log_std = np.log2(df[sample_columns].std(axis=1) + 1e-5)

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.scatterplot(x=log_mean, y=log_std, ax=ax, color="blue", edgecolor="black", s=60, alpha=0.6)
    ax.set_xlabel("Average Log-Expression (log2)", fontsize=20)
    ax.set_ylabel("Log2(Standard Deviation)", fontsize=20)

    return fig_to_base64(fig)


@callback(
    Output("density-plot", "src"),
    Input("upload_status", "data"),
)
def generate_density_plot(status_data):
    """Generate density plot as a static image."""
    if not status_data:
        raise PreventUpdate

    file_path = status_data['uploaded_files'][0]
    df = load_data(file_path)
    sample_columns = df.columns[1:]

    fig, ax = plt.subplots(figsize=(10, 8))
    for sample in sample_columns:
        sns.kdeplot(df[sample], ax=ax, label=sample, linewidth=1.5)
    ax.set_xlabel("Gene Expression", fontsize=20)
    ax.set_ylabel("Density", fontsize=20)
    return fig_to_base64(fig)


@callback(
    Output("qq-plot", "src"),
    Input("upload_status", "data"),
)
def generate_qq_plot(status_data):
    """Generate QQ plot as a static image."""
    if not status_data:
        raise PreventUpdate

    file_path = status_data['uploaded_files'][0]
    df = load_data(file_path)
    sample_columns = df.columns[1:]

    fig, ax = plt.subplots(figsize=(10, 8))
    for sample in sample_columns:
        stats.probplot(df[sample], dist="norm", plot=ax)
    ax.set_xlabel("Theoretical Quantiles", fontsize=20)
    ax.set_ylabel("Sample Quantiles", fontsize=20)

    return fig_to_base64(fig)
