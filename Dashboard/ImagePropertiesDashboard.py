from dash import Dash, dcc, html, Input, Output, State, dash_table
import pandas as pd
import plotly.express as px
import base64
import io
import tempfile
import nibabel as nib
import numpy as np
import os
import scipy.stats as stats
import sys

# Function to load NIfTI files and compute numerical properties

# Function to load NIfTI files and compute numerical properties
def load_nifti_data(contents, filename):
    try:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        file_extension = '.nii.gz' if filename.endswith('.gz') else '.nii'
        with tempfile.NamedTemporaryFile(delete=False, suffix=file_extension) as tmp:
            tmp.write(decoded)
            tmp.flush()
            tmp.close()
            img = nib.load(tmp.name)
        img_data = img.get_fdata()
        img_data_nonzero = img_data[img_data > 0]  # consider only non-zero voxels
        
        mean_intensity = np.mean(img_data)
        std_intensity = np.std(img_data)
        max_intensity = np.max(img_data)
        min_intensity = np.min(img_data)
        entropy = stats.entropy(np.histogram(img_data_nonzero, bins=100)[0])
        skewness = stats.skew(img_data_nonzero)
        kurtosis = stats.kurtosis(img_data_nonzero)
        contrast_ratio = np.max(img_data_nonzero) / np.mean(img_data_nonzero)
        percentile_10 = np.percentile(img_data_nonzero, 10)
        percentile_90 = np.percentile(img_data_nonzero, 90)

        os.unlink(tmp.name)
        return {
            'File Name': filename,
            'Mean Intensity': mean_intensity,
            'Standard Deviation': std_intensity,
            'Max Intensity': max_intensity,
            'Min Intensity': min_intensity,
            'Entropy': entropy,
            'Skewness': skewness,
            'Kurtosis': kurtosis,
            'Contrast Ratio': contrast_ratio,
            '10th Percentile Intensity': percentile_10,
            '90th Percentile Intensity': percentile_90
        }
    except Exception as e:
        print(f"Error processing file {filename}: {e}")
        return None
app = Dash(__name__)

app.layout = html.Div([
    html.H1("View NIfTI Image Properties", style={'text-align': 'center', 'margin-bottom': '20px'}),
    dcc.Upload(
        id='upload-data',
        children=html.Button('Upload Files', style={'width': '100%', 'padding': '10px', 'font-size': '20px'}),
        multiple=True
    ),
    html.Div(id='output-data-upload', style={'margin-top': '20px', 'margin-bottom': '20px'}),
    dcc.Graph(id='line-plot', style={'margin-bottom': '20px'}),
    dash_table.DataTable(
        id='data-table',
        page_size=10,
        style_table={'height': '300px', 'overflowY': 'auto', 'clear': 'both'},
        sort_action="native",  # Enable sorting on all columns
        style_cell_conditional=[
            {'if': {'column_id': 'pagination'},
             'textAlign': 'left'}
        ],
        css=[{"selector": ".previous-page, .next-page", "rule": "float: left; margin-right: 10px;"}]
    )
], style={'padding': '20px', 'font-family': 'Arial, sans-serif'})

@app.callback(
    [Output('output-data-upload', 'children'),
     Output('line-plot', 'figure'),
     Output('data-table', 'data'),
     Output('data-table', 'columns')],
    [Input('upload-data', 'contents')],
    [State('upload-data', 'filename')]
)
def update_output(list_of_contents, list_of_names):
    if list_of_contents is not None:
        data = []
        for c, n in zip(list_of_contents, list_of_names):
            data.append(load_nifti_data(c, n))
        df = pd.DataFrame(data)
        if df.empty:
            return 'No NIfTI files found.', {}, [], []
        df_melted = df.melt(id_vars=['File Name'], var_name='Property', value_name='Value')
        line_fig = px.line(df_melted, x='Property', y='Value', color='File Name', title='Properties of NIfTI Files')
        line_fig.update_layout(showlegend=False)
        columns = [{"name": i, "id": i} for i in df.columns]
        data = df.to_dict('records')
        return 'Data loaded successfully', line_fig, data, columns
    return 'Upload NIfTI files to view properties.', {}, [], []

if __name__ == "__main__":
    # Read the port from the command line
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8050
    # Start the app with the specified port
    app.run_server(host="0.0.0.0", port=port, debug=True)

