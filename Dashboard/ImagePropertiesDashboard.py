from dash import Dash, dcc, html, Input, Output, State, dash_table
import dash
import pandas as pd
import plotly.express as px
import base64
import io
import tempfile
import nibabel as nib
import numpy as np
import os
import scipy.stats as stats

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

        # First-order statistics
        mean_intensity = np.mean(img_data)
        median_intensity = np.median(img_data)
        std_intensity = np.std(img_data)
        variance_intensity = np.var(img_data)
        max_intensity = np.max(img_data)
        min_intensity = np.min(img_data)
        range_intensity = max_intensity - min_intensity

        # Histogram for non-zero voxels
        hist_nonzero, bin_edges_nonzero = np.histogram(img_data_nonzero, bins=100)
        entropy = stats.entropy(hist_nonzero, base=2)
        skewness = stats.skew(img_data_nonzero)
        kurtosis = stats.kurtosis(img_data_nonzero)
        contrast_ratio = np.max(img_data_nonzero) / np.mean(img_data_nonzero)

        percentile_10 = np.percentile(img_data_nonzero, 10)
        percentile_90 = np.percentile(img_data_nonzero, 90)

        # Mean Absolute Deviation
        mean_abs_dev = np.mean(np.abs(img_data_nonzero - np.mean(img_data_nonzero)))

        # Robust Mean Absolute Deviation
        robust_mad = np.median(np.abs(img_data_nonzero - np.median(img_data_nonzero)))

        # Root Mean Square
        rms = np.sqrt(np.mean(img_data_nonzero ** 2))

        # Uniformity
        probabilities = hist_nonzero / np.sum(hist_nonzero)
        uniformity = np.sum(probabilities ** 2)

        # Mode (approximate using histogram bin with max count)
        mode_index = np.argmax(hist_nonzero)
        mode_intensity = (bin_edges_nonzero[mode_index] + bin_edges_nonzero[mode_index + 1]) / 2

        os.unlink(tmp.name)
        return {
            'File Name': filename,
            'Mean Intensity': mean_intensity,
            'Median Intensity': median_intensity,
            'Mode Intensity': mode_intensity,
            'Standard Deviation': std_intensity,
            'Variance': variance_intensity,
            'Range': range_intensity,
            'Mean Absolute Deviation': mean_abs_dev,
            'Robust Mean Absolute Deviation': robust_mad,
            'Max Intensity': max_intensity,
            'Min Intensity': min_intensity,
            'Entropy': entropy,
            'Skewness': skewness,
            'Kurtosis': kurtosis,
            'Root Mean Square': rms,
            'Uniformity': uniformity,
            'Contrast Ratio': contrast_ratio,
            '10th Percentile Intensity': percentile_10,
            '90th Percentile Intensity': percentile_90
        }
    except Exception as e:
        print(f"Error processing file {filename}: {e}")
        return None

app = Dash(__name__)

app.layout = html.Div([
    html.H1("View NIfTI Image Properties", style={
        'text-align': 'center',
        'margin-bottom': '20px'
    }),
    dcc.Upload(
        id='upload-data',
        children=html.Button('Select Files', style={
            'width': '100%',
            'height': '60px',
            'fontSize': '20px',
            'backgroundColor': 'steelblue',  # Changed color to a softer blue
            'color': 'white',
            'border': 'none',
            'borderRadius': '5px',
            'cursor': 'pointer',
            'margin-bottom': '20px'
        }),
        multiple=True
    ),
    # ===== ADDED: toggle for original vs z-score =====
    html.Div([
        html.Span("Value scale: ", style={'margin-right': '10px', 'fontWeight': 'bold'}),
        dcc.RadioItems(
            id='value-scale',
            options=[
                {'label': 'Original', 'value': 'original'},
                {'label': 'Z-score', 'value': 'zscore'}
            ],
            value='original',
            inline=True
        )
    ], style={'margin': '10px 0 20px 0'}),
    # =================================================
    html.Div(id='output-data-upload', style={
        'margin-top': '20px',
        'margin-bottom': '20px'
    }),
    dcc.Graph(id='line-plot', style={'height': '500px'}),
    dash_table.DataTable(
        id='data-table',
        page_size=10,
        style_table={'height': '500px', 'overflowY': 'auto'},
        sort_action="native",
        style_header={
            'backgroundColor': 'rgb(230, 230, 230)',
            'fontWeight': 'bold'
        },
        style_cell={
            'textAlign': 'left',
            'padding': '5px',
            'font-family': 'Arial, sans-serif',
            'fontSize': '12px'
        },
        css=[{
            "selector": ".previous-page, .next-page",
            "rule": "float: left; margin-right: 10px;"
        }]
    ),
    html.Button("Download CSV", id="btn_csv", n_clicks=0, style={
        'margin-top': '10px',
        'backgroundColor': 'green',
        'color': 'white',
        'fontSize': '16px',
        'padding': '10px 20px',
        'border': 'none',
        'borderRadius': '5px',
        'cursor': 'pointer'
    }),
    dcc.Download(id="download-dataframe-csv")
], style={
    'backgroundColor': 'lightgrey',  # Changed background to light grey
    'padding': '20px',
    'font-family': 'Arial, sans-serif'
})

@app.callback(
    [Output('output-data-upload', 'children'),
     Output('line-plot', 'figure'),
     Output('data-table', 'data'),
     Output('data-table', 'columns')],
    [Input('upload-data', 'contents'),
     Input('value-scale', 'value')],  # ADDED
    [State('upload-data', 'filename')]
)
def update_output(list_of_contents, scale_mode, list_of_names):  # ADDED scale_mode
    if list_of_contents is not None:
        data = []
        for c, n in zip(list_of_contents, list_of_names):
            stats_data = load_nifti_data(c, n)
            if stats_data:
                data.append(stats_data)
        df = pd.DataFrame(data)
        if df.empty:
            return 'No NIfTI files found.', {}, [], []

        # ===== ADDED: apply z-score per property across files when selected =====
        if scale_mode == 'zscore':
            numeric_cols = [c for c in df.columns if c != 'File Name']
            for col in numeric_cols:
                col_std = df[col].std()
                if pd.notna(col_std) and col_std != 0:
                    df[col] = (df[col] - df[col].mean()) / col_std
                else:
                    df[col] = 0.0
        # ======================================================================

        df_melted = df.melt(
            id_vars=['File Name'],
            var_name='Property',
            value_name='Value'
        )
        line_fig = px.line(
            df_melted,
            x='Property',
            y='Value',
            color='File Name',
            markers=True,  # Added markers to the plot
            template='plotly_white'  # Improved plot appearance
        )
        line_fig.update_layout(
            title='Properties of NIfTI Files' + (' (Z-score)' if scale_mode == 'zscore' else ''),
            xaxis_title='Property',
            yaxis_title='Z-score' if scale_mode == 'zscore' else 'Value',
            xaxis=dict(showgrid=True),  # Show gridlines on x-axis
            yaxis=dict(showgrid=True),  # Show gridlines on y-axis
            font=dict(
                family='Arial, sans-serif',
                size=12,
                color='#333'
            ),
            showlegend=False  # Legend remains removed as per previous instruction
        )
        columns = [{"name": i, "id": i} for i in df.columns]
        data = df.to_dict('records')
        return 'Data loaded successfully', line_fig, data, columns
    return 'Upload NIfTI files to view properties.', {}, [], []

# Callback to handle CSV download
@app.callback(
    Output("download-dataframe-csv", "data"),
    Input("btn_csv", "n_clicks"),
    State("data-table", "derived_virtual_data"),
    prevent_initial_call=True
)
def download_csv(n_clicks, derived_virtual_data):
    if n_clicks:
        df = pd.DataFrame(derived_virtual_data)
        return dcc.send_data_frame(df.to_csv, "nifti_data.csv")
    return dash.no_update

if __name__ == '__main__':
    import sys
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8050
    host = sys.argv[2] if len(sys.argv) > 2 else '127.0.0.1'
    app.run_server(debug=True, host=host, port=port)
