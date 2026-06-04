"""
Image Properties QC Dashboard for NIfTI MRI/PET (Dash)

- Metrics are computed only on the non-zero voxel region.
- Dashboard includes interpretive guidelines and artifact mapping (entropy/skewness/kurtosis etc.)
- Provides actionable QC flags using robust z-scores across the uploaded batch:
    WARN if |robust_z| >= 3
    FAIL if |robust_z| >= 5

Notes:
- Absolute "good/bad" thresholds are NOT stable across sites, tracers, sequences, and preprocessing.
  The dashboard therefore emphasizes within-batch robust outlier detection + a few "hard failure" checks
  (NaN/Inf, extreme zero fraction, extreme mask coverage, clipping).

Dependencies:
- dash, plotly, pandas, numpy, nibabel, scipy, scikit-image
"""

from dash import Dash, dcc, html, Input, Output, State, dash_table, ALL
import dash
import pandas as pd
import numpy as np
import nibabel as nib
import base64
import tempfile
import os
import io
import json
import scipy.stats as stats
from PIL import Image

from typing import Optional

try:
    import pydicom
except Exception:
    pydicom = None

from skimage.filters import threshold_otsu
from skimage.morphology import ball, binary_closing, binary_opening
from skimage.measure import label


# -----------------------------
# Robust statistics helpers
# -----------------------------
def robust_mad(x: np.ndarray) -> float:
    """Median absolute deviation (MAD)."""
    x = np.asarray(x)
    med = np.median(x)
    return np.median(np.abs(x - med))


def robust_zscore(x: pd.Series) -> pd.Series:
    """
    Robust z-score using MAD:
        z = 0.6745 * (x - median) / MAD
    If MAD==0 -> zeros.
    """
    med = x.median()
    mad = robust_mad(x.dropna().values)
    if mad == 0 or np.isnan(mad):
        return pd.Series(np.zeros(len(x)), index=x.index)
    return 0.6745 * (x - med) / mad


def qc_flag_from_rz(rz: float) -> str:
    ar = abs(rz)
    if ar >= 5:
        return "FAIL"
    if ar >= 3:
        return "WARN"
    return "OK"


# -----------------------------
# Metric computation
# -----------------------------
def shannon_entropy_from_hist(vals: np.ndarray, bins: int = 256) -> float:
    vals = vals[np.isfinite(vals)]
    if vals.size < 10:
        return np.nan
    hist, _ = np.histogram(vals, bins=bins)
    hist = hist.astype(np.float64)
    if hist.sum() == 0:
        return np.nan
    p = hist / hist.sum()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def clipping_fraction(vals: np.ndarray) -> float:
    vals = vals[np.isfinite(vals)]
    if vals.size < 10:
        return np.nan
    vmax = np.max(vals)
    if vmax == 0:
        return 0.0
    tol = max(1e-6, 1e-4 * abs(vmax))
    return float(np.mean(np.abs(vals - vmax) <= tol))


def compute_metrics(vals: np.ndarray) -> dict:
    vals = vals[np.isfinite(vals)]
    if vals.size < 50:
        return {k: np.nan for k in [
            "mean", "median", "std", "var", "p10", "p90", "p99",
            "mad_mean", "mad_median", "rms", "entropy",
            "skew", "kurtosis", "uniformity", "snr_proxy",
            "cv", "range", "clip_frac"
        ]}

    mean = float(np.mean(vals))
    median = float(np.median(vals))
    std = float(np.std(vals))
    var = float(np.var(vals))
    p10 = float(np.percentile(vals, 10))
    p90 = float(np.percentile(vals, 90))
    p99 = float(np.percentile(vals, 99))
    mad_mean = float(np.mean(np.abs(vals - mean)))
    mad_median = float(np.median(np.abs(vals - median)))
    rms = float(np.sqrt(np.mean(vals ** 2)))
    ent = float(shannon_entropy_from_hist(vals, bins=256))
    skew = float(stats.skew(vals)) if vals.size >= 3 else np.nan
    kurt = float(stats.kurtosis(vals)) if vals.size >= 4 else np.nan

    hist, _ = np.histogram(vals, bins=256)
    if hist.sum() == 0:
        uniformity = np.nan
    else:
        p = hist / hist.sum()
        uniformity = float(np.sum(p ** 2))

    snr_proxy = float(mean / std) if std > 0 else np.nan
    cv = float(std / mean) if mean != 0 else np.nan
    vmin, vmax = float(np.min(vals)), float(np.max(vals))
    rng = vmax - vmin
    clip_frac = float(clipping_fraction(vals))

    return {
        "mean": mean,
        "median": median,
        "std": std,
        "var": var,
        "p10": p10,
        "p90": p90,
        "p99": p99,
        "mad_mean": mad_mean,
        "mad_median": mad_median,
        "rms": rms,
        "entropy": ent,
        "skew": skew,
        "kurtosis": kurt,
        "uniformity": uniformity,
        "snr_proxy": snr_proxy,
        "cv": cv,
        "range": rng,
        "clip_frac": clip_frac,
    }


def load_nifti_from_upload(contents: str, filename: str) -> np.ndarray:
    _, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)

    suffix = ".nii.gz" if filename.endswith(".gz") else ".nii"
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        tmp.write(decoded)
        tmp.flush()
        tmp_path = tmp.name

    try:
        img = nib.load(tmp_path)
        data = img.get_fdata(dtype=np.float32)
        if data.ndim == 4 and data.shape[-1] > 1:
            data = np.nanmean(data, axis=-1).astype(np.float32)
        return data
    finally:
        try:
            os.unlink(tmp_path)
        except Exception:
            pass


def extract_qc_row(contents: str, filename: str) -> Optional[dict]:
    try:
        data = load_nifti_from_upload(contents, filename)
        finite = np.isfinite(data)

        naninf_frac = float(np.mean(~finite))
        data2 = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)

        whole_mask = (data2 != 0)
        whole_vals = data2[whole_mask]

        zero_frac = float(np.mean(data2 == 0))
        neg_frac_whole = float(np.mean(whole_vals < 0)) if whole_vals.size else np.nan
        whole_coverage = float(whole_mask.mean())

        wm = compute_metrics(whole_vals)

        row = {
            "File": filename,
            "NaNInf_Frac": naninf_frac,
            "Zero_Frac": zero_frac,
            "WHOLE_Coverage": whole_coverage,
            "WHOLE_NegFrac": neg_frac_whole,
            "WHOLE_Nvox": int(whole_vals.size),
        }

        for k, v in wm.items():
            row[f"WHOLE_{k}"] = v

        return row
    except Exception as e:
        print(f"[ERROR] {filename}: {e}")
        return None




def pretty_metric_name(name: str) -> str:
    if not isinstance(name, str):
        return name
    return name.replace("WHOLE_", "")


# -----------------------------
# Orthogonal slice preview helpers
# -----------------------------
def _normalize_slice_for_png(sl: np.ndarray) -> np.ndarray:
    sl = np.asarray(sl, dtype=np.float32)
    sl = np.nan_to_num(sl, nan=0.0, posinf=0.0, neginf=0.0)
    finite = np.isfinite(sl)
    vals = sl[finite]
    if vals.size == 0:
        return np.zeros(sl.shape, dtype=np.uint8)
    lo, hi = np.percentile(vals, [1, 99])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(vals.min()), float(vals.max())
    if hi <= lo:
        return np.zeros(sl.shape, dtype=np.uint8)
    x = np.clip((sl - lo) / (hi - lo), 0, 1)
    return (x * 255).astype(np.uint8)


def _slice_to_data_uri(sl: np.ndarray) -> str:
    arr = _normalize_slice_for_png(sl)
    img = Image.fromarray(arr, mode="L")
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("utf-8")


def _extract_orthogonal_views(data: np.ndarray) -> dict:
    x = np.asarray(data, dtype=np.float32)
    if x.ndim != 3:
        raise ValueError("Expected a 3D volume after loading")
    i = x.shape[0] // 2
    j = x.shape[1] // 2
    k = x.shape[2] // 2

    axial = np.rot90(x[:, :, k])
    sagittal = np.rot90(x[i, :, :])
    coronal = np.rot90(x[:, j, :])

    return {
        "Axial": _slice_to_data_uri(axial),
        "Sagittal": _slice_to_data_uri(sagittal),
        "Coronal": _slice_to_data_uri(coronal),
    }


def _parse_heatmap_click(click_data: dict) -> tuple:
    if not click_data or "points" not in click_data or not click_data["points"]:
        return None, None, None
    pt = click_data["points"][0]
    metric = pt.get("x")
    row_label = pt.get("y", "")
    rz_val = pt.get("z")
    filename = row_label.split(": ", 1)[1] if isinstance(row_label, str) and ": " in row_label else row_label
    return filename, metric, rz_val


# -----------------------------
# Metadata extraction helpers
# -----------------------------
def _metadata_json_safe(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, bytes):
        return f"<bytes: {len(value)}>"
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_metadata_json_safe(v) for v in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [_metadata_json_safe(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _metadata_json_safe(v) for k, v in value.items()}
    return str(value)


def _decode_upload_to_tempfile(contents: str, filename: str) -> str:
    _, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    suffix = os.path.splitext(filename)[1] or '.bin'
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        tmp.write(decoded)
        tmp.flush()
        return tmp.name


def extract_nifti_metadata_from_upload(contents: str, filename: str) -> dict:
    tmp_path = _decode_upload_to_tempfile(contents, filename)
    try:
        img = nib.load(tmp_path)
        hdr = img.header
        meta = {
            'filename': filename,
            'file_type': 'NIfTI',
            'metadata': {
                'shape': _metadata_json_safe(img.shape),
                'voxel_sizes': _metadata_json_safe(hdr.get_zooms()),
                'datatype': str(hdr.get_data_dtype()),
                'affine': _metadata_json_safe(np.asarray(img.affine)),
            }
        }
        for key in hdr.keys():
            try:
                meta['metadata'][str(key)] = _metadata_json_safe(hdr[key])
            except Exception:
                continue
        return meta
    finally:
        try:
            os.unlink(tmp_path)
        except Exception:
            pass


def extract_dicom_metadata_from_upload(contents: str, filename: str) -> dict:
    if pydicom is None:
        return {
            'filename': filename,
            'file_type': 'DICOM',
            'metadata': {'error': 'pydicom is required to read DICOM metadata in this tab.'}
        }

    tmp_path = _decode_upload_to_tempfile(contents, filename)
    try:
        ds = pydicom.dcmread(tmp_path, stop_before_pixels=True, force=True)
        meta = {'filename': filename, 'file_type': 'DICOM', 'metadata': {}}
        for elem in ds.iterall():
            if elem.keyword == 'PixelData' or elem.VR == 'SQ':
                continue
            key = elem.keyword or str(elem.tag)
            try:
                meta['metadata'][key] = _metadata_json_safe(elem.value)
            except Exception:
                meta['metadata'][key] = str(elem.value)
        return meta
    finally:
        try:
            os.unlink(tmp_path)
        except Exception:
            pass


def extract_metadata_from_upload(contents: str, filename: str) -> dict:
    lower = (filename or '').lower()
    if lower.endswith('.nii') or lower.endswith('.nii.gz'):
        return extract_nifti_metadata_from_upload(contents, filename)
    return extract_dicom_metadata_from_upload(contents, filename)


def metadata_records_for_table(metadata_payload: list) -> list:
    rows = []
    for item in metadata_payload or []:
        filename = item.get('filename', '')
        file_type = item.get('file_type', '')
        metadata = item.get('metadata', {}) or {}
        for key, value in metadata.items():
            rows.append({
                'File': filename,
                'Type': file_type,
                'Field': key,
                'Value': json.dumps(value) if isinstance(value, (dict, list)) else str(value)
            })
    return rows


# -----------------------------
# Interpretive guidelines text
# -----------------------------
GUIDE_MD = r"""
### What the dashboard computes

For every uploaded NIfTI, each metric is computed from the non-zero voxel region, defined here as **all non-zero voxels** in the image.

---

### How to use the metrics in practice

Because MRI/PET intensity scales vary widely across scanners, protocols, tracers, reconstruction settings, and preprocessing steps, **absolute thresholds for entropy/skewness/kurtosis are not reliable** across studies.  
This dashboard therefore uses a **batch-based robust outlier rule**:

- **WARN** if **|robust z| ≥ 3**  
- **FAIL** if **|robust z| ≥ 5**

Robust z-score is computed using **median and MAD** across the uploaded set, per metric.

Use this workflow:

1) Upload a batch from the *same modality/protocol*.
2) Look at the **QC Heatmap** and **Flagged Table**.
3) Inspect files that are outliers on multiple metrics.
4) Use the orthogonal slice preview to visually inspect flagged scans.

---

### Interpretation: what entropy / skewness / kurtosis usually mean (first-order)

These are computed from the **intensity histogram** of the non-zero voxel region.

#### Entropy (higher = broader / more uniform intensity distribution)
**Often elevated when:**
- Increased noise
- Motion + blurring broadens the intensity distribution
- Poor normalization / wrong scaling spreads values

**Often very low when:**
- Near-empty/blank images, heavy thresholding, or extreme clipping/saturation
- Wrong file (mask, label map, mostly zeros)

#### Skewness (asymmetry of histogram)
**Positive skew (long right tail)** can indicate:
- Hot voxels / spikes
- Mis-scaling that creates a few very large values
- Partial FOV / cropping where only high-intensity structures remain

**Negative skew** can indicate:
- Unexpected negative values
- Over-aggressive bias correction or intensity shifting

#### Kurtosis (tail heaviness / outliers)
**High kurtosis** often indicates:
- Presence of outliers
- Quantization / clipping artifacts that pile up mass + tails

**Very low kurtosis** can happen in overly uniform/flattened distributions.

---

### Practical “hard” QC checks included

- **NaNInf_Frac**: any non-finite values are a pipeline failure.
- **Zero_Frac / Coverage**: extreme zeros or coverage suggests cropping, padding, wrong orientation, or wrong file.
- **clip_frac**: large fraction of voxels at the maximum value suggests saturation/clipping/rescale problems.
- **NegFrac**: many negatives in PET-like data is suspicious.
"""


# -----------------------------
# Dash app
# -----------------------------
app = Dash(__name__)
app.config.suppress_callback_exceptions = True
app.title = "NIfTI QC: Image Properties Dashboard"

app.layout = html.Div(
    style={"fontFamily": "Arial, sans-serif", "backgroundColor": "#f4f4f4", "padding": "20px"},
    children=[
        html.H2("NIfTI Image Properties QC Dashboard (MRI/PET)", style={"textAlign": "center"}),

        dcc.Upload(
            id="upload-data",
            children=html.Button(
                "Select NIfTI files (.nii / .nii.gz)",
                style={
                    "width": "100%", "height": "60px", "fontSize": "18px",
                    "backgroundColor": "steelblue", "color": "white",
                    "border": "none", "borderRadius": "6px", "cursor": "pointer"
                },
            ),
            multiple=True
        ),

        html.Div(style={"marginTop": "10px", "display": "flex", "gap": "16px", "flexWrap": "wrap"}, children=[
            html.Div(style={"minWidth": "260px"}, children=[
                html.Label("Metric view", style={"fontWeight": "bold"}),
                dcc.RadioItems(
                    id="metric-view",
                    options=[
                        {"label": "Raw values", "value": "raw"},
                        {"label": "Robust z-scores (recommended)", "value": "rz"},
                    ],
                    value="rz",
                    inline=True
                ),
            ]),
            html.Div(style={"minWidth": "360px"}, children=[
                html.Label("Outlier thresholds (robust z)", style={"fontWeight": "bold"}),
                dcc.Slider(
                    id="rz-threshold",
                    min=2, max=6, step=0.5, value=3,
                    marks={2: "2", 3: "3 (WARN)", 5: "5 (FAIL)", 6: "6"},
                    tooltip={"placement": "bottom", "always_visible": False},
                ),
            ]),
            html.Div(style={"minWidth": "360px"}, children=[
                html.Label("Pick a metric for the trend plot", style={"fontWeight": "bold"}),
                dcc.Dropdown(id="metric-dropdown", placeholder="Upload files first…"),
            ]),
        ]),

        dcc.Tabs(id="tabs", value="tab-qc", children=[
            dcc.Tab(label="QC Overview", value="tab-qc"),
            dcc.Tab(label="Interpretation & Guidelines", value="tab-guide"),
            dcc.Tab(label="Full Table", value="tab-table"),
            dcc.Tab(label="Metadata", value="tab-metadata"),
        ]),

        # ADDED: global top-right download button (exists in initial layout)
        html.Div(
            style={"display": "flex", "justifyContent": "flex-end", "marginTop": "10px"},
            children=[
                html.Button(
                    "Download CSV",
                    id="btn-download-fulltable",
                    style={
                        "display": "none",  # will be enabled only on Full Table tab
                        "backgroundColor": "steelblue",
                        "color": "white",
                        "border": "none",
                        "borderRadius": "6px",
                        "padding": "10px 14px",
                        "cursor": "pointer",
                        "fontWeight": "bold"
                    }
                )
            ]
        ),

        html.Div(id="tab-content", style={"marginTop": "14px"}),

        dcc.Store(id="store-df"),
        dcc.Store(id="store-uploads"),
        dcc.Store(id="metadata-store"),
        dcc.Download(id="download-fulltable-csv"),
        dcc.Download(id="download-metadata-json"),
    ]
)


# -----------------------------
# Callbacks
# -----------------------------
@app.callback(
    Output("store-df", "data"),
    Output("store-uploads", "data"),
    Output("metric-dropdown", "options"),
    Output("metric-dropdown", "value"),
    Input("upload-data", "contents"),
    State("upload-data", "filename"),
    prevent_initial_call=True
)
def compute_all_files(list_of_contents, list_of_names):
    if not list_of_contents or not list_of_names:
        return None, None, [], None

    rows = []
    for c, n in zip(list_of_contents, list_of_names):
        r = extract_qc_row(c, n)
        if r is not None:
            rows.append(r)

    if not rows:
        return None, None, [], None

    df = pd.DataFrame(rows)

    numeric_cols = [c for c in df.columns if c != "File" and pd.api.types.is_numeric_dtype(df[c])]
    default_metric = "WHOLE_entropy" if "WHOLE_entropy" in numeric_cols else (numeric_cols[0] if numeric_cols else None)

    options = [{"label": pretty_metric_name(col), "value": col} for col in numeric_cols]
    uploads_payload = [{"filename": n, "contents": c} for c, n in zip(list_of_contents, list_of_names)]
    return df.to_json(date_format="iso", orient="split"), uploads_payload, options, default_metric


# ADDED: show/hide the download button depending on the selected tab
@app.callback(
    Output("btn-download-fulltable", "style"),
    Input("tabs", "value"),
)
def toggle_download_button(tab_value):
    base = {
        "backgroundColor": "steelblue",
        "color": "white",
        "border": "none",
        "borderRadius": "6px",
        "padding": "10px 14px",
        "cursor": "pointer",
        "fontWeight": "bold"
    }
    if tab_value == "tab-table":
        base["display"] = "inline-block"
    else:
        base["display"] = "none"
    return base


@app.callback(
    Output("tab-content", "children"),
    Input("tabs", "value"),
    Input("store-df", "data"),
    Input("store-uploads", "data"),
    Input("metadata-store", "data"),
    Input("metric-view", "value"),
    Input("metric-dropdown", "value"),
    Input("rz-threshold", "value"),
    Input({"type": "qc-heatmap", "index": ALL}, "clickData"),
)
def render_tabs(tab, df_json, uploads_payload, metadata_payload, metric_view, metric_col, rz_thr, heatmap_click_list):
    if tab == "tab-guide":
        return html.Div(
            style={"backgroundColor": "white", "padding": "16px", "borderRadius": "10px"},
            children=[dcc.Markdown(GUIDE_MD)]
        )

    if tab == "tab-metadata":
        table_rows = metadata_records_for_table(metadata_payload)
        return html.Div(
            style={"backgroundColor": "white", "padding": "16px", "borderRadius": "10px"},
            children=[
                html.H4("Metadata Viewer"),
                html.Div(
                    "Import one or more DICOM files to inspect any available metadata fields and optionally export them as JSON.",
                    style={"marginBottom": "12px"}
                ),
                dcc.Upload(
                    id={"type": "metadata-upload", "index": 0},
                    children=html.Button(
                        "Select DICOM files for metadata review",
                        style={
                            "width": "100%", "height": "52px", "fontSize": "16px",
                            "backgroundColor": "steelblue", "color": "white",
                            "border": "none", "borderRadius": "6px", "cursor": "pointer"
                        },
                    ),
                    multiple=True
                ),
                html.Div(
                    style={"display": "flex", "justifyContent": "flex-end", "marginTop": "10px", "marginBottom": "10px"},
                    children=[
                        html.Button(
                            "Download metadata JSON",
                            id={"type": "metadata-download-btn", "index": 0},
                            style={
                                "backgroundColor": "steelblue",
                                "color": "white",
                                "border": "none",
                                "borderRadius": "6px",
                                "padding": "10px 14px",
                                "cursor": "pointer",
                                "fontWeight": "bold"
                            }
                        )
                    ]
                ),
                dash_table.DataTable(
                    data=table_rows,
                    columns=[
                        {"name": "File", "id": "File"},
                        {"name": "Type", "id": "Type"},
                        {"name": "Field", "id": "Field"},
                        {"name": "Value", "id": "Value"},
                    ],
                    page_size=15,
                    sort_action="native",
                    filter_action="native",
                    style_table={"overflowX": "auto", "height": "600px", "overflowY": "auto"},
                    style_cell={"textAlign": "left", "fontSize": "11px", "padding": "6px", "minWidth": "120px", "whiteSpace": "normal", "height": "auto"},
                    style_header={"fontWeight": "bold", "backgroundColor": "#eeeeee"},
                ) if table_rows else html.Div(
                    "No metadata loaded yet.",
                    style={"marginTop": "8px", "fontSize": "12px", "color": "#444"}
                )
            ]
        )

    if df_json is None:
        return html.Div(
            style={"backgroundColor": "white", "padding": "16px", "borderRadius": "10px"},
            children=[
                html.H4("Upload NIfTI files to compute QC metrics."),
                html.Div("Tip: upload a batch from the same protocol/modality for meaningful outlier detection.")
            ]
        )

    df = pd.read_json(df_json, orient="split")

    num_cols = [c for c in df.columns if c != "File" and pd.api.types.is_numeric_dtype(df[c])]
    rz_df = df.copy()
    for c in num_cols:
        rz_df[c] = robust_zscore(df[c])

    abs_rz = rz_df[num_cols].abs()
    flagged_mask = (abs_rz >= float(rz_thr)).any(axis=1)
    flagged = df.loc[flagged_mask, ["File"] + num_cols].copy()

    reasons = []
    for idx in flagged.index:
        s = abs_rz.loc[idx].sort_values(ascending=False).head(5)
        parts = [f"{pretty_metric_name(m)} (rz={rz_df.loc[idx, m]:.2f})" for m in s.index]
        reasons.append("; ".join(parts))
    if len(flagged) > 0:
        flagged.insert(1, "Top_Outlier_Metrics", reasons)

    if metric_col is None or metric_col not in num_cols:
        metric_col = "WHOLE_entropy" if "WHOLE_entropy" in num_cols else num_cols[0]

    plot_df = df[["File", metric_col]].copy()
    plot_df["robust_z"] = rz_df[metric_col]

    y_col = "robust_z" if metric_view == "rz" else metric_col
    display_metric_col = pretty_metric_name(metric_col)
    y_title = f"{display_metric_col} (robust z)" if metric_view == "rz" else f"{display_metric_col} (raw)"

    import plotly.express as px

    # CHANGED: trend plot is a line (connected) instead of dots only
    fig_trend = px.line(
        plot_df,
        x="File",
        y=y_col,
        hover_data=["File", metric_col, "robust_z"],
        title=f"Per-file trend: {y_title}",
        markers=True
    )
    fig_trend.update_layout(
        xaxis_title="File",
        yaxis_title=y_title,
        xaxis={"tickangle": 45},
        height=420,
        margin={"l": 40, "r": 20, "t": 60, "b": 120},
        showlegend=False
    )
    if metric_view == "rz":
        fig_trend.add_hline(y=float(rz_thr), line_dash="dash")
        fig_trend.add_hline(y=-float(rz_thr), line_dash="dash")
        fig_trend.add_hline(y=5, line_dash="dot")
        fig_trend.add_hline(y=-5, line_dash="dot")

    heat_cols = [c for c in num_cols if c.startswith("WHOLE_")]
    heat_cols = heat_cols[:40]

    preview_panel = html.Div(
        "Click a z-score cell in the heatmap to view axial, sagittal, and coronal slices.",
        style={"marginTop": "12px", "fontSize": "12px", "color": "#444"}
    )

    if uploads_payload:
        heatmap_click = None
        if isinstance(heatmap_click_list, list) and len(heatmap_click_list) > 0:
            heatmap_click = heatmap_click_list[0]
        sel_filename, sel_metric, sel_rz = _parse_heatmap_click(heatmap_click)
        if sel_filename:
            match = next((u for u in uploads_payload if u.get("filename") == sel_filename), None)
            if match is not None:
                try:
                    vol = load_nifti_from_upload(match["contents"], match["filename"])
                    views = _extract_orthogonal_views(vol)
                    rz_text = f"{float(sel_rz):.2f}" if sel_rz is not None and np.isfinite(sel_rz) else str(sel_rz)
                    preview_panel = html.Div(
                        style={"marginTop": "12px"},
                        children=[
                            html.Div([
                                html.Div(f"Selected file: {sel_filename}", style={"fontWeight": "bold"}),
                                html.Div(f"Metric: {pretty_metric_name(sel_metric)}"),
                                html.Div(f"Robust z: {rz_text}"),
                            ], style={"marginBottom": "10px"}),
                            html.Div(
                                style={"display": "grid", "gridTemplateColumns": "repeat(3, minmax(0, 1fr))", "gap": "12px"},
                                children=[
                                    html.Div(
                                        style={"border": "1px solid #dddddd", "borderRadius": "8px", "padding": "10px", "backgroundColor": "#fafafa"},
                                        children=[
                                            html.Div(name, style={"fontWeight": "bold", "marginBottom": "8px", "textAlign": "center"}),
                                            html.Img(src=src, style={"width": "100%", "display": "block"})
                                        ]
                                    )
                                    for name, src in views.items()
                                ]
                            )
                        ]
                    )
                except Exception:
                    preview_panel = html.Div(
                        "Unable to generate orthogonal views for the selected file.",
                        style={"marginTop": "12px", "color": "crimson"}
                    )

    heat_display_cols = [pretty_metric_name(c) for c in heat_cols]

    hm = px.imshow(
        rz_df[heat_cols].values,
        labels=dict(x="Metric", y="File index", color="robust z"),
        x=heat_display_cols,
        y=[f"{i}: {f}" for i, f in enumerate(df["File"].tolist())],
        title="QC heatmap (robust z-scores) — scan for rows with many extremes"
    )
    hm.update_layout(height=520, margin={"l": 40, "r": 20, "t": 60, "b": 80})

    if tab == "tab-qc":
        return html.Div(children=[
            html.Div(
                style={"backgroundColor": "white", "padding": "16px", "borderRadius": "10px", "marginBottom": "12px"},
                children=[
                    html.H4("QC Overview"),
                    html.Div([
                        html.Div(f"Files uploaded: {len(df)}"),
                        html.Div(f"Flag rule: |robust z| ≥ {rz_thr} (WARN); |robust z| ≥ 5 (FAIL)")
                    ], style={"marginBottom": "10px"}),

                    dcc.Graph(id={"type": "qc-heatmap", "index": 0}, figure=hm, clear_on_unhover=True),

                    preview_panel,

                    html.H5("Flagged files (any metric beyond threshold)"),
                    dash_table.DataTable(
                        data=flagged.to_dict("records") if len(flagged) else [],
                        columns=[{"name": pretty_metric_name(c), "id": c} for c in flagged.columns] if len(flagged) else [],
                        page_size=8,
                        sort_action="native",
                        style_table={"overflowX": "auto"},
                        style_cell={"textAlign": "left", "fontSize": "12px", "padding": "6px"},
                        style_header={"fontWeight": "bold", "backgroundColor": "#eeeeee"},
                    ),
                ]
            ),

            html.Div(
                style={"backgroundColor": "white", "padding": "16px", "borderRadius": "10px"},
                children=[
                    html.H4("Metric trend plot"),
                    dcc.Graph(figure=fig_trend),
                    html.Div(
                        "Tip: Outliers across multiple metrics are more likely to reflect true image quality, scaling, clipping, or preprocessing problems.",
                        style={"fontSize": "12px", "color": "#333"}
                    )
                ]
            )
        ])

    if tab == "tab-table":
        full = df.copy()
        for c in num_cols:
            full[f"{c}__rz"] = rz_df[c]
            full[f"{c}__flag"] = [qc_flag_from_rz(v) if np.isfinite(v) else "NA" for v in rz_df[c].values]

        return html.Div(
            style={"backgroundColor": "white", "padding": "16px", "borderRadius": "10px"},
            children=[
                html.H4("Full metrics table (raw + robust z + flags)"),
                dash_table.DataTable(
                    data=full.to_dict("records"),
                    columns=[{"name": pretty_metric_name(c), "id": c} for c in full.columns],
                    page_size=12,
                    sort_action="native",
                    filter_action="native",
                    style_table={"overflowX": "auto", "height": "600px", "overflowY": "auto"},
                    style_cell={"textAlign": "left", "fontSize": "11px", "padding": "6px", "minWidth": "110px"},
                    style_header={"fontWeight": "bold", "backgroundColor": "#eeeeee"},
                ),
                html.Div(
                    "Use the column filters to quickly find FAIL/WARN flags (e.g., type 'FAIL' into a __flag column).",
                    style={"marginTop": "8px", "fontSize": "12px"}
                )
            ]
        )

    return html.Div()


@app.callback(
    Output("metadata-store", "data"),
    Input({"type": "metadata-upload", "index": ALL}, "contents"),
    State({"type": "metadata-upload", "index": ALL}, "filename"),
    prevent_initial_call=True
)
def load_metadata_files(contents_list, filenames_list):
    if not contents_list or not filenames_list:
        return None

    flat_contents = contents_list[0] if isinstance(contents_list, list) and len(contents_list) == 1 and isinstance(contents_list[0], list) else contents_list
    flat_names = filenames_list[0] if isinstance(filenames_list, list) and len(filenames_list) == 1 and isinstance(filenames_list[0], list) else filenames_list

    if not flat_contents or not flat_names:
        return None

    payload = []
    for c, n in zip(flat_contents, flat_names):
        try:
            payload.append(extract_metadata_from_upload(c, n))
        except Exception as e:
            payload.append({
                'filename': n,
                'file_type': 'Unknown',
                'metadata': {'error': str(e)}
            })
    return payload


@app.callback(
    Output("download-metadata-json", "data"),
    Input({"type": "metadata-download-btn", "index": ALL}, "n_clicks"),
    State("metadata-store", "data"),
    prevent_initial_call=True
)
def download_metadata_json(n_clicks_list, metadata_payload):
    if not metadata_payload:
        return dash.no_update
    total_clicks = sum([int(x or 0) for x in (n_clicks_list or [])])
    if total_clicks <= 0:
        return dash.no_update
    return {
        "content": json.dumps(metadata_payload, indent=2),
        "filename": "image_metadata.json"
    }


# Download callback (now safe because the button exists in the initial layout)
@app.callback(
    Output("download-fulltable-csv", "data"),
    Input("btn-download-fulltable", "n_clicks"),
    State("store-df", "data"),
    prevent_initial_call=True
)
def download_full_table(n_clicks, df_json):
    if not df_json:
        return dash.no_update

    df = pd.read_json(df_json, orient="split")
    num_cols = [c for c in df.columns if c != "File" and pd.api.types.is_numeric_dtype(df[c])]

    rz_df = df.copy()
    for c in num_cols:
        rz_df[c] = robust_zscore(df[c])

    full = df.copy()
    for c in num_cols:
        full[f"{c}__rz"] = rz_df[c]
        full[f"{c}__flag"] = [qc_flag_from_rz(v) if np.isfinite(v) else "NA" for v in rz_df[c].values]

    return dcc.send_data_frame(full.to_csv, "nifti_qc_full_table.csv", index=False)


if __name__ == "__main__":
    import sys
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8050
    app.run_server(debug=True, port=port)
