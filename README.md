# CurtainUtils

A utility package for preprocessing and uploading processed and analyzed mass spectrometry-based proteomics data to [Curtain](https://curtain.proteo.info) and [CurtainPTM](https://curtainptm.proteo.info) visualization platforms.

> **What is Curtain?** Curtain is a web-based visualization tool for proteomics data that allows interactive exploration of protein expression data.

> **What is CurtainPTM?** CurtainPTM extends Curtain's functionality to visualize post-translational modifications (PTMs) in proteomics data.

## Installation

### Requirements

- Python 3.6 or higher
- pip package manager

### Install from PyPI

```bash
pip install curtainutils
```

### Install from source

```bash
pip install git+https:///github.com/noatgnu/curtainutils.git
```

## Conversion to CurtainPTM upload format

### Convert MSFragger PTM single site output

```Bash
msf-curtainptm -f msfragger_output.txt -i "Index" -o curtainptm_input.txt -p "Peptide" -a proteome.fasta
```

<table>
<tr><td>Parameter</td><td>Description</td></tr>
<tr><td>-f</td><td>MSFragger PTM output file containing differential analysis</td></tr>
<tr><td>-i</td><td>Column name containing site information (with accession ID and PTM position)</td></tr>
<tr><td>-o</td><td>Output file name for CurtainPTM format</td></tr>
<tr><td>-p</td><td>Column name containing peptide sequences</td></tr>
<tr><td>-a</td><td>FASTA file for protein sequence reference</td></tr>
</table>

### Convert DIA-NN PTM output

```Bash
diann-curtainptm -p diann_differential.txt -r diann_report.txt -o curtainptm_input.txt -m "Phospho"
```

<table>
<tr><td>Parameter</td><td>Description</td></tr>
<tr><td>-p</td><td>Differential analysis file containing Modified.Sequence, Precursor.Id, Protein.Group</td></tr>
<tr><td>-r</td><td>DIA-NN report file containing protein sequences</td></tr>
<tr><td>-o</td><td>Output file name for CurtainPTM format</td></tr>
<tr><td>-m</td><td>Modification type (e.g., Phospho, Acetyl, Methyl, etc.)</td></tr>
</table>

### Convert Spectronaut output

```Bash
spn-curtainptm -f spectronaut_data.txt -o curtain_input.txt
```

<table>
<tr><td>Parameter</td><td>Description</td></tr>
<tr><td>-f</td><td>Spectronaut output file containing differential analysis</td></tr>
<tr><td>-o</td><td>Output file name for CurtainPTM format</td></tr>
</table>

## API Intergration

### Upload to Curtain backend

```py
from curtainutils.client import CurtainClient, add_imputation_map, create_imputation_map, add_uniprot_data

# Initialize client
client = CurtainClient("https://your-curtain-server.com") # Optional api_key parameters

# Define parameters
de_file = "differential_data.txt"
raw_file = "raw_data.txt"
fc_col = "log2FC"
p_col = "p_value"
primary_id_de_col = "Protein"
primary_id_raw_col = "Protein"
sample_cols = ["Sample1.1", "Sample1.2", "Sample1.3", "Sample2.1", "Sample2.2", "Sample2.3"]
description = "My protein analysis"
# Create payload
payload = client.create_curtain_session_payload(
    de_file=de_file,
    raw_file=raw_file,
    fc_col=fc_col,
    transform_fc=False,  # Set to True if fold change needs log transformation
    transform_significant=False,  # Set to True if p-values need -log10 transformation
    reverse_fc=False,  # Set to True to reverse fold change direction
    p_col=p_col,
    comp_col="",  # Optional comparison column
    comp_select=[],  # Optional comparison values to select
    primary_id_de_col=primary_id_de_col,
    primary_id_raw_col=primary_id_raw_col,
    sample_cols=sample_cols,
    description=description
)

# Optional: Add uniprot data
add_uniprot_data(payload, raw_file)

# Optional: Add imputation map
imputation_file = "imputed_data.txt" 
imputation_map = create_imputation_map(imputation_file, primary_id_raw_col, sample_cols)
add_imputation_map(payload, imputation_map)

# Submit to server
package = {
    "enable": "True",
    "description": description,
    "curtain_type": "TP",
  "permanent": "False",
}
link_id = client.post_curtain_session(package, payload)
print(f"Access your visualization at: https:/frontend/#/{link_id}")
```

### Upload to CurtainPTM backend (PTM-specific)

```py
from curtainutils.client import CurtainClient, add_imputation_map, create_imputation_map, add_uniprot_data_ptm

# Initialize client for CurtainPTM
client = CurtainClient("https://your-curtain-server.com") # Optional api_key parameter

# Define PTM-specific parameters
de_file = "ptm_differential_data.txt"  # Must contain PTM-specific columns
raw_file = "ptm_raw_data.txt"
fc_col = "log2FoldChange"
p_col = "adj.P.Val"
primary_id_de_col = "Unique identifier"  # Primary ID in differential data
primary_id_raw_col = "T: Unique identifier"  # Primary ID in raw data (PTM format)
sample_cols = ["Control.1", "Control.2", "Control.3", "Treatment.1", "Treatment.2", "Treatment.3"]

# PTM-specific columns
peptide_col = "Phospho (STY) Probabilities"  # Peptide sequences
acc_col = "Protein"  # UniProt accessions
position_col = "Position"  # PTM positions in protein
position_in_peptide_col = "Position in peptide"  # PTM positions in peptide
sequence_window_col = "Sequence window"  # Protein sequence windows
score_col = "Localization prob"  # PTM localization scores

description = "My PTM analysis"

# Create CurtainPTM payload
payload = client.create_curtain_ptm_session_payload(
    de_file=de_file,
    raw_file=raw_file,
    fc_col=fc_col,
    transform_fc=False,  # Set to True if fold change needs log2 transformation
    transform_significant=False,  # Set to True if p-values need -log10 transformation
    reverse_fc=False,  # Set to True to reverse fold change direction
    p_col=p_col,
    primary_id_de_col=primary_id_de_col,
    primary_id_raw_col=primary_id_raw_col,
    sample_cols=sample_cols,
    peptide_col=peptide_col,
    acc_col=acc_col,
    position_col=position_col,
    position_in_peptide_col=position_in_peptide_col,
    sequence_window_col=sequence_window_col,
    score_col=score_col,
    description=description
)

# Add PTM-specific UniProt data and mappings
add_uniprot_data_ptm(payload, raw_file, de_file)

# Optional: Add imputation map
imputation_map = create_imputation_map(raw_file, primary_id_raw_col, sample_cols)
add_imputation_map(payload, imputation_map)

# Submit to CurtainPTM server
package = {
    "enable": "True",
    "description": description,
    "curtain_type": "PTM",  # Important: Set to PTM for CurtainPTM
    "permanent": "False",
}
link_id = client.post_curtain_session(package, payload)
print(f"Access your PTM visualization at: https://curtainptm.proteo.info/#/{link_id}")
```

### CurtainPTM-specific parameters

<table>
<tr><td>Parameter</td><td>Description</td></tr>
<tr><td>peptide_col</td><td>Column name containing peptide sequences</td></tr>
<tr><td>acc_col</td><td>Column name containing UniProt accession IDs</td></tr>
<tr><td>position_col</td><td>Column name containing PTM positions in protein sequence</td></tr>
<tr><td>position_in_peptide_col</td><td>Column name containing PTM positions in peptide sequence</td></tr>
<tr><td>sequence_window_col</td><td>Column name containing protein sequence windows</td></tr>
<tr><td>score_col</td><td>Column name containing PTM localization scores</td></tr>
<tr><td>curtain_type</td><td>Must be set to "PTM" for CurtainPTM submissions</td></tr>
</table>

**Important Notes for CurtainPTM:**
- Use `add_uniprot_data_ptm()` instead of `add_uniprot_data()` for PTM data
- The function requires both raw and differential data files for proper accession mapping
- Raw data primary ID column typically has "T: " prefix (e.g., "T: Unique identifier")
- PTM-specific columns are required for proper visualization
- Set `curtain_type` to "PTM" in the submission package

### Common API payload creation parameters

<table>
<tr><td>Parameter</td><td>Description</td></tr>
<tr><td>de_file</td><td>Path to differential expression file</td></tr>
<tr><td>raw_file</td><td>Path to raw data file</td></tr>
<tr><td>fc_col</td><td>Column name containing fold change values</td></tr>
<tr><td>transform_fc</td><td>Whether fold change values need log transformation</td></tr>
<tr><td>p_col</td><td>Column name containing significance/p-values</td></tr>
<tr><td>primary_id_de_col</td><td>ID column name in differential expression file</td></tr>
<tr><td>primary_id_raw_col</td><td>ID column name in raw data file</td></tr>
<tr><td>sample_cols</td><td>List of column names containing sample data</td></tr>
</table>

## Plot Customization

CurtainUtils provides comprehensive functions to customize the appearance and behavior of plots in both Curtain and CurtainPTM. All customization functions work by modifying the payload's settings before submission.

### Volcano Plot Customization

Use `configure_volcano_plot()` to customize volcano plot appearance:

```python
from curtainutils.client import configure_volcano_plot

# Basic volcano plot customization
payload = configure_volcano_plot(payload,
    title="My Experiment: Treatment vs Control",
    x_title="Log2 Fold Change",
    y_title="-log10(adjusted p-value)",
    width=1000,
    height=800
)

# Advanced axis configuration
payload = configure_volcano_plot(payload,
    x_min=-5, x_max=5,
    y_min=0, y_max=10,
    x_tick_interval=1,
    y_tick_interval=2,
    show_x_grid=True,
    show_y_grid=False
)

# Custom margins and legend positioning
payload = configure_volcano_plot(payload,
    margin_left=120,
    margin_right=80,
    margin_bottom=100,
    margin_top=80,
    legend_x=0.8,
    legend_y=0.9,
    marker_size=8
)
```

#### Volcano Plot Parameters

<table>
<tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
<tr><td>x_min, x_max, y_min, y_max</td><td>Axis ranges</td><td>Auto</td></tr>
<tr><td>x_title, y_title</td><td>Axis titles</td><td>"Log2FC", "-log10(p-value)"</td></tr>
<tr><td>x_tick_interval, y_tick_interval</td><td>Tick spacing</td><td>Auto</td></tr>
<tr><td>x_tick_length, y_tick_length</td><td>Tick mark length</td><td>5</td></tr>
<tr><td>width, height</td><td>Plot dimensions (pixels)</td><td>800, 1000</td></tr>
<tr><td>margin_left, margin_right, margin_bottom, margin_top</td><td>Plot margins</td><td>Auto</td></tr>
<tr><td>show_x_grid, show_y_grid</td><td>Grid line visibility</td><td>True</td></tr>
<tr><td>title</td><td>Main plot title</td><td>""</td></tr>
<tr><td>legend_x, legend_y</td><td>Legend position (0-1)</td><td>Auto</td></tr>
<tr><td>marker_size</td><td>Point marker size</td><td>10</td></tr>
<tr><td>additional_shapes</td><td>Custom plot annotations</td><td>[]</td></tr>
</table>

### Bar Chart Customization

Use `configure_bar_chart()` to customize bar charts and violin plots:

```python
from curtainutils.client import configure_bar_chart

# Basic bar chart sizing
payload = configure_bar_chart(payload,
    bar_chart_width=50,
    average_bar_chart_width=40,
    violin_plot_width=45
)

# Custom condition colors
condition_colors = {
    'Control': '#4477AA',
    'Treatment_A': '#EE6677', 
    'Treatment_B': '#228833'
}

payload = configure_bar_chart(payload,
    condition_colors=condition_colors,
    violin_point_position=-2  # Show points at an offset on violin plot
)
```

#### Bar Chart Parameters

<table>
<tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
<tr><td>bar_chart_width</td><td>Width per column in individual bar chart</td><td>0 (auto)</td></tr>
<tr><td>average_bar_chart_width</td><td>Width per column in average bar chart</td><td>0 (auto)</td></tr>
<tr><td>violin_plot_width</td><td>Width per column in violin plot</td><td>0 (auto)</td></tr>
<tr><td>profile_plot_width</td><td>Width per column in profile plot</td><td>0 (auto)</td></tr>
<tr><td>condition_colors</td><td>Dict mapping condition names to colors</td><td>{}</td></tr>
<tr><td>violin_point_position</td><td>Point position relative to violin</td><td>-2</td></tr>
</table>

### General Plot Settings

Use `configure_general_plot_settings()` for settings that affect all visualizations:

```python
from curtainutils.client import configure_general_plot_settings

# Font and significance thresholds
payload = configure_general_plot_settings(payload,
    font_family="Arial",
    p_cutoff=0.01,
    fc_cutoff=1.0
)

# Color palette and condition management
default_colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A']
condition_colors = {
    'Control': '#808080',
    'Treatment': '#FF6B6B'
}

payload = configure_general_plot_settings(payload,
    default_colors=default_colors,
    condition_colors=condition_colors,
    condition_order=['Control', 'Treatment'],
    sample_visibility={'Sample1': True, 'Sample2': False}
)
```

#### General Plot Parameters

<table>
<tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
<tr><td>font_family</td><td>Font for all plots</td><td>"Arial"</td></tr>
<tr><td>p_cutoff</td><td>P-value significance threshold</td><td>0.05</td></tr>
<tr><td>fc_cutoff</td><td>Fold change significance threshold</td><td>0.6</td></tr>
<tr><td>default_colors</td><td>Default color palette (list)</td><td>Curtain default</td></tr>
<tr><td>condition_colors</td><td>Condition-specific colors (dict)</td><td>{}</td></tr>
<tr><td>condition_order</td><td>Order of conditions in plots</td><td>[]</td></tr>
<tr><td>sample_visibility</td><td>Show/hide specific samples</td><td>{}</td></tr>
</table>

### PTM-Specific Customization (CurtainPTM Only)

Use `configure_ptm_specific_settings()` for PTM analysis features:

```python
from curtainutils.client import configure_ptm_specific_settings

# PTM database integration
custom_ptm_data = {
    'phosphorylation': {'P12345': {'P12345-1': [{'position': 11, 'residue': 'S'}, {'position': 12, 'residue': 'T'}]}},
}

payload = configure_ptm_specific_settings(payload,
    custom_ptm_data=custom_ptm_data,
    variant_corrections={'P12345': 'P12345-6'},
    custom_sequences={'custom_seq_1': 'PEPTIDESEQUENCE'}
)
```

#### PTM-Specific Parameters

<table>
<tr><th>Parameter</th><th>Description</th><th>Usage</th></tr>
<tr><td>custom_ptm_data</td><td>Custom PTM database annotations</td><td>External PTM database integration</td></tr>
<tr><td>variant_corrections</td><td>PTM position corrections for variants</td><td>Handle protein isoforms</td></tr>
<tr><td>custom_sequences</td><td>Custom peptide sequences</td><td>Non-UniProt sequences</td></tr>
</table>

### Color Schemes

CurtainUtils includes several predefined color schemes:

```python
# Scientific publication colors
scientific_colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
    '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'
]

# Colorblind-friendly palette  
colorblind_friendly = [
    '#0173B2', '#DE8F05', '#029E73', '#CC78BC',
    '#CA9161', '#FBAFE4', '#949494', '#ECE133'
]

# High contrast
high_contrast = [
    '#000000', '#E69F00', '#56B4E9', '#009E73',
    '#F0E442', '#0072B2', '#D55E00', '#CC79A7'
]

payload = configure_general_plot_settings(payload,
    default_colors=colorblind_friendly
)
```