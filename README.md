# CurtainUtils

A utility package for preprocessing and uploading processed and analyzed mass spectrometry-based proteomics data to [Curtain](https://curtain.proteo.info) and [CurtainPTM](https://curtainptm.proteo.info) visualization platforms.

> **What is Curtain?** Curtain is a web-based visualization tool for proteomics data that allows interactive exploration of protein expression data.

> **What is CurtainPTM?** CurtainPTM extends Curtain's functionality to visualize post-translational modifications (PTMs) in proteomics data.

## Installation

### Requirements

- Python 3.10 or higher
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

client = CurtainClient("https://your-curtain-server.com")

payload = client.create_curtain_session_payload(
    de_file="differential_data.txt",
    raw_file="raw_data.txt",
    fc_col="log2FC",
    transform_fc=False,
    transform_significant=False,
    reverse_fc=False,
    p_col="p_value",
    comp_col="",
    comp_select=[],
    primary_id_de_col="Protein",
    primary_id_raw_col="Protein",
    sample_cols=["Sample1.1", "Sample1.2", "Sample1.3", "Sample2.1", "Sample2.2", "Sample2.3"],
    description="My protein analysis"
)

add_uniprot_data(payload, "raw_data.txt")

imputation_map = create_imputation_map("imputed_data.txt", "Protein", sample_cols)
add_imputation_map(payload, imputation_map)

package = {
    "enable": "True",
    "description": "My protein analysis",
    "curtain_type": "TP",
    "permanent": "False",
}
link_id = client.post_curtain_session(package, payload)
print(f"Access your visualization at: https:/frontend/#/{link_id}")
```

### Upload to CurtainPTM backend

```py
from curtainutils.client import CurtainClient, add_imputation_map, create_imputation_map, add_uniprot_data_ptm

client = CurtainClient("https://your-curtain-server.com")

payload = client.create_curtain_ptm_session_payload(
    de_file="ptm_differential_data.txt",
    raw_file="ptm_raw_data.txt",
    fc_col="log2FoldChange",
    transform_fc=False,
    transform_significant=False,
    reverse_fc=False,
    p_col="adj.P.Val",
    primary_id_de_col="Unique identifier",
    primary_id_raw_col="T: Unique identifier",
    sample_cols=["Control.1", "Control.2", "Control.3", "Treatment.1", "Treatment.2", "Treatment.3"],
    peptide_col="Phospho (STY) Probabilities",
    acc_col="Protein",
    position_col="Position",
    position_in_peptide_col="Position in peptide",
    sequence_window_col="Sequence window",
    score_col="Localization prob",
    description="My PTM analysis"
)

add_uniprot_data_ptm(payload, "ptm_raw_data.txt", "ptm_differential_data.txt")
imputation_map = create_imputation_map("ptm_raw_data.txt", "T: Unique identifier", sample_cols)
add_imputation_map(payload, imputation_map)

package = {
    "enable": "True",
    "description": "My PTM analysis",
    "curtain_type": "PTM",
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

**Notes for CurtainPTM:**
- Use `add_uniprot_data_ptm()` instead of `add_uniprot_data()` for PTM data
- Raw data primary ID column typically has "T: " prefix (e.g., "T: Unique identifier")
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

## Search Groups and Data Point Selection

CurtainUtils provides search group functionality to create custom selections of proteins or data points and assign them colors for visualization. Search groups allow you to highlight specific proteins of interest across all visualizations.

### Creating Search Groups

#### Basic Search Groups

```python
from curtainutils.client import add_search_group

add_search_group(
    payload=payload,
    group_name="Kinases of Interest", 
    gene_names=["AKT1", "MTOR", "PIK3CA", "PTEN"],
    color="#FF6B6B",
    curtain_type="TP"
)

add_search_group(
    payload=payload,
    group_name="Key Targets",
    protein_ids=["P31749", "P42345", "Q9Y243"],
    color="#4ECDC4",
    curtain_type="TP"
)

add_search_group(
    payload=payload,
    group_name="Mixed Selection",
    gene_names=["TP53", "MYC"], 
    protein_ids=["P04637"],
    color="#45B7D1",
    curtain_type="TP"
)
```

#### Advanced Search Groups

```python
# Regex-based selection from raw data
payload = create_search_group_from_regex(
    payload=payload,
    group_name="Histone Proteins",
    pattern="HIST.*",
    raw_df=raw_df,
    search_column="Gene_Name", 
    primary_id_col="Primary.IDs",
    color="#9B59B6"
)

# Conditional selection based on data values
payload = create_search_group_from_conditions(
    payload=payload,
    group_name="High FC Proteins",
    raw_df=raw_df,
    primary_id_col="Primary.IDs",
    conditions={
        "Log2_FC": lambda x: abs(x) > 2.0,
        "P_Value": lambda x: x < 0.001
    },
    color="#E74C3C"
)
```

#### Comparison-Specific Groups

```python
add_search_group(
    payload=payload,
    group_name="Treatment Response",
    gene_names=["EGFR", "KRAS", "PIK3CA"],
    color="#F39C12",
    curtain_type="TP"
)

add_search_group(
    payload=payload,
    group_name="PTM Targets",
    protein_ids=["P31749", "Q9Y243"],
    color="#8E44AD",
    curtain_type="PTM"
)
```

### Search Group Parameters

<table>
<tr><th>Parameter</th><th>Description</th><th>Required</th></tr>
<tr><td>group_name</td><td>Name for the search group</td><td>Yes</td></tr>
<tr><td>gene_names</td><td>List of gene names to search for</td><td>No*</td></tr>
<tr><td>protein_ids</td><td>List of protein/accession IDs</td><td>No*</td></tr>
<tr><td>color</td><td>Hex color code for highlighting</td><td>No (auto-assigned)</td></tr>
<tr><td>specific_comparison_label</td><td>Target specific comparison group</td><td>No</td></tr>
</table>

*At least one of `gene_names` or `protein_ids` must be provided.

### Managing Multiple Search Groups

```python
groups = [
    {"name": "Oncogenes", "genes": ["MYC", "KRAS", "PIK3CA"], "color": "#E74C3C"},
    {"name": "Tumor Suppressors", "genes": ["TP53", "PTEN", "RB1"], "color": "#3498DB"},  
    {"name": "DNA Repair", "genes": ["BRCA1", "BRCA2", "ATM"], "color": "#27AE60"},
    {"name": "Cell Cycle", "genes": ["CDK1", "CDK2", "CCND1"], "color": "#F39C12"}
]

for group in groups:
    add_search_group(
        payload=payload,
        group_name=group["name"],
        gene_names=group["genes"], 
        color=group["color"],
        curtain_type="TP"
    )
```

### Search Group Colors

CurtainUtils automatically assigns colors if not specified, but you can use predefined color schemes:

```python
# Distinct color palette for multiple groups
group_colors = [
    "#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A",
    "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9"
]

pathways = {
    "MAPK Signaling": ["MAPK1", "MAPK3", "RAF1"],
    "PI3K/AKT": ["PIK3CA", "AKT1", "MTOR"], 
    "p53 Pathway": ["TP53", "MDM2", "CDKN1A"],
    "Cell Cycle": ["CDK1", "CDK2", "CCNB1"]
}

for i, (pathway, genes) in enumerate(pathways.items()):
    add_search_group(
        payload=payload,
        group_name=pathway,
        gene_names=genes,
        color=group_colors[i % len(group_colors)],
        curtain_type="TP"
    )
```

### Working with Raw Data

When working with custom datasets, you can reference raw data for more flexible searching:

```python
import pandas as pd

# Load your raw data
raw_df = pd.read_csv("proteomics_data.csv")

high_abundance = raw_df[raw_df['Intensity'] > raw_df['Intensity'].quantile(0.9)]
add_search_group(
    payload=payload,
    group_name="High Abundance",
    protein_ids=high_abundance['Primary.IDs'].tolist(),
    color="#E67E22",
    curtain_type="TP"
)

membrane_proteins = raw_df[raw_df['Subcellular_Location'].str.contains('membrane', na=False)]
add_search_group(
    payload=payload, 
    group_name="Membrane Proteins",
    protein_ids=membrane_proteins['Primary.IDs'].tolist(),
    color="#9B59B6",
    curtain_type="TP"
)
```

### Best Practices

1. **Meaningful Names**: Use descriptive names that will be clear in the visualization
2. **Color Coordination**: Choose colors that provide good contrast and are colorblind-friendly
3. **Group Size**: Keep groups reasonably sized (10-50 proteins) for clear visualization
4. **Biological Relevance**: Group proteins by pathway, function, or biological significance
5. **Multiple Comparisons**: For Curtain datasets with multiple comparisons, consider whether you want comparison-specific groups

### Integration with Visualizations

Search groups automatically appear in:
- **Volcano plots**: Highlighted points with group colors
- **Bar charts**: Grouped and colored bars
- **Profile plots**: Colored traces for group members
- **Data tables**: Searchable and filterable by group

The groups become interactive elements in the web interface, allowing users to:
- Toggle group visibility
- Filter data by group membership
- Export group-specific results
- Modify group colors and names

## Custom Sample-Condition Mapping

By default, CurtainUtils expects sample names to follow the `condition.replicate` format (e.g., "Control.1", "Treatment.2"). However, you can override this with custom sample-condition mappings for datasets with non-standard naming conventions.

### Basic Sample Mapping

```python
from curtainutils.client import configure_sample_conditions

# Define custom sample-to-condition mapping
sample_mapping = {
    'Sample_001': 'Control',
    'Sample_002': 'Control', 
    'Sample_003': 'Control',
    'Sample_004': 'Treatment_A',
    'Sample_005': 'Treatment_A',
    'Sample_006': 'Treatment_A',
    'Sample_007': 'Treatment_B',
    'Sample_008': 'Treatment_B',
    'Sample_009': 'Treatment_B'
}

# Configure custom colors for conditions
condition_colors = {
    'Control': '#808080',
    'Treatment_A': '#FF6B6B', 
    'Treatment_B': '#4ECDC4'
}

# Apply the mapping to your payload
payload = configure_sample_conditions(
    payload=payload,
    sample_condition_map=sample_mapping,
    condition_colors=condition_colors,
    condition_order=['Control', 'Treatment_A', 'Treatment_B']
)
```

### Automatic Pattern Detection

CurtainUtils can analyze your sample names and suggest condition mappings:

```python
from curtainutils.client import detect_sample_patterns, validate_sample_mapping
import pandas as pd

# Load your raw data
raw_df = pd.read_csv("proteomics_data.csv")

# Get sample columns (exclude non-sample columns)
sample_columns = [col for col in raw_df.columns 
                 if col not in ['Primary.IDs', 'Gene.Name', 'Description']]

# Detect potential condition patterns
suggested_mapping = detect_sample_patterns(raw_df, sample_columns)
print("Suggested mapping:", suggested_mapping)

# Validate the mapping
validation_results = validate_sample_mapping(raw_df, suggested_mapping)
if validation_results['warnings']:
    for warning in validation_results['warnings']:
        print(f"Warning: {warning}")

# Apply if satisfactory
payload = configure_sample_conditions(payload, suggested_mapping)
```

### Supported Naming Patterns

The pattern detection recognizes common sample naming conventions:

| Pattern | Example | Detected Condition |
|---------|---------|-------------------|
| Standard format | `Control.1`, `Treatment.2` | `Control`, `Treatment` |
| Underscore format | `Control_rep1`, `Treatment_A_rep2` | `Control`, `Treatment A` |
| Common prefixes | `Ctrl01`, `Trt02`, `Drug03` | `Control`, `Treatment`, `Drug` |
| Numeric grouping | `Sample1`, `Sample2`, `Sample3` | `Group_1` (groups by 3s) |
| Custom alphanumeric | `CondA1`, `CondA2`, `CondB1` | `CondA`, `CondB` |

### Complex Sample Organization

For datasets with multiple experimental factors:

```python
# Multi-factor experiment: Treatment × Time Point
sample_mapping = {
    # Control samples
    'Ctrl_4h_Rep1': 'Control_4h',
    'Ctrl_4h_Rep2': 'Control_4h',
    'Ctrl_4h_Rep3': 'Control_4h',
    'Ctrl_24h_Rep1': 'Control_24h', 
    'Ctrl_24h_Rep2': 'Control_24h',
    'Ctrl_24h_Rep3': 'Control_24h',
    
    # Treatment samples  
    'Drug_4h_Rep1': 'Treatment_4h',
    'Drug_4h_Rep2': 'Treatment_4h',
    'Drug_4h_Rep3': 'Treatment_4h',
    'Drug_24h_Rep1': 'Treatment_24h',
    'Drug_24h_Rep2': 'Treatment_24h', 
    'Drug_24h_Rep3': 'Treatment_24h'
}

# Define condition colors
condition_colors = {
    'Control_4h': '#E8F4F8',
    'Control_24h': '#4A90A4', 
    'Treatment_4h': '#FADBD8',
    'Treatment_24h': '#E74C3C'
}

# Set condition order for logical grouping
condition_order = ['Control_4h', 'Control_24h', 'Treatment_4h', 'Treatment_24h']

payload = configure_sample_conditions(
    payload=payload,
    sample_condition_map=sample_mapping,
    condition_colors=condition_colors,
    condition_order=condition_order
)
```

### Validation and Quality Control

The validation function checks for common issues:

```python
# Validate sample mapping against raw data
validation = validate_sample_mapping(raw_df, sample_mapping)

# Check validation results
print("Validation Results:")
print(f"Missing samples: {validation['missing_samples']}")
print(f"Extra samples: {validation['extra_samples']}")  
print(f"Unbalanced conditions: {validation['unbalanced_conditions']}")

# Address any warnings
if validation['warnings']:
    print("Warnings found:")
    for warning in validation['warnings']:
        print(f"  - {warning}")
```

### Integration with Search Groups

Custom sample mappings work seamlessly with search group functionality:

```python
# First configure sample conditions
payload = configure_sample_conditions(payload, sample_mapping, condition_colors)

add_search_group(
    payload=payload,
    group_name="Key Proteins", 
    gene_names=["TP53", "MYC", "KRAS"],
    color="#FF6B6B",
    curtain_type="TP"
)
```

### Parameters Reference

#### `configure_sample_conditions()`

<table>
<tr><th>Parameter</th><th>Type</th><th>Description</th><th>Required</th></tr>
<tr><td>payload</td><td>Dict</td><td>The payload to modify</td><td>Yes</td></tr>
<tr><td>sample_condition_map</td><td>Dict[str, str]</td><td>Sample name to condition mapping</td><td>Yes</td></tr>
<tr><td>condition_colors</td><td>Dict[str, str]</td><td>Condition to hex color mapping</td><td>No</td></tr>
<tr><td>condition_order</td><td>List[str]</td><td>Order of conditions in visualizations</td><td>No</td></tr>
</table>

#### `detect_sample_patterns()`

<table>
<tr><th>Parameter</th><th>Type</th><th>Description</th><th>Required</th></tr>
<tr><td>raw_df</td><td>pd.DataFrame</td><td>Raw data containing sample columns</td><td>Yes</td></tr>
<tr><td>sample_columns</td><td>List[str]</td><td>List of sample column names</td><td>Yes</td></tr>
</table>

