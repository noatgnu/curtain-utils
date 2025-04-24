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