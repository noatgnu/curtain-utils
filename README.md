# CurtainUtils
A utility package for converting different MS output files into a format usable by Curtain (https://curtain.proteo.info) and CurtainPTM (https://curtainptm.proteo.info)

## Installation
The package can be installed using the following command.
`pip install curtainutils`

## Convert MSFragger PTM single site output to CurtainPTM input
This script should be used to convert a differential analysis file that contains the index column and peptide column. The index column should also be the original index column output by MS-Fragger that contains both the Accession ID as well as position of the PTM within the protein sequence.

Usage:
```bash
msf-curtainptm -f <MSFragger PTM single site output file> -i <index column with site information> -o <output file> -p <peptide column> -a <fasta file> 
```

## Convert DIA-NN PTM output to CurtainPTM input
This script should be used to convert a differential analysis file that contains the following column "Modified.Sequence", "Precursor.Id", "Protein.Group" from the pr report file by combining the file with the Report file which contains the column "PTM.Site.Confidence".

```bash
diann-curtainptm -p <differential analysis file> -r <report file> -o <output file> -m <modification_of_interests from the Modified.Sequence column>
```

## Convert Spectronaut output to Curtain input
This script should be used to convert a differential analysis file that contains the "PTM_collapse_key" and "PEP.StrippedSequence" columns from the original Spectronaut output.

```bash
spn-curtainptm -f <differential analysis file> -o <output file>
```

## Submit data to a Curtain server
Usage example:

```python
from curtainutils.client import CurtainClient

de_file = r"differential-file-path"
raw_file = r"raw-file-path"

fc_col = "foldchange-column-name"
transform_fc = False
transform_significant = False
reverse_fc = False
p_col = "significance-column-name"

comp_col = ""  # Leave empty if no comparison column is used
comp_select = []  # Leave empty if no comparison column is used

primary_id_de_col = "primary-id-column-name-in-differential-file"
primary_id_raw_col = "primary-id-column-name-in-raw-file"

sample_cols = ["4Hr-AGB1.01", "4Hr-AGB1.02", "4Hr-AGB1.03", "4Hr-AGB1.04", "4Hr-AGB1.05", "24Hr-AGB1.01",
               "24Hr-AGB1.02", "24Hr-AGB1.03", "24Hr-AGB1.04", "24Hr-AGB1.05", "4Hr-Cis.01", "4Hr-Cis.02", "4Hr-Cis.03",
               "24Hr-Cis.01", "24Hr-Cis.02", "24Hr-Cis.03"]
c = CurtainClient("curtain-backend-url")
payload = c.create_curtain_session_payload(
    de_file,
    raw_file,
    fc_col,
    transform_fc,
    transform_significant,
    reverse_fc,
    p_col,
    comp_col,
    comp_select,
    primary_id_de_col,
    primary_id_raw_col,
    sample_cols
)

package = {
    "enable": "True",
    "description": payload["settings"]["description"],
    "curtain_type": "TP",
}

result = c.post_curtain_session(package, payload)
print(result)
```

## Submit data to a CurtainPTM server
```python
from curtainutils.client import CurtainClient

de_file = r"differential-file-path"
raw_file = r"raw-file-path"

fc_col = "foldchange-column-name"
transform_fc = False
transform_significant = False
reverse_fc = False
p_col = "significance-column-name"
comp_col = ""  # Leave empty if no comparison column is used
comp_select = []  # Leave empty if no comparison column is used
primary_id_de_col = "primary-id-column-name-in-differential-file"
primary_id_raw_col = "primary-id-column-name-in-raw-file"
sample_cols = []
peptide_col = "peptide-sequence-column-name"
acc_col = "protein-accession-column-name"
position_col = "position-in-protein-column-name"
position_in_peptide_col = "position-in-peptide-column-name"
sequence_window_col = "sequence-window-column-name"
score_col = "score-column-name"

c = CurtainClient("curtain-backend-url")

payload = c.create_curtain_ptm_session_payload(
    de_file,
    raw_file,
    fc_col,
    transform_fc,
    transform_significant,
    reverse_fc,
    p_col,
    comp_col,
    comp_select,
    primary_id_de_col,
    primary_id_raw_col,
    sample_cols,
    peptide_col,
    acc_col,
    position_col,
    position_in_peptide_col,
    sequence_window_col,
    score_col
)

package = {
    "enable": "True",
    "description": payload["settings"]["description"],
    "curtain_type": "PTM",
}

result = c.post_curtain_session(package, payload)
```