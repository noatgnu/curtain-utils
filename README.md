# CurtainUtils
A utility package for converting different MS output files into a format usable by Curtain (https://curtain.proteo.info) and CurtainPTM (https://curtainptm.proteo.info)

## Installation
The package can be installed using the following command.
`pip install curtainutils`

## Convert MSFragger PTM single site output to CurtainPTM input

Usage:
```bash
msf-curtainptm -f <MSFragger PTM single site output file> -i <index column with site information> -o <output file> -p <peptide column> -a <fasta file> 
```

## Submit data to a Curtain server

Usage example:
```python
import curtainutils.common as common

de_file = r"differential-file-path"
raw_file = r"raw-file-path"

fc_col = "foldchange-column-name"
transform_fc = False
transform_significant = False
reverse_fc = False
p_col = "significance-column-name"

comp_col = "" # Leave empty if no comparison column is used
comp_select = [] # Leave empty if no comparison column is used

primary_id_de_col = "primary-id-column-name-in-differential-file"
primary_id_raw_col = "primary-id-column-name-in-raw-file"

sample_cols = ["4Hr-AGB1.01", "4Hr-AGB1.02", "4Hr-AGB1.03", "4Hr-AGB1.04", "4Hr-AGB1.05", "24Hr-AGB1.01", "24Hr-AGB1.02", "24Hr-AGB1.03", "24Hr-AGB1.04", "24Hr-AGB1.05", "4Hr-Cis.01", "4Hr-Cis.02", "4Hr-Cis.03", "24Hr-Cis.01", "24Hr-Cis.02", "24Hr-Cis.03"]

payload = common.create_curtain_session_payload(
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

result = common.post_curtain_session(package, payload, "back-end-url")
print(result)
```