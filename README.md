# CurtainUtils
A utility package for converting different MS output files into a format usable by Curtain (https://curtain.proteo.info) and CurtainPTM (https://curtainptm.proteo.info)

## Installation
The package can be installed using the following command.
`pip install curtainutils`

## Convert MSFragger PTM single site output to CurtainPTM input

```bash
msf-curtainptm -f <MSFragger PTM single site output file> -i <index column with site information> -o <output file> -p <peptide column> -a <fasta file> 
```