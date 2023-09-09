import pandas as pd
from uniprotparser.betaparser import UniprotSequence


def read_fasta(fasta_file: str) -> pd.DataFrame:
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        current_acc = ""
        for line in f:
            if line.startswith('>'):
                acc = UniprotSequence(line.strip(), True)

                if acc.accession:
                    fasta_dict[str(acc)] = ""
                    current_acc = str(acc)
                else:
                    fasta_dict[line.strip().replace(">", "")] = ""
                    current_acc = line.strip().replace(">", "")

            else:
                fasta_dict[current_acc] += line.strip()
    return pd.DataFrame([[k, fasta_dict[k]] for k in fasta_dict], columns=["Entry", "Sequence"])
