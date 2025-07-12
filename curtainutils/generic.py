import io

import click
import pandas as pd
from uniprotparser.betaparser import UniprotParser, UniprotSequence
from sequal.sequence import Sequence

from curtainutils.common import read_fasta

if __name__ == "__main__":
    fasta_file = r"C:\Users\Toan Phung\Downloads\20250113_Mouse_Reviewed.fasta"
    input_file = r"C:\Users\Toan Phung\Downloads\20250414_All-phosphoPeptide.txt"
    accession_col = "PrimaryID"
    fasta_df = read_fasta(fasta_file)
    df = pd.read_csv(input_file, sep="\t")
    fin_df = []
    for i, row in df.iterrows():
        accession = UniprotSequence(row[accession_col], parse_acc=True)
        acc = accession.accession
        if accession.isoform:
            acc += "-" + accession.isoform

        sequence = row[accession_col].split(";")[1]
        stripped_seq = ""
        # remove numbers and special characters
        sequence = sequence[:-1]
        if sequence.startswith("("):
            modded_seq = "_" + sequence
            seq = Sequence(modded_seq)
            seq = seq[1:]
            seq2 = ""
            for i in seq:
                if i.mods:
                    seq2 += i.value + "(" + i.mods[0].value + ")"
                else:
                    seq2 += i.value
            seq = Sequence(seq2)
            stripped_seq = seq.to_stripped_string()
        else:
            seq = Sequence(sequence)
            stripped_seq = seq.to_stripped_string()
        print(sequence)
        if acc in fasta_df["Entry"].values:
            matched_acc_row = fasta_df[fasta_df["Entry"].str.contains(acc)]
            if not matched_acc_row.empty:
                for i in seq:
                    if any(mod.value == "UniMod:21" for mod in i.mods):
                        row["Position.in.peptide"] = i.position + 1
                        row["Peptide"] = sequence
                        row["Variant"] = acc
                        for _, row2 in matched_acc_row.iterrows():
                            protein_seq = row2["Sequence"]
                            peptide_seq = stripped_seq
                            try:
                                peptide_position = protein_seq.index(peptide_seq)
                            except ValueError:
                                try:
                                    peptide_position = protein_seq.replace(
                                        "I", "L"
                                    ).index(peptide_seq.replace("I", "L"))
                                    row["Comment"] = "I replaced by L"
                                except ValueError:
                                    print("Error", acc, peptide_seq)

                                    continue
                            if peptide_position >= 0:
                                position_in_protein = i.position + peptide_position
                                row["Position"] = position_in_protein + 1
                                row["Variant"] = row2["Entry"]

                                start = max(0, position_in_protein - 11)
                                end = min(len(protein_seq), position_in_protein + 10)
                                sequence_window = (
                                    protein_seq[start : position_in_protein - 1]
                                    + protein_seq[position_in_protein - 1]
                                    + protein_seq[position_in_protein:end]
                                )

                                if start == 0:
                                    sequence_window = (
                                        "_" * (10 - (position_in_protein - 1))
                                        + sequence_window
                                    )
                                if end == len(protein_seq):
                                    sequence_window += "_" * (21 - len(sequence_window))

                                row["Sequence.window"] = sequence_window
                                break
                        fin_df.append(row)
                        break
    fin_df = pd.DataFrame(fin_df)
    fin_df.to_csv(input_file + "output.txt", sep="\t", index=False)
