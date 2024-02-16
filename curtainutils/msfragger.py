import io

import click
import pandas as pd
from uniprotparser.betaparser import UniprotSequence, UniprotParser

from curtainutils.common import read_fasta
import re

reg_positon_residue = re.compile("_(\w)(\d+)")

def lambda_function_for_msfragger_ptm_single_site(row: pd.Series, index_col: str, peptide_col: str, fasta_df: pd.DataFrame, parse_from_peptide_col: bool = False) -> pd.Series:
    match = reg_positon_residue.search(row[index_col])
    print(f"Processing {row[index_col]}")
    if match:
        position = int(match.group(2))
        row["Position"] = position
        row["Residue"] = match.group(1)
        row["Peptide_Sequence"] = row[peptide_col].split(";")[0]
        if row["PrimaryID"] in fasta_df["Entry"].values:
            matched_acc_row = fasta_df[fasta_df["Entry"].str.contains(row["PrimaryID"])]

            if len(matched_acc_row) > 0:
                if not parse_from_peptide_col:
                    for i2, row2 in matched_acc_row.iterrows():
                        seq = row2["Sequence"]
                        if pd.notnull(seq):
                            peptide_seq = row["Peptide_Sequence"].split(";")[0].upper()
                            try:
                                peptide_position = seq.index(peptide_seq)
                            except ValueError:
                                peptide_position = seq.replace("I", "L").index(
                                    peptide_seq.replace("I", "L"))
                                row["Comment"] = "I replaced by L"
                            if peptide_position >= -1:

                                position_in_peptide = position - peptide_position
                                row["Position.in.peptide"] = position_in_peptide
                                row["Variant"] = row2["Entry"]
                                sequence_window = ""
                                if row["Position"] - 1 - 10 >= 0:
                                    sequence_window += seq[row["Position"] - 1 - 10:row["Position"] - 1]
                                else:
                                    sequence_window += seq[:row["Position"] - 1]
                                    if len(sequence_window) < 10:
                                        sequence_window = "_" * (10 - len(sequence_window)) + sequence_window
                                sequence_window += seq[row["Position"] - 1]
                                if row["Position"] + 10 <= len(seq):
                                    sequence_window += seq[row["Position"]:row["Position"] + 10]
                                else:
                                    sequence_window += seq[row["Position"]:]
                                    if len(sequence_window) < 21:
                                        sequence_window += "_" * (21 - len(sequence_window))

                                row["Sequence.window"] = sequence_window
                                if "Protein names" in row2:
                                    row["Protein.Name"] = row2["Protein names"]
                                break
                else:
                    # get index of character in lower case
                    for i, aa in enumerate(row["Peptide_Sequence"]):
                        if aa.islower() and row["Residue"] == aa.upper():
                            row["Position.in.peptide"] = i + 1
                            break

    return row


def process_msfragger_ptm_single_site(file_path: str, index_col: str, peptide_col: str, output_file: str, fasta_file: str = "", get_position_from_peptide_column: bool = False, columns: str = "accession,id,sequence,protein_name"):
    df = pd.read_csv(file_path, sep="\t")
    print(df.shape)
    df["PrimaryID"] = df[index_col].apply(lambda x: str(UniprotSequence(x, parse_acc=True)) if UniprotSequence(x, parse_acc=True).accession else x)
    if fasta_file:
        fasta_df = read_fasta(fasta_file)
    else:
        parser = UniprotParser(columns=columns, include_isoform=True)
        fasta_df = []
        for i in parser.parse(df["PrimaryID"].unique().tolist(), 500):
            fasta_df.append(pd.read_csv(io.StringIO(i), sep="\t"))
        if len(fasta_df) == 1:
            fasta_df = fasta_df[0]
        else:
            fasta_df = pd.concat(fasta_df, ignore_index=True)
    print(df.shape)
    df = df.apply(lambda x: lambda_function_for_msfragger_ptm_single_site(x, index_col, peptide_col, fasta_df, get_position_from_peptide_column),
                  axis=1)
    print(df.shape)

    df.to_csv(output_file, sep="\t", index=False)


@click.command()
@click.option("--file_path", "-f", help="Path to the file to be processed")
@click.option("--index_col", "-i", help="Name of the index column", default="Index")
@click.option("--peptide_col", "-p", help="Name of the peptide column", default="Peptide")
@click.option("--output_file", "-o", help="Path to the output file")
@click.option("--fasta_file", "-a", help="Path to the fasta file")
@click.option("--columns", "-c", help="UniProt data columns to be included", default="accession,id,sequence,protein_name")
def main(file_path: str, index_col: str, peptide_col: str, output_file: str, fasta_file: str, columns: str):
    process_msfragger_ptm_single_site(file_path, index_col, peptide_col, output_file, fasta_file, columns=columns)
