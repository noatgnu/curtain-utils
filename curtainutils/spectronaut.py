import io

import click
import pandas as pd
from uniprotparser.betaparser import UniprotSequence, UniprotParser
import re
from curtainutils.common import read_fasta
reg_pattern = re.compile("_\w(\d+)_")
protein_name_pattern = re.compile("(\w+_\w+)")
def lambda_function_for_spectronaut_ptm(row: pd.Series, index_col: str, peptide_col: str, fasta_df: pd.DataFrame) -> pd.Series:
    d = row[index_col].split("_")
    row["Position"] = int(d[-2][1:])
    if row["UniprotID"] in fasta_df["Entry"].values:
        matched_acc_row = fasta_df[fasta_df["Entry"].str.contains(row["UniprotID"])]
        reformat_seq = row[peptide_col].split(";")[0].upper()
        if len(matched_acc_row) > 0:
            for i2, row2 in matched_acc_row.iterrows():
                row2["PeptideSequence"] = reformat_seq[:len(reformat_seq)-2]
                seq = row2["Sequence"]
                try:
                    peptide_position = seq.index(row2["PeptideSequence"])
                except ValueError:
                    peptide_position = seq.replace("I", "L").index(
                        row2["PeptideSequence"].replace("I", "L"))
                    row["Comment"] = "I replaced by L"
                if peptide_position >= -1:
                    if "Protein names" in row2:
                        row["Protein.Name"] = row2["Protein names"]
                    position_in_peptide = row["Position"] - peptide_position
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
                    break
    return row

@click.command()
@click.option("--file_path", "-f", help="Path to the file to be processed")
@click.option("--index_col", "-i", help="Name of the index column", default="PTM_collapse_key")
@click.option("--peptide_col", "-p", help="Name of the peptide column", default="PEP.StrippedSequence")
@click.option("--output_file", "-o", help="Path to the output file")
@click.option("--fasta_file", "-a", help="Path to the fasta file")
@click.option("--columns", "-c", help="UniProt data columns to be included", default="accession,id,sequence,protein_name")
def process_spectronaut_ptm(
        file_path: str,
        index_col: str,
        peptide_col: str,
        output_file: str,
        fasta_file: str = "", columns: str = "accession,id,sequence,protein_name"):
    df = pd.read_csv(file_path, sep="\t")
    df["UniprotID"] = df[index_col].apply(lambda x: str(UniprotSequence(x, parse_acc=True)) if UniprotSequence(x, parse_acc=True).accession else x)
    if fasta_file:
        fasta_df = read_fasta(fasta_file)
    else:
        parser = UniprotParser(columns=columns, include_isoform=True)
        fasta_df = []
        for i in parser.parse(df["UniprotID"].unique().tolist()):
            fasta_df.append(pd.read_csv(io.StringIO(i), sep="\t"))
        if len(fasta_df) == 1:
            fasta_df = fasta_df[0]
        else:
            fasta_df = pd.concat(fasta_df, ignore_index=True)
    # peptide_col="PEP.StrippedSequence"
    # index_col="PTM_collapse_key"
    df = df.apply(lambda x: lambda_function_for_spectronaut_ptm(x, index_col, peptide_col, fasta_df), axis=1)
    df.to_csv(output_file, sep="\t", index=False)

@click.command()
@click.option("--file_path", "-f", help="Path to the file to be processed")
@click.option("--index_col", "-i", help="Name of the index column", default="PTM_collapse_key")
@click.option("--peptide_col", "-p", help="Name of the peptide column", default="PEP.StrippedSequence")
@click.option("--output_file", "-o", help="Path to the output file")
@click.option("--fasta_file", "-a", help="Path to the fasta file")
def main(file_path: str, index_col: str, peptide_col: str, output_file: str, fasta_file: str):
    process_spectronaut_ptm(file_path, index_col, peptide_col, output_file, fasta_file)