import io

import click
import pandas as pd
from uniprotparser.betaparser import UniprotSequence, UniprotParser
import re
from curtainutils.common import read_fasta
reg_pattern = re.compile("_\w(\d+)_")
protein_name_pattern = re.compile("(\w+_\w+)")
# def lambda_function_for_spectronaut_ptm(row: pd.Series, index_col: str, peptide_col: str, fasta_df: pd.DataFrame) -> pd.Series:
#     d = row[index_col].split("_")
#     row["Position"] = int(d[-2][1:])
#     if row["UniprotID"] in fasta_df["Entry"].values:
#         matched_acc_row = fasta_df[fasta_df["Entry"].str.contains(row["UniprotID"])]
#         reformat_seq = row[peptide_col].split(";")[0].upper()
#         if len(matched_acc_row) > 0:
#             for i2, row2 in matched_acc_row.iterrows():
#                 row2["PeptideSequence"] = reformat_seq[:len(reformat_seq)-2]
#                 seq = row2["Sequence"]
#                 try:
#                     peptide_position = seq.index(row2["PeptideSequence"])
#                 except ValueError:
#                     peptide_position = seq.replace("I", "L").index(
#                         row2["PeptideSequence"].replace("I", "L"))
#                     row["Comment"] = "I replaced by L"
#                 if peptide_position >= -1:
#                     if "Protein names" in row2:
#                         row["Protein.Name"] = row2["Protein names"]
#                     position_in_peptide = row["Position"] - peptide_position
#                     row["Position.in.peptide"] = position_in_peptide
#                     row["Variant"] = row2["Entry"]
#                     sequence_window = ""
#                     if row["Position"] - 1 - 10 >= 0:
#                         sequence_window += seq[row["Position"] - 1 - 10:row["Position"] - 1]
#                     else:
#                         sequence_window += seq[:row["Position"] - 1]
#                         if len(sequence_window) < 10:
#                             sequence_window = "_" * (10 - len(sequence_window)) + sequence_window
#                     sequence_window += seq[row["Position"] - 1]
#                     if row["Position"] + 10 <= len(seq):
#                         sequence_window += seq[row["Position"]:row["Position"] + 10]
#                     else:
#                         sequence_window += seq[row["Position"]:]
#                         if len(sequence_window) < 21:
#                             sequence_window += "_" * (21 - len(sequence_window))
#
#                     row["Sequence.window"] = sequence_window
#                     break
#     return row
def lambda_function_for_spectronaut_ptm(row: pd.Series, index_col: str, peptide_col: str, fasta_df: pd.DataFrame) -> pd.Series:
    """
    Process a row of Spectronaut PTM data to extract and calculate various fields.

    Args:
        row (pd.Series): A row from the Spectronaut PTM DataFrame.
        index_col (str): The name of the index column in the DataFrame.
        peptide_col (str): The name of the peptide column in the DataFrame.
        fasta_df (pd.DataFrame): A DataFrame containing FASTA sequences.

    Returns:
        pd.Series: The processed row with additional fields.
    """
    # Extract position from the index column
    d = row[index_col].split("_")
    row["Position"] = int(d[-2][1:])

    # Check if the UniprotID exists in the FASTA DataFrame
    uniprot_id = row["UniprotID"]
    if uniprot_id in fasta_df["Entry"].values:
        matched_acc_row = fasta_df[fasta_df["Entry"].str.contains(uniprot_id)]
        reformat_seq = row[peptide_col].split(";")[0].upper()

        if not matched_acc_row.empty:
            for _, row2 in matched_acc_row.iterrows():
                peptide_seq = reformat_seq[:len(reformat_seq)-2]
                seq = row2["Sequence"]

                # Find the position of the peptide sequence in the protein sequence
                try:
                    peptide_position = seq.index(peptide_seq)
                except ValueError:
                    try:
                        peptide_position = seq.replace("I", "L").index(peptide_seq.replace("I", "L"))
                        row["Comment"] = "I replaced by L"
                    except ValueError:
                        print(uniprot_id, peptide_seq)
                        variants = fasta_df[fasta_df["Entry"].str.contains(uniprot_id)]
                        print(variants)
                        for _, variant in variants.iterrows():
                            if "Sequence" in variant:
                                seq = variant["Sequence"]
                                try:
                                    peptide_position = seq.index(peptide_seq)
                                except ValueError:
                                    try:
                                        peptide_position = seq.replace("I", "L").index(peptide_seq.replace("I", "L"))
                                        row["Comment"] = "I replaced by L"
                                    except ValueError:
                                        continue
                                if peptide_position >= 0:
                                    break

                if peptide_position >= 0:
                    # Populate additional fields in the row
                    row["Protein.Name"] = row2.get("Protein names", "")
                    position_in_peptide = row["Position"] - peptide_position
                    row["Position.in.peptide"] = position_in_peptide
                    row["Variant"] = row2["Entry"]

                    # Calculate the sequence window
                    start = max(0, row["Position"] - 11)
                    end = min(len(seq), row["Position"] + 10)
                    sequence_window = seq[start:row["Position"] - 1] + seq[row["Position"] - 1] + seq[row["Position"]:end]

                    # Pad the sequence window if necessary
                    if start == 0:
                        sequence_window = "_" * (10 - (row["Position"] - 1)) + sequence_window
                    if end == len(seq):
                        sequence_window += "_" * (21 - len(sequence_window))

                    row["Sequence.window"] = sequence_window
                    break
    return row

def process_spectronaut_ptm(
        file_path: str,
        index_col: str,
        peptide_col: str,
        output_file: str,
        fasta_file: str = "", columns: str = "accession,id,sequence,protein_name"):
    """
    Process a Spectronaut PTM file to extract and calculate various fields, and save the processed data to an output file.

    Args:
        file_path (str): Path to the Spectronaut PTM file to be processed.
        index_col (str): Name of the index column in the DataFrame.
        peptide_col (str): Name of the peptide column in the DataFrame.
        output_file (str): Path to the output file where processed data will be saved.
        fasta_file (str, optional): Path to the FASTA file. If not provided, UniProt data will be fetched.
        columns (str, optional): UniProt data columns to be included. Defaults to "accession,id,sequence,protein_name".

    Returns:
        None
    """
    # Read the input file into a DataFrame
    df = pd.read_csv(file_path, sep="\t")

    # Extract UniprotID from the index column
    df["UniprotID"] = df[index_col].apply(lambda x: str(UniprotSequence(x, parse_acc=True)) if UniprotSequence(x, parse_acc=True).accession else x)

    # Read or fetch the FASTA data
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

    # Apply the lambda function to process each row
    df = df.apply(lambda x: lambda_function_for_spectronaut_ptm(x, index_col, peptide_col, fasta_df), axis=1)

    # Save the processed DataFrame to the output file
    df.to_csv(output_file, sep="\t", index=False)

@click.command()
@click.option("--file_path", "-f", help="Path to the file to be processed")
@click.option("--index_col", "-i", help="Name of the index column", default="PTM_collapse_key")
@click.option("--peptide_col", "-p", help="Name of the peptide column", default="PEP.StrippedSequence")
@click.option("--output_file", "-o", help="Path to the output file")
@click.option("--fasta_file", "-a", help="Path to the fasta file")
def main(file_path: str, index_col: str, peptide_col: str, output_file: str, fasta_file: str):
    process_spectronaut_ptm(file_path, index_col, peptide_col, output_file, fasta_file)