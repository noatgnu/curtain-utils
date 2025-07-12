import io

import click
import pandas as pd
from sequal.sequence import Sequence
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
def lambda_function_for_spectronaut_ptm(
    row: pd.Series, index_col: str, peptide_col: str, fasta_df: pd.DataFrame
) -> pd.Series:
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
        reformat_seq = row[peptide_col].split(";")[0].upper().replace("_", "")
        if not matched_acc_row.empty:
            for _, row2 in matched_acc_row.iterrows():
                peptide_seq = reformat_seq[: len(reformat_seq) - 2]
                peptide_seq = Sequence(peptide_seq).to_stripped_string()
                seq = row2["Sequence"]
                if pd.isnull(seq):
                    continue
                # Find the position of the peptide sequence in the protein sequence
                try:
                    peptide_position = seq.index(peptide_seq)
                except ValueError:
                    try:
                        peptide_position = seq.replace("I", "L").index(
                            peptide_seq.replace("I", "L")
                        )
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
                                        peptide_position = seq.replace("I", "L").index(
                                            peptide_seq.replace("I", "L")
                                        )
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
                    row["PeptideSequence"] = peptide_seq

                    # Calculate the sequence window
                    start = max(0, row["Position"] - 11)
                    end = min(len(seq), row["Position"] + 10)
                    sequence_window = (
                        seq[start : row["Position"] - 1]
                        + seq[row["Position"] - 1]
                        + seq[row["Position"] : end]
                    )

                    # Pad the sequence window if necessary
                    if start == 0:
                        sequence_window = (
                            "_" * (10 - (row["Position"] - 1)) + sequence_window
                        )
                    if end == len(seq):
                        sequence_window += "_" * (21 - len(sequence_window))

                    row["Sequence.window"] = sequence_window
                    break
    return row


def lambda_function_for_spectronaut_ptm_mode_2(
    row: pd.Series,
    index_col: str,
    peptide_col: str,
    fasta_df: pd.DataFrame,
    modification: str,
) -> pd.Series:
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
    if row[peptide_col].startswith("(") or row[peptide_col].startswith("["):
        row[peptide_col] = "_" + row[peptide_col].split(".")[0]
        seq = Sequence(row[peptide_col])
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
        seq = Sequence(row[peptide_col].split(".")[0])
        stripped_seq = seq.to_stripped_string()
    # Check if the UniprotID exists in the FASTA DataFrame
    entry = row["UniprotID"]
    if entry in fasta_df["Entry"].values:
        matched_acc_row = fasta_df[fasta_df["Entry"].str.contains(entry)]

        if not matched_acc_row.empty:
            for i in seq:
                if any(mod.value == modification for mod in i.mods):
                    row["Position.in.peptide"] = i.position
                    row["Residue"] = i.value
                    row["Variant"] = row["PG.ProteinGroups"]

                    for _, row2 in matched_acc_row.iterrows():
                        protein_seq = row2["Sequence"]
                        peptide_seq = stripped_seq
                        try:
                            peptide_position = protein_seq.index(peptide_seq)
                        except ValueError:
                            try:
                                peptide_position = protein_seq.replace("I", "L").index(
                                    peptide_seq.replace("I", "L")
                                )
                                row["Comment"] = "I replaced by L"
                            except ValueError:
                                print("Error", entry, peptide_seq)

                                continue
                                # for _, variant in variants.iterrows():
                                #    if "Sequence" in variant:
                                #        seq = variant["Sequence"]
                                #        try:
                                #            peptide_position = seq.index(peptide_seq)
                                #        except ValueError:
                                #            try:
                                #                peptide_position = seq.replace("I", "L").index(
                                #                    peptide_seq.replace("I", "L"))
                                #                row["Comment"] = "I replaced by L"
                                #            except ValueError:
                                #                continue
                                #        if peptide_position >= 0:
                                #            break
                        if peptide_position >= 0:
                            position_in_protein = i.position + peptide_position
                            row["Position"] = position_in_protein
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
                    break
    return row


def process_spectronaut_ptm(
    file_path: str,
    index_col: str,
    peptide_col: str,
    output_file: str,
    fasta_file: str = "",
    mode: str = "1",
    modification: str = "Phospho (STY)",
    columns: str = "accession,id,sequence,protein_name",
):
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
    df["UniprotID"] = df["Uniprot"].apply(
        lambda x: (
            str(UniprotSequence(x, parse_acc=True))
            if UniprotSequence(x, parse_acc=True).accession
            else x
        )
    )

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

    if mode == "1":
        # Apply the lambda function to process each row
        df = df.apply(
            lambda x: lambda_function_for_spectronaut_ptm(
                x, index_col, peptide_col, fasta_df
            ),
            axis=1,
        )
    elif mode == "2":
        ptm_group_cols = [i for i in df.columns if i.endswith("PTM.Group")]
        site_prob_cols = [i for i in df.columns if i.endswith("PTM.SiteProbability")]
        new_df = []
        for i, row in df.iterrows():
            if row["PTM.ModificationTitle"] == modification:
                new_row = {
                    "UniprotID": row["UniprotID"],
                    index_col: row[index_col],
                    "PEP.SiteProbability": "",
                    "PTM.ModificationTitle": row["PTM.ModificationTitle"],
                    "PG.ProteinGroups": row["PG.ProteinGroups"],
                    "PG.Genes": row["PG.Genes"],
                }
                for col in ptm_group_cols:
                    if row[col] != "Filtered" and pd.notnull(row[col]):
                        new_row[peptide_col] = row[col].replace("_", "").split(".")[0]
                        break
                site_probs = []
                for col in site_prob_cols:
                    if row[col] != "Filtered" and row[col] != "":
                        site_probs.append(float(row[col]))
                if len(site_probs) > 0:
                    try:
                        new_row["PEP.SiteProbability"] = max(site_probs)
                    except ValueError or TypeError:
                        new_row["PEP.SiteProbability"] = 0
                    new_df.append(new_row)
        if len(new_df) > 0:
            df = pd.DataFrame(new_df)
            df = df.apply(
                lambda x: lambda_function_for_spectronaut_ptm_mode_2(
                    x, index_col, peptide_col, fasta_df, modification
                ),
                axis=1,
            )
    # Save the processed DataFrame to the output file
    df.to_csv(output_file, sep="\t", index=False)


@click.command()
@click.option("--file_path", "-f", help="Path to the file to be processed")
@click.option(
    "--index_col", "-i", help="Name of the index column", default="PTM_collapse_key"
)
@click.option(
    "--peptide_col",
    "-p",
    help="Name of the peptide column",
    default="PEP.StrippedSequence",
)
@click.option("--output_file", "-o", help="Path to the output file")
@click.option("--fasta_file", "-a", help="Path to the fasta file")
@click.option("--mode", "-m", help="Mode of operation", default="1")
@click.option(
    "--modification", "-d", help="Modification to be processed", default="Phospho (STY)"
)
def main(
    file_path: str,
    index_col: str,
    peptide_col: str,
    output_file: str,
    fasta_file: str,
    mode: str,
    modification: str,
):
    process_spectronaut_ptm(
        file_path, index_col, peptide_col, output_file, fasta_file, mode, modification
    )
