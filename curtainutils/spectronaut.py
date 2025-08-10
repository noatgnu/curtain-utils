"""
Spectronaut PTM Data Processing Module.

This module provides functionality for processing post-translational modification (PTM) data
from Spectronaut, a mass spectrometry data analysis software. It handles mapping of peptide
sequences to protein sequences, position calculation, and sequence window generation around
modification sites.

The module supports three different processing modes:
- Mode 1: Basic PTM processing with position extraction from index columns
- Mode 2: PTM group processing with site probability calculations
- Mode 3: ProForma sequence processing with enhanced position validation

Key Features:
- UniProt ID parsing and FASTA sequence retrieval
- Peptide-to-protein sequence mapping with I/L ambiguity handling
- Modification site position calculation and validation
- Sequence window generation around modification sites
- Support for multiple PTM data formats and notations

Typical Usage:
    >>> from curtainutils.spectronaut import process_spectronaut_ptm
    >>> process_spectronaut_ptm(
    ...     "spectronaut_data.tsv",
    ...     "PTM_collapse_key",
    ...     "PEP.StrippedSequence", 
    ...     "processed_output.tsv"
    ... )

Command Line Usage:
    $ python -m curtainutils.spectronaut -f data.tsv -o output.tsv
"""

import io

import click
import pandas as pd
from sequal.sequence import Sequence
from uniprotparser.betaparser import UniprotSequence, UniprotParser
import re
from curtainutils.common import read_fasta

reg_pattern = re.compile("_\w(\d+)_")
"""Regex pattern for extracting position information from index columns."""

protein_name_pattern = re.compile("(\w+_\w+)")
"""Regex pattern for extracting protein names from identifiers."""


def _find_peptide_in_protein(peptide_seq: str, protein_seq: str) -> tuple[int, str]:
    """
    Find the position of a peptide sequence within a protein sequence.
    
    Handles I/L amino acid ambiguity by attempting substitution if direct match fails.
    
    Args:
        peptide_seq (str): The peptide sequence to find.
        protein_seq (str): The protein sequence to search in.
        
    Returns:
        tuple[int, str]: A tuple containing:
            - peptide_position (int): 0-based position of peptide in protein, -1 if not found
            - comment (str): Comment about any substitutions made (e.g., "I replaced by L")
    """
    try:
        return protein_seq.index(peptide_seq), ""
    except ValueError:
        try:
            position = protein_seq.replace("I", "L").index(peptide_seq.replace("I", "L"))
            return position, "I replaced by L"
        except ValueError:
            return -1, ""


def _generate_sequence_window(protein_seq: str, position: int, window_size: int = 21) -> str:
    """
    Generate a sequence window around a modification position.
    
    Creates a window of specified size centered on the modification site,
    padding with underscores when near protein termini.
    
    Args:
        protein_seq (str): The full protein sequence.
        position (int): 1-based position of the modification in the protein.
        window_size (int, optional): Size of the window to generate. Defaults to 21.
        
    Returns:
        str: A sequence window with the modification site in the center.
    """
    modification_pos_0based = position - 1
    half_window = window_size // 2
    start = max(0, modification_pos_0based - half_window)
    end = min(len(protein_seq), modification_pos_0based + half_window + 1)
    
    sequence_window = (
        protein_seq[start:modification_pos_0based]
        + protein_seq[modification_pos_0based]
        + protein_seq[modification_pos_0based + 1:end]
    )
    
    if start == 0:
        padding = half_window - modification_pos_0based
        sequence_window = "_" * padding + sequence_window
    
    if end == len(protein_seq):
        current_length = len(sequence_window)
        if current_length < window_size:
            sequence_window += "_" * (window_size - current_length)
    
    return sequence_window


def _process_peptide_sequence(peptide_str: str, mode: str) -> tuple:
    """
    Process peptide sequence based on the processing mode.
    
    Args:
        peptide_str (str): Raw peptide sequence string.
        mode (str): Processing mode ("1", "2", or "3").
        
    Returns:
        tuple: Contains (seq_object, stripped_sequence) where seq_object is the
               parsed Sequence object and stripped_sequence is the clean amino acid sequence.
    """
    if mode == "1":
        reformat_seq = peptide_str.split(";")[0].upper().replace("_", "")
        peptide_seq = reformat_seq[: len(reformat_seq) - 2]
        peptide_seq = Sequence(peptide_seq).to_stripped_string()
        return None, peptide_seq
    
    elif mode in ["2", "3"]:
        if peptide_str.startswith("(") or peptide_str.startswith("["):
            if mode == "3":
                seq = Sequence.from_proforma(peptide_str)
            else:
                peptide_str = "_" + peptide_str.split(".")[0]
                seq = Sequence(peptide_str)
            seq = seq[1:]
            seq2 = ""
            for i in seq:
                if i.mods:
                    seq2 += i.value + "(" + i.mods[0].value + ")"
                else:
                    seq2 += i.value
            seq = Sequence(seq2)
        else:
            seq = Sequence(peptide_str.split(".")[0] if mode == "2" else peptide_str)
        
        return seq, seq.to_stripped_string()
    
    return None, ""


def lambda_function_for_spectronaut_ptm(
    row: pd.Series, index_col: str, peptide_col: str, fasta_df: pd.DataFrame
) -> pd.Series:
    """
    Process a row of Spectronaut PTM data to map peptides to protein sequences and calculate position information.

    This function extracts modification position information from the index column, maps peptide sequences
    to their corresponding protein sequences in the FASTA database, and calculates various position-related
    fields including sequence windows around modification sites.

    Args:
        row (pd.Series): A row from the Spectronaut PTM DataFrame containing peptide and protein information.
        index_col (str): The name of the index column containing position information in the format "*_*X123_*".
        peptide_col (str): The name of the column containing the peptide sequence.
        fasta_df (pd.DataFrame): A DataFrame containing FASTA sequences with 'Entry' and 'Sequence' columns.

    Returns:
        pd.Series: The processed row with additional fields:
            - Position: Absolute position of modification in protein (1-based)
            - Position.in.peptide: Position of modification within the peptide
            - Protein.Name: Name of the protein from FASTA data
            - Variant: UniProt accession number
            - PeptideSequence: Clean peptide sequence without modifications
            - Sequence.window: 21-amino acid window centered on modification site
            - Comment: Notes about sequence matching (e.g., "I replaced by L")

    Note:
        - Handles I/L amino acid ambiguity by attempting substitution
        - Creates 21-residue sequence windows with padding ('_') when near protein termini
        - Processes semicolon-separated peptide sequences by taking the first entry
    """
    d = row[index_col].split("_")
    row["Position"] = int(d[-2][1:])

    uniprot_id = row["UniprotID"]
    if uniprot_id in fasta_df["Entry"].values:
        matched_acc_row = fasta_df[fasta_df["Entry"].str.contains(uniprot_id)]
        _, peptide_seq = _process_peptide_sequence(row[peptide_col], "1")
        
        if not matched_acc_row.empty:
            for _, row2 in matched_acc_row.iterrows():
                seq = row2["Sequence"]
                if pd.isnull(seq):
                    continue
                
                peptide_position, comment = _find_peptide_in_protein(peptide_seq, seq)
                if comment:
                    row["Comment"] = comment
                
                if peptide_position == -1:
                    variants = fasta_df[fasta_df["Entry"].str.contains(uniprot_id)]
                    for _, variant in variants.iterrows():
                        if "Sequence" in variant:
                            seq = variant["Sequence"]
                            peptide_position, comment = _find_peptide_in_protein(peptide_seq, seq)
                            if comment:
                                row["Comment"] = comment
                            if peptide_position >= 0:
                                break
                
                if peptide_position >= 0:
                    row["Protein.Name"] = row2.get("Protein names", "")
                    position_in_peptide = row["Position"] - (peptide_position + 1) + 1
                    row["Position.in.peptide"] = position_in_peptide
                    row["Variant"] = row2["Entry"]
                    row["PeptideSequence"] = peptide_seq
                    row["Sequence.window"] = _generate_sequence_window(seq, row["Position"])
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
    Process Spectronaut PTM data in mode 2, handling modification-specific peptide sequences.

    This function processes peptide sequences that may contain modification information in
    parentheses or brackets, identifies specific modification sites, and maps them to
    protein sequences to calculate absolute positions and sequence windows.

    Args:
        row (pd.Series): A row from the Spectronaut PTM DataFrame.
        index_col (str): The name of the index column containing position information.
        peptide_col (str): The name of the column containing the peptide sequence with modifications.
        fasta_df (pd.DataFrame): A DataFrame containing FASTA sequences with 'Entry' and 'Sequence' columns.
        modification (str): The specific modification type to search for (e.g., "Phospho (STY)").

    Returns:
        pd.Series: The processed row with additional fields:
            - Position: Absolute position of modification in protein (1-based)
            - Position.in.peptide: Position of modification within the peptide
            - Residue: The amino acid residue that is modified
            - Variant: UniProt accession number
            - Sequence.window: 21-amino acid window centered on modification site
            - Comment: Notes about sequence matching (e.g., "I replaced by L")

    Note:
        - Handles peptide sequences starting with '(' or '[' by preprocessing
        - Uses the sequal library to parse modification information
        - Only processes residues that have the specified modification type
    """
    d = row[index_col].split("_")
    row["Position"] = int(d[-2][1:])
    
    seq, stripped_seq = _process_peptide_sequence(row[peptide_col], "2")
    
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
                        peptide_position, comment = _find_peptide_in_protein(stripped_seq, protein_seq)
                        if comment:
                            row["Comment"] = comment
                                
                        if peptide_position >= 0:
                            position_in_protein = (peptide_position + 1) + (i.position - 1) + 1
                            row["Position"] = position_in_protein
                            row["Variant"] = row2["Entry"]
                            row["Sequence.window"] = _generate_sequence_window(protein_seq, position_in_protein)
                            break
                    break
    return row

def lambda_function_for_spectronaut_ptm_mode_3(
    row: pd.Series,
    index_col: str,
    peptide_col: str,
    fasta_df: pd.DataFrame,
    modification: str,
) -> pd.Series:
    """
    Process Spectronaut PTM data in mode 3, with enhanced ProForma sequence parsing and position validation.

    This function handles ProForma-formatted peptide sequences and performs strict position validation
    to ensure modification sites match between the index column and the peptide sequence. It supports
    both standard peptide sequences and ProForma notation with modifications.

    Args:
        row (pd.Series): A row from the Spectronaut PTM DataFrame.
        index_col (str): The name of the index column containing position information.
        peptide_col (str): The name of the column containing the peptide sequence (may be ProForma format).
        fasta_df (pd.DataFrame): A DataFrame containing FASTA sequences with 'Entry' and 'Sequence' columns.
        modification (str): The specific modification type to search for (e.g., "Phospho (STY)").

    Returns:
        pd.Series: The processed row with additional fields:
            - Position: Absolute position of modification in protein (1-based, validated)
            - Position.in.peptide: Position of modification within the peptide (1-based)
            - Residue: The amino acid residue that is modified
            - Variant: UniProt accession number
            - Sequence.window: 21-amino acid window centered on modification site
            - Protein.Name: Name of the protein from FASTA data
            - Comment: Notes about sequence matching (e.g., "I replaced by L")

    Note:
        - Supports ProForma sequence notation starting with '(' or '['
        - Performs position validation to ensure consistency between index and peptide
        - Returns early for null or empty peptide sequences
        - Includes debug output for the index column
    """
    d = row[index_col].split("_")
    row["Position"] = int(d[-2][1:])
    
    if pd.isnull(row[peptide_col]) or row[peptide_col].strip() == "":
        return row
    
    seq, stripped_seq = _process_peptide_sequence(row[peptide_col], "3")
    
    entry = row["UniprotID"]
    if entry in fasta_df["Entry"].values:
        matched_acc_row = fasta_df[fasta_df["Entry"].str.contains(entry)]
        if not matched_acc_row.empty:
            for i in seq:
                if any(mod.value == modification for mod in i.mods):
                    row["Position.in.peptide"] = str(i.position + 1)
                    row["Residue"] = i.value

                    for _, row2 in matched_acc_row.iterrows():
                        protein_seq = row2["Sequence"]
                        peptide_position, comment = _find_peptide_in_protein(stripped_seq, protein_seq)
                        if comment:
                            row["Comment"] = comment
                        
                        if peptide_position >= 0:
                            position_in_protein = i.position + peptide_position
                            if row["Position"] != position_in_protein + 1:
                                continue
                            row["Position"] = position_in_protein + 1
                            row["Variant"] = row2["Entry"]
                            row["Sequence.window"] = _generate_sequence_window(protein_seq, position_in_protein + 1)
                            row["Protein.Name"] = row2.get("Protein names", "")
                            break
                    if "Sequence.window" not in row:
                        continue
                    break
    return row

def process_spectronaut_ptm(
    file_path: str,
    index_col: str,
    peptide_col: str,
    output_file: str,
    fasta_file: str = "",
    uniprot_id_col: str = "UniprotID",
    mode: str = "1",
    modification: str = "Phospho (STY)",
    columns: str = "accession,id,sequence,protein_name",
):
    """
    Process a Spectronaut PTM file to map modifications to protein sequences and calculate position information.

    This function serves as the main entry point for processing Spectronaut PTM data files.
    It reads the input file, processes UniProt IDs, fetches or reads FASTA sequences,
    and applies different processing modes depending on the data format and requirements.

    Args:
        file_path (str): Path to the tab-separated Spectronaut PTM file to be processed.
        index_col (str): Name of the index column containing position information (e.g., "PTM_collapse_key").
        peptide_col (str): Name of the column containing peptide sequences (e.g., "PEP.StrippedSequence").
        output_file (str): Path where the processed tab-separated output file will be saved.
        fasta_file (str, optional): Path to a local FASTA file. If empty, UniProt data will be fetched online.
        uniprot_id_col (str, optional): Name of the column containing UniProt IDs. Defaults to "UniprotID".
        mode (str, optional): Processing mode - "1" for basic processing, "2" for PTM group processing,
            "3" for ProForma sequence processing. Defaults to "1".
        modification (str, optional): Target modification type for modes 2 and 3. Defaults to "Phospho (STY)".
        columns (str, optional): Comma-separated list of UniProt columns to retrieve.
            Defaults to "accession,id,sequence,protein_name".

    Returns:
        None: The function saves results directly to the output file.

    Raises:
        FileNotFoundError: If the input file cannot be found.
        ValueError: If the mode is not recognized or required columns are missing.

    Note:
        - Mode 1: Basic PTM processing with position extraction from index column
        - Mode 2: Advanced processing for PTM groups with site probability calculation
        - Mode 3: ProForma sequence processing with position validation
        - UniProt data is fetched online if no FASTA file is provided
        - All modes output tab-separated files with position and sequence window information

    Example:
        >>> process_spectronaut_ptm(
        ...     "data.tsv",
        ...     "PTM_collapse_key", 
        ...     "PEP.StrippedSequence",
        ...     "output.tsv",
        ...     mode="1"
        ... )
    """
    df = pd.read_csv(file_path, sep="\t")
    
    if uniprot_id_col:
        df["UniprotID"] = df[uniprot_id_col].apply(
            lambda x: (
                str(UniprotSequence(x, parse_acc=True))
                if UniprotSequence(x, parse_acc=True).accession
                else x
            )
        )
    else:
        df[index_col] = df["Uniprot"].apply(
            lambda x: (
                str(UniprotSequence(x, parse_acc=True))
                if UniprotSequence(x, parse_acc=True).accession
                else x.split("_")[0]
            )
        )

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
    elif mode == "3":
        for i, row in df.iterrows():
            df.at[i, "ProcessedSequence"] = row[peptide_col].replace("_", "").split(".")[0]
            df.at[i, "PG.ProteinGroups"] = row[uniprot_id_col]
        df = df.apply(
            lambda x: lambda_function_for_spectronaut_ptm_mode_3(
                x, index_col, "ProcessedSequence", fasta_df, modification
            ),
            axis=1,
        )
    
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
@click.option(
    "--uniprot_id_col",
    "-u",
    help="Column name for Uniprot ID",
    default="",
)
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
    uniprot_id_col: str,
    mode: str,
    modification: str,
):
    """
    Command-line interface for processing Spectronaut PTM data files.

    This CLI tool provides access to the process_spectronaut_ptm function with configurable
    parameters for different processing modes and data formats.

    Args:
        file_path (str): Path to the input Spectronaut PTM file (tab-separated).
        index_col (str): Name of the column containing position information.
        peptide_col (str): Name of the column containing peptide sequences.
        output_file (str): Path for the output processed file.
        fasta_file (str): Optional path to local FASTA file (fetches from UniProt if empty).
        uniprot_id_col (str): Name of the column containing UniProt IDs.
        mode (str): Processing mode ("1", "2", or "3").
        modification (str): Target modification type for modes 2 and 3.

    Example:
        $ python spectronaut.py -f data.tsv -i PTM_collapse_key -p PEP.StrippedSequence -o output.tsv -m 1
    """
    process_spectronaut_ptm(
        file_path, index_col, peptide_col, output_file, fasta_file, uniprot_id_col, mode, modification
    )
