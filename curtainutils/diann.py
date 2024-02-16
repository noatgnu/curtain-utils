import io

import click
import pandas as pd
from uniprotparser.betaparser import UniprotParser, UniprotSequence
from sequal.sequence import Sequence

from curtainutils.common import read_fasta


def lambda_function_for_diann_ptm_single_site(row: pd.Series, modified_seq_col: str, entry_col: str, fasta_df: pd.DataFrame, modification_of_interests: str) -> pd.Series:
    seq = Sequence(row[modified_seq_col])
    stripped_seq = seq.to_stripped_string()
    print(f"Processing {row[modified_seq_col]}, {row[entry_col]}")
    if row[entry_col] in fasta_df["Entry"].values:
        for i in seq:
            if len(i.mods) > 0:
                for mod in i.mods:
                    if mod.value == modification_of_interests:
                        row["Position.in.peptide"] = i.position
                        row["Residue"] = i.value
                        row["Variant"] = row["Protein.Group"]

                        matched_acc_row = fasta_df[fasta_df["Entry"].str.contains(row[entry_col])]
                        if len(matched_acc_row) > 0:
                            for i2, row2 in matched_acc_row.iterrows():
                                seq = row2["Sequence"]
                                peptide_seq = stripped_seq
                                try:
                                    peptide_position = seq.index(peptide_seq)
                                except ValueError:
                                    peptide_position = seq.replace("I", "L").index(
                                        peptide_seq.replace("I", "L"))
                                    row["Comment"] = "I replaced by L"
                                print(peptide_position)
                                if peptide_position >= -1:

                                    position_in_protein = i.position + peptide_position
                                    row["Position"] = position_in_protein
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
                                    if "Protein names" in row2:
                                        row["Protein.Name"] = row2["Protein names"]
                                    row["Sequence.window"] = sequence_window
                                    break
                        break
    return row

def process_diann_ptm(pr_file_path: str, report_file_path: str, output_file: str, modification_of_interests: str = "UniMod:21", columns: str = "accession,id,sequence,protein_name", fasta_file: str = None):
    df = pd.read_csv(pr_file_path, sep="\t")
    localization_score_col = "PTM.Site.Confidence"
    modified_seq_col = "Modified.Sequence"
    index_col = "Precursor.Id"
    pr_file_col = "File.Name"
    protein_group_col = "Protein.Group"
    df_meta = pd.read_csv(report_file_path, sep="\t")
    df = df[df[modified_seq_col].str.contains(modification_of_interests)]
    parser = UniprotParser(columns=columns, include_isoform=True)
    fasta_df = []
    df["parse_id"] = df[protein_group_col].apply(lambda x: str(UniprotSequence(x, parse_acc=True)) if UniprotSequence(x, parse_acc=True).accession else x)
    if fasta_file:
        fasta_df = read_fasta(fasta_file)
    else:
        for i in parser.parse(df["parse_id"].unique().tolist()):
            fasta_df.append(pd.read_csv(io.StringIO(i), sep="\t"))
        if len(fasta_df) == 1:
            fasta_df = fasta_df[0]
        else:
            fasta_df = pd.concat(fasta_df, ignore_index=True)

    df = df.apply(lambda x: lambda_function_for_diann_ptm_single_site(x, modified_seq_col, "parse_id", fasta_df,
                                                                      modification_of_interests), axis=1)
    df.set_index(index_col, inplace=True)
    for i, g in df_meta.groupby(index_col):
        if i in df.index:
            highest_score = g[localization_score_col].max()
            df.loc[i, localization_score_col] = highest_score
    df.reset_index(inplace=True)
    df.to_csv(output_file, sep="\t",
              index=False)

@click.command()
@click.option("--pr_file_path", "-p", help="Path to the PR file to be processed")
@click.option("--report_file_path", "-r", help="Path to the report file to be processed")
@click.option("--output_file", "-o", help="Path to the output file")
@click.option("--modification_of_interests", "-m", help="Modification of interests", default="UniMod:21")
@click.option("--columns", "-c", help="UniProt data columns to be included", default="accession,id,sequence,protein_name")
@click.option("--fasta_file", "-f", help="Path to the fasta file", default=None)
def main(pr_file_path: str, report_file_path: str, output_file: str, modification_of_interests: str, columns: str, fasta_file: str):
    process_diann_ptm(pr_file_path, report_file_path, output_file, modification_of_interests, columns=columns, fasta_file=fasta_file)


if __name__ == '__main__':
    df = pd.read_csv(r"Z:\ALESSI\Toan\CURTAIN_Rebuttal\TiO2_diaPASEF_Test\Reports.pr_matrix.tsv", sep="\t")
    localization_score_col = "PTM.Site.Confidence"
    modified_seq_col = "Modified.Sequence"
    index_col = "Precursor.Id"
    pr_file_col = "File.Name"
    protein_group_col = "Protein.Group"
    df_meta = pd.read_csv(r"Z:\ALESSI\Toan\CURTAIN_Rebuttal\TiO2_diaPASEF_Test\Reports.tsv", sep="\t")
    df = df[df[modified_seq_col].str.contains("UniMod:21")]
    parser = UniprotParser(columns="accession,id,sequence", include_isoform=True)
    fasta_df = []
    df["parse_id"] = df[protein_group_col].apply(lambda x: str(UniprotSequence(x, parse_acc=True)) if UniprotSequence(x, parse_acc=True).accession else x)

    for i in parser.parse(df["parse_id"].unique().tolist()):
        fasta_df.append(pd.read_csv(io.StringIO(i), sep="\t"))
    if len(fasta_df) == 1:
        fasta_df = fasta_df[0]
    else:
        fasta_df = pd.concat(fasta_df, ignore_index=True)

    df = df.apply(lambda x: lambda_function_for_diann_ptm_single_site(x, modified_seq_col, "parse_id", fasta_df, "UniMod:21"), axis=1)
    df.set_index(index_col, inplace=True)
    for i, g in df_meta.groupby(index_col):
        if i in df.index:
            highest_score = g[localization_score_col].max()
            df.loc[i, localization_score_col] = highest_score
    df.reset_index(inplace=True)
    df.to_csv(r"Z:\ALESSI\Toan\CURTAIN_Rebuttal\TiO2_diaPASEF_Test\Reports.pr_matrix_processed.tsv", sep="\t", index=False)
