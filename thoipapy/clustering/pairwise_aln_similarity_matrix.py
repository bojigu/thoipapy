import sys

import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
from Bio import SeqIO
import numpy as np
from pathlib import Path


def create_identity_matrix_from_protein_set(s, logging):
    set_number = s["set_number"]
    input_fasta = Path(f"/home/mark/Dropbox/tm_homodimer_dropbox/sets/fasta/set{set_number:02d}_full.fas")
    output_align = Path(f"/home/mark/Dropbox/tm_homodimer_dropbox/sets/fasta/{set_number:02d}_sim_matrix_alignments.txt")
    sim_matrix_xlsx = Path(f"/home/mark/Dropbox/tm_homodimer_dropbox/sets/fasta/{set_number:02d}_sim_matrix.xlsx")
    settings_csv = Path(f"/home/mark/Dropbox/tm_homodimer_dropbox/sets/fasta/{set_number:02d}_full_sim_matrix_settings.csv")
    gap_open = -15
    gap_extend = -0.5
    create_identity_matrix_using_pairwise_alignments(input_fasta, output_align, sim_matrix_xlsx, settings_csv, gap_open, gap_extend, logging)


def create_identity_matrix_using_pairwise_alignments(input_fasta, output_align, sim_matrix_xlsx, settings_csv, gap_open, gap_extend, logging, matrix = matlist.blosum62, aln_cutoff=15):
    logging.info('~~~~~~~~~~~~                 starting create_identity_matrix_using_pairwise_alignments              ~~~~~~~~~~~~')
    acc_pairs_analysed = []
    acc_pairs_above_cutoff = []
    df_out = pd.DataFrame()
    with open(output_align, "w") as f:
        for seq_record in SeqIO.parse(input_fasta, "fasta"):
            acc = seq_record.id
            query = seq_record.seq
            for target_record in SeqIO.parse(input_fasta, "fasta"):
                target_acc = target_record.id
                target_seq = target_record.seq
                if acc == target_acc:
                    df_out.at[acc, target_acc] = 0
                    continue
                if (acc, target_acc) not in acc_pairs_analysed and (target_acc, acc) not in acc_pairs_analysed:
                    alignments = pairwise2.align.globalds(query, target_seq, matrix, gap_open, gap_extend)
                    a = alignments[0]
                    q = np.array(list(a[0]))
                    m = np.array(list(a[1]))
                    b = q == m
                    same = b.sum()
                    ident = same / len(b) * 100
                    df_out.at[acc, target_acc] = ident
                    df_out.at[target_acc, acc] = ident
                    acc_pairs_analysed.append((acc, target_acc))

                    if ident > aln_cutoff:
                        acc_pairs_above_cutoff.append((acc, target_acc))
                        description = "acc={}, target_acc={}, ident={:0.2f}, n_alignments={}".format(acc, target_acc, ident, len(alignments))
                        alignment = format_alignment(*alignments[0])
                        result = "{}, alignment=\n{}\n".format(description, alignment)
                        f.write(result)
                        f.write("\n------------------------------------\n")
                        sys.stdout.write(".")
                        sys.stdout.flush()

    # save the similarity matrix
    writer = pd.ExcelWriter(sim_matrix_xlsx)
    df_out.to_excel(writer, sheet_name="sim_matrix")

    df_pairs_above_cutoff = pd.DataFrame()
    df_pairs_above_cutoff["pairs"] = acc_pairs_above_cutoff
    df_pairs_above_cutoff.to_excel(writer, sheet_name="pairs_above_cutoff")

    # save the settings used
    df_settings = pd.DataFrame()
    df_settings.at["input_fasta", "value"] = input_fasta
    df_settings.at["alignment_file", "value"] = output_align
    df_settings.at["sim_matrix_csv", "value"] = sim_matrix_xlsx
    df_settings.at["matrix", "value"] = "matlist.blosum62"
    df_settings.at["gap_open", "value"] = gap_open
    df_settings.at["gap_extend", "value"] = gap_extend
    df_settings.at["acc_pairs_analysed", "value"] = acc_pairs_analysed
    df_settings.to_excel(writer, sheet_name="settings")

    writer.save()
    writer.close()

    logging.info('~~~~~~~~~~~~                 finished create_identity_matrix_using_pairwise_alignments              ~~~~~~~~~~~~')