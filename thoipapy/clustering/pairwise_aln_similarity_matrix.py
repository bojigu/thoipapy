import sys

import numpy as np
import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
from Bio import SeqIO
from ast import literal_eval as make_tuple
from pathlib import Path
from typing import List, Set, Union


def create_identity_matrix_from_protein_set(s, logging):
    setname = s["setname"]
    protein_set_full_seq_fasta = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/clusters/{setname}_full_seqs.fas"
    output_align = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/clusters/{setname}_sim_matrix_alignments.txt"
    sim_matrix_xlsx = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/clusters/{setname}_sim_matrix.xlsx"
    if not sim_matrix_xlsx.parent.is_dir():
        sim_matrix_xlsx.parent.mkdir(parents=True)
    gap_open = -40.0
    gap_extend = -0.5
    aln_cutoff = 15
    create_identity_matrix_using_pairwise_alignments(protein_set_full_seq_fasta, output_align, sim_matrix_xlsx, gap_open, gap_extend, logging, aln_cutoff=aln_cutoff)


def create_identity_matrix_using_pairwise_alignments(protein_set_full_seq_fasta: Union[Path, str], output_align: Union[Path, str], ident_matrix_xlsx: Union[Path, str], gap_open: float, gap_extend: float, logging, matrix = matlist.blosum62, aln_cutoff: float = 15.0):
    """Create identity matrix using pairwise alignments.

    This is a clustering method used to complement cd-hit, which is not designed for clustering at very low levels of identity.

    The method is slow, as it utilises the pairwise alignment method from biopython.

    The aln_cutoff determines the pairwise identity that admits the sequence into the cluster.
    Clusters comprise any proteins that could possibly be homologous.
    Not all proteins in a cluster may show identity in the pairwise alignments.
    For example, in cluster of three proteins [a,b,c], the pairwise alignments of a-b and b-c might
    have identity above the threshold, even if the alignment of a-c does not.
    """
    logging.info('~~~~~~~~~~~~                 starting create_identity_matrix_using_pairwise_alignments              ~~~~~~~~~~~~')
    acc_pairs_analysed = []
    acc_pairs_above_cutoff = []
    df_ident = pd.DataFrame()
    with open(output_align, "w") as f:
        for seq_record in SeqIO.parse(protein_set_full_seq_fasta, "fasta"):
            acc = seq_record.id
            query = seq_record.seq
            for target_record in SeqIO.parse(protein_set_full_seq_fasta, "fasta"):
                target_acc = target_record.id
                target_seq = target_record.seq
                if acc == target_acc:
                    df_ident.at[acc, target_acc] = 0
                    continue
                if (acc, target_acc) not in acc_pairs_analysed and (target_acc, acc) not in acc_pairs_analysed:
                    alignments = pairwise2.align.globalds(query, target_seq, matrix, gap_open, gap_extend)
                    a = alignments[0]
                    q = np.array(list(a[0]))
                    m = np.array(list(a[1]))
                    b = q == m
                    same = b.sum()
                    ident = same / len(b) * 100
                    df_ident.at[acc, target_acc] = ident
                    df_ident.at[target_acc, acc] = ident
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

    """ df_ident looks like this:
        ...              0-P20963-NMR  1-P21709-NMR  2-O43914-NMR  3-P05067-NMR  4-P09619-NMR  5-P22607-NMR  6-O15455-NMR
    0-P20963-NMR      0.000000      2.254098     14.024390      3.636364      2.260398      2.605459      2.876106
    1-P21709-NMR      2.254098      0.000000      2.459016      8.773679     15.926893     15.779468      9.705043
    2-O43914-NMR     14.024390      2.459016      0.000000      1.683938      1.808318      2.605459      2.212389
    3-P05067-NMR      3.636364      8.773679      1.683938      0.000000      8.543165      6.923838      9.172483
    4-P09619-NMR      2.260398     15.926893      1.808318      8.543165      0.000000     18.412162     10.162254
    """

    # save the similarity matrix
    writer = pd.ExcelWriter(ident_matrix_xlsx)
    df_ident.to_excel(writer, sheet_name="sim_matrix")

    df_pairs_above_cutoff = pd.DataFrame()
    df_pairs_above_cutoff["pairs"] = acc_pairs_above_cutoff
    df_pairs_above_cutoff.to_excel(writer, sheet_name="pairs_above_cutoff")

    """ df_pairs_above_cutoff looks like this:
    0  (1-P21709-NMR, 4-P09619-NMR)
    1  (1-P21709-NMR, 5-P22607-NMR)
    2  (4-P09619-NMR, 5-P22607-NMR)
    """

    # create a dictionary with each TMD_name as index, and the cluster of TMDs that share parwise similarity
    TMD_names = df_ident.index.to_list()
    cluster_TMD_dict = {}
    for TMD_name in TMD_names:
        cluster = [TMD_name]
        for pair_tuple in df_pairs_above_cutoff["pairs"]:
            if isinstance(pair_tuple, str):
                pair_tuple = make_tuple(pair_tuple)
            if TMD_name in pair_tuple:
                for paired_TMD_name in pair_tuple:
                    if paired_TMD_name != TMD_name:
                        cluster.append(paired_TMD_name)
        cluster_TMD_dict[TMD_name] = cluster

    """ cluster_TMD_dict looks like this
    {'0-P20963-NMR': ['0-P20963-NMR'], 
    '1-P21709-NMR': ['1-P21709-NMR', '4-P09619-NMR', '5-P22607-NMR']...
    """

    df_clusters = pd.DataFrame()
    for TMD_name, cluster in cluster_TMD_dict.items():
        df_clusters.at[TMD_name, "cluster"] = str(tuple(cluster))
    all_clusters = df_clusters.cluster.unique()
    all_clusters = [make_tuple(t) for t in all_clusters]
    flattened_names_in_all_clusters = [item for sublist in all_clusters for item in sublist]
    # confirm that all TMDs are represented in the list of clusters
    assert set(TMD_names) == set(flattened_names_in_all_clusters)

    df_reduced_clusters = pd.DataFrame()
    reduced_clusters = reduce_clusters_based_on_common_elements(all_clusters, df_ident.index.to_list())
    df_reduced_clusters["reduced_clusters"] = reduced_clusters
    df_reduced_clusters["reduced_cluster_length"] = [len(x) for x in reduced_clusters]
    df_reduced_clusters.to_excel(writer, sheet_name="reduced_clusters")

    """df_reduced_clusters looks like this:
                                 reduced_clusters  reduced_cluster_length
    0                              {0-P20963-NMR}                       1
    1  {1-P21709-NMR, 5-P22607-NMR, 4-P09619-NMR}                       3
    2                              {2-O43914-NMR}                       1
    3                              {3-P05067-NMR}                       1
    4  {1-P21709-NMR, 5-P22607-NMR, 4-P09619-NMR}                       3
    """

    logging.info(f"reduced_cluster_length: {df_reduced_clusters['reduced_cluster_length'].to_list()}")

    # save the settings used
    df_settings = pd.DataFrame()
    df_settings.at["input_fasta", "value"] = protein_set_full_seq_fasta
    df_settings.at["alignment_file", "value"] = output_align
    df_settings.at["sim_matrix_csv", "value"] = ident_matrix_xlsx
    df_settings.at["matrix", "value"] = "matlist.blosum62"
    df_settings.at["gap_open", "value"] = gap_open
    df_settings.at["gap_extend", "value"] = gap_extend
    df_settings.at["acc_pairs_analysed", "value"] = acc_pairs_analysed
    df_settings.to_excel(writer, sheet_name="settings")

    writer.save()
    writer.close()

    logging.info('~~~~~~~~~~~~                 finished create_identity_matrix_using_pairwise_alignments              ~~~~~~~~~~~~')


def reduce_clusters_based_on_common_elements(all_clusters: List[List], TMD_name_list: list):
    n_iterations = 0
    n_clusters_is_constant = False
    output_list: List[Set] = []
    while not n_clusters_is_constant:
        output_list = conduct_one_round_of_cluster_reduction(all_clusters)
        if len(output_list) == len(all_clusters):
            n_clusters_is_constant = True
        else:
            all_clusters = output_list
            n_iterations += 1

    flattened_names_in_all_clusters = [item for sublist in all_clusters for item in sublist]
    # confirm that all TMDs are represented in the list of clusters
    assert set(TMD_name_list) == set(flattened_names_in_all_clusters)
    assert pd.Series(flattened_names_in_all_clusters).value_counts().iloc[0] == 1

    return output_list

def conduct_one_round_of_cluster_reduction(input_clusters):
    output_list: List[Set] = []
    cluster_indices_already_examined: List[int] = []
    for n, outer_cluster in enumerate(input_clusters):
        joined = False
        if n in cluster_indices_already_examined:
            continue
        outer_cluster_set = set(outer_cluster)
        cluster_indices_already_examined.append(n)
        for m, inner_cluster in enumerate(input_clusters):
            if m in cluster_indices_already_examined:
                continue
            inner_cluster_set = set(inner_cluster)
            if outer_cluster_set == inner_cluster_set:
                cluster_indices_already_examined.append(m)
                continue
            intersection = outer_cluster_set.intersection(inner_cluster_set)
            if len(intersection) >= 1:
                new_combined_set_from_2_orig_sets = outer_cluster_set.union(inner_cluster_set)
                output_list.append(new_combined_set_from_2_orig_sets)
                cluster_indices_already_examined.append(m)
                joined = True
                break
        # if none of the elements in outer_set are found in the other sets, then simply append the original outer_set to output file
        if joined == False:
            output_list.append(outer_cluster_set)

    return output_list