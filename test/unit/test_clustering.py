import pytest
from typing import List, Set

from thoipapy.clustering.pairwise_aln_similarity_matrix import reduce_clusters_based_on_common_elements


def test_reduce_clusters_based_on_common_elements():
    test_clusters: List[List[int]] = [[1, 2], [2, 3], [4], [5], [2, 6], [3, 7]]
    output_list: List[Set] = reduce_clusters_based_on_common_elements(test_clusters, [1,2,3,4,5,6,7])
    assert output_list == [{1, 2, 3, 6, 7}, {4}, {5}]

    test_clusters: List[List[int]] = [[1, 2, 3]]
    output_list: List[List] = reduce_clusters_based_on_common_elements(test_clusters, [1,2,3])
    assert output_list == [{1, 2, 3}]
