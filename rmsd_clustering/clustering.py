
import numpy as np


def cluster_cutoff(rmsd_matrix:np.ndarray, cutoff:float) -> list:
    """
    Perform a simple clustering over a rmsd matrix using a cutoff

    Parameters
    ----------
    rmsd_matrix : ndarray, optional
        RMSD matrix to use instead of calculating (def: None)
    cutoff : float
        RMSD cutoff for clustering

    Returns
    -------
    list[list[int]]
        List of cluster's indexes
    """
    clusters = []
    for i in range(rmsd_matrix.shape[0]):
        for cluster in clusters:
            if all(rmsd_matrix[i,j] < cutoff for j in cluster):
                cluster.append(i)
                break
        else:
            clusters.append([i])
    return clusters

def cluster_n(rmsd_matrix:np.ndarray, n:int, cutoff:float=0.1) -> tuple:
    """
    Perform a simple clustering over a rmsd matrix to obtain n clusters

    Parameters
    ----------
    rmsd_matrix : ndarray, optional
        RMSD matrix to use instead of calculating (def: None)
    n : int
        Number of clusters to obtain
    cutoff : float, optional
        Initial RMSD cutoff for clustering (def: 0.1)


    Returns
    -------
    (float, list[list[int]])
        Final RMSD cutoff used for clustering
        List of cluster's indexes
    """
    if not 1 <= n <= rmsd_matrix.shape[0]:
        raise ValueError("'n' must be between 1 and the number of molecules")
    clusters = cluster_cutoff(rmsd_matrix, cutoff)
    while len(clusters) != n:
        if len(clusters) > n:
            cutoff += cutoff * 0.5
        else:
            cutoff -= cutoff * 0.5
        clusters = cluster_cutoff(cutoff=cutoff, rmsd_matrix=rmsd_matrix)
    return cutoff, clusters
