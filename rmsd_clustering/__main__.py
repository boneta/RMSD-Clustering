#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np

from . import __version__
from . import clustering, matrix_tools
from .molec_group import MolecGroup


def main():

    parser = argparse.ArgumentParser(
        prog="rmsd_clustering",
        description="Clustering based on RMSD",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"RMSD Clustering - v{__version__}")
    parser.add_argument(
        "input", metavar="INPUT", type=str, nargs="+",
        help="Input file(s) to be processed")
    parser.add_argument(
        "-i", "--in-type", metavar="<>", type=str,
        help="Input file format to apply to all inputs, otherwise inferred\nAvailable formats can be check with 'ase info --formats'")
    parser.add_argument(
        "-r", "--reorder", metavar="<>", type=str, default="none",
        choices=['none', 'hungarian', 'intertia-hungarian', 'brute', 'distance'],
        help="method for reordering the atoms (none/hungarian/inertia-hungarian/brute/distance) (def: %(default)s)")
    parser.add_argument(
        "-c", "--cutoff", metavar="#", type=float, default=1.0,
        help="RMSD cutoff for clustering (default: %(default)s)")
    parser.add_argument(
        "-n", "--n-clusters", metavar="#", type=int, default=None,
        help="Number of clusters to generate, overrides cutoff")
    parser.add_argument(
        "--heatmap", metavar="<>", type=str, default=None,
        help="save the RMSD matrix as a heatmap image to a file")
    parser.add_argument(
        "--rmsd-matrix", metavar="<>", type=str, default=None,
        help="read/write the RMSD matrix from/to a .csv file")
    args = parser.parse_args()
    input_files = args.input
    input_type = args.in_type
    reorder_method = args.reorder if args.reorder != 'none' else None
    cutoff = args.cutoff
    n_clusters = args.n_clusters
    heatmap = args.heatmap
    rmsd_matrix_file = args.rmsd_matrix

    ## read input files
    m = MolecGroup()
    print(f'\n ## Files:\n', flush=True)
    for input in input_files:
        n_molecs_prev = m.n_molecs
        m.read(input, format=input_type)
        n_molecs_last = m.n_molecs - n_molecs_prev
        ending_molec = '' if n_molecs_last == 1 else f' - {m.n_molecs}'
        print(f'    {input} ({n_molecs_last}):  {n_molecs_prev + 1}{ending_molec}', flush=True)
    m.set_centroid()

    ## RMSD matrix
    print(f'\n ## RMSD matrix', end="", flush=True)
    if rmsd_matrix_file and Path(rmsd_matrix_file).exists():
        print(f' (loaded from \'{rmsd_matrix_file}\')', end="", flush=True)
        rmsd_matrix = np.loadtxt(rmsd_matrix_file, delimiter=',')
        if not m.n_molecs == rmsd_matrix.shape[0] == rmsd_matrix.shape[1]:
            raise ValueError(f'Number of molecules in the input files ({m.n_molecs}) does not match the RMSD matrix ({rmsd_matrix.shape[0]}x{rmsd_matrix.shape[1]})')
    else:
        print(f' (reorder method = {reorder_method})', end="", flush=True)
        rmsd_matrix = m.rmsd_matrix(align=True, reorder_method=reorder_method)
        if rmsd_matrix_file:
            print(f' - saved to \'{rmsd_matrix_file}\'', end="", flush=True)
            np.savetxt(rmsd_matrix_file, rmsd_matrix, delimiter=',')
    rmsd_matrix_stat = matrix_tools.stat(rmsd_matrix)
    print(f':  {rmsd_matrix_stat[0]:.4f} +/- {rmsd_matrix_stat[1]:.4f}     (mean +/- std  intra-RMSD)')
    if heatmap:
        print(f'\n ## Heatmap image saved to \'{heatmap}\'\n', flush=True)
        matrix_tools.heatmap(rmsd_matrix, filename=heatmap)

    ## clustering
    print(f'\n ## Clusters', end="", flush=True)
    if n_clusters:
        cutoff, clusters = clustering.cluster_n(rmsd_matrix, n_clusters, cutoff)
    else:
        clusters = clustering.cluster_cutoff(rmsd_matrix, cutoff)
    print(f' (cutoff = {cutoff:.4f}):\n', flush=True)
    for i, cluster_ndx in enumerate(clusters):
        cluster_rmsd_matrix = matrix_tools.submatrix(rmsd_matrix, cluster_ndx)
        cluster_stat = matrix_tools.stat(cluster_rmsd_matrix)
        cluster_ndx = ', '.join([str(ndx+1) for ndx in cluster_ndx])
        print(f'{i+1:6d}:  {cluster_ndx:50s}     {cluster_stat[0]:.4f} +/- {cluster_stat[1]:.4f}')
    print()


if __name__ == "__main__":
    main()
