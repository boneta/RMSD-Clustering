
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def fmt(matrix:np.ndarray, header_init:int=1) -> str:
    """
    Format a matrix as a string

    Parameters
    ----------
    matrix : ndarray
        Matrix to format
    header_init : int, optional
        Initial value for the header (def: 1)

    Returns
    -------
    str
        Matrix as a formatted string
    """
    s = ' '*7 + ' '.join([f'{i:8d}' for i in range(header_init, header_init + matrix.shape[0])])
    for i in range(matrix.shape[0]):
        s += f'\n{i+header_init:6d} ' + ' '.join([f' {matrix[i,j]:7.4f}' for j in range(matrix.shape[1])])
    return s

def stat(matrix:np.ndarray) -> tuple:
    """
    Calculate statistics of a matrix

    Parameters
    ----------
    matrix : ndarray
        Matrix to calculate the statistics for

    Returns
    -------
    tuple
        Mean of upper triangle (without diagonal), standard deviation of upper triangle (without diagonal)
    """
    if matrix.size == 1:
        mean = 0.0
        std = 0.0
    else:
        mean = np.mean(matrix[np.triu_indices(matrix.shape[0], k=1)])
        std = np.std(matrix[np.triu_indices(matrix.shape[0], k=1)])
    return mean, std

def submatrix(matrix:np.ndarray, ndxs:np.ndarray) -> np.ndarray:
    """
    Extract a submatrix from a matrix

    Parameters
    ----------
    matrix : ndarray
        Matrix to extract the submatrix from
    ndxs : ndarray
        Indices of the submatrix

    Returns
    -------
    ndarray
        Submatrix
    """
    return matrix[np.ix_(ndxs, ndxs)]

def heatmap(matrix:np.ndarray, filename:str, header_init:int=1) -> None:
    """
    Generate a heatmap of a matrix and save it to a file

    Parameters
    ----------
    matrix : ndarray
        The matrix to format
    filename : str
        Filename to save the heatmap image to. File extension is added automatically if not present (def: 'png')
    header_init : int, optional
        The initial value for the header (def: 1)
    """
    filename = filename if Path(filename).suffix[1:] in plt.gcf().canvas.get_supported_filetypes() else filename + '.png'
    #TODO: better handlig of tick labels
    sns.heatmap(
        matrix,
        xticklabels=range(header_init, header_init + matrix.shape[0]),
        yticklabels=range(header_init, header_init + matrix.shape[0]),
        linewidth=0.0,
        cmap='viridis',
        cbar_kws={'label': 'RMSD'},
        square=True
        )
    plt.savefig(filename, dpi=300)
