
import itertools
import multiprocessing as mp
from copy import deepcopy

import ase
import ase.io
import numpy as np
import rmsd

from . import rmsd_fort


class MolecGroup:
    """
    Class to hold a group of molecules

    Attributes
    ----------
    Atoms : list
        List of ASE Atoms objects
    """

    reorder_methods = {
        'hungarian': rmsd.reorder_hungarian,
        'intertia-hungarian': rmsd.reorder_inertia_hungarian,
        'brute': rmsd.reorder_brute,
        'distance': rmsd.reorder_distance}

    def __init__(self):
        self.Atoms: list[ase.Atoms] = []

    def read(self, filename, format=None) -> None:
        """
        Read a file and add it to the molecule group

        Parameters
        ----------
        filename : str
            File to read
        format : str, optional
            File format to use, otherwise inferred from file extension
        """
        n_atoms = self.n_atoms
        for f in ase.io.read(filename, index=":", format=format):
            if f.get_global_number_of_atoms() != n_atoms != 0:
                raise ValueError("Number of atoms in new file does not match previous files")
            self.Atoms.append(f)

    @property
    def n_molecs(self) -> int:
        """Number of molecules in the molecule group"""
        return len(self.Atoms)

    @property
    def n_atoms(self) -> int:
        """Number of atoms in every molecule of the molecule group"""
        return self.Atoms[0].get_global_number_of_atoms() if self.Atoms else 0

    @property
    def coords(self) -> np.ndarray:
        """Coordinates of all atoms in the molecule group"""
        return np.array([m.get_positions() for m in self.Atoms], dtype=np.float64)

    @property
    def zatoms(self) -> np.ndarray:
        """Atomic numbers of all atoms in the molecule group"""
        return np.array([m.get_atomic_numbers() for m in self.Atoms], dtype=np.int32)

    def copy_ndx(self, ndx: list) -> 'MolecGroup':
        """
        Partial copy of the self object based on a list of molecules

        Parameters
        ----------
        ndx : list[int]
            List of indexes to copy

        Returns
        -------
        MolecGroup
            New molecule group with the selected molecules
        """
        new_cluster = MolecGroup()
        new_cluster.Atoms = deepcopy([self.Atoms[i] for i in ndx])
        return new_cluster

    def set_centroid(self) -> None:
        """Move all the molecules in the molecule group to their centroid"""
        for m in self.Atoms:
            m.positions -= m.get_center_of_mass()

    @property
    def consistent_ordered(self) -> bool:
        """Check if the atoms in each molecule of the molecule group are equally ordered"""
        return np.all(np.all(self.zatoms, axis=1))

    def sort_atoms(self) -> None:
        """Sort the atoms in each molecule of the molecule group by atomic number"""
        for m in self.Atoms:
            ndx_sort = np.argsort(m.get_atomic_numbers())
            m.set_atomic_numbers(m.get_atomic_numbers()[ndx_sort])
            m.set_positions(m.get_positions()[ndx_sort])

    def rmsd_ndx(self, ndx1:int, ndx2:int, align:bool=True, reorder_method=None) -> float:
        """
        Calculate the RMSD between two molecules in the molecule group

        Parameters
        ----------
        ndx1 : int
            Index of the first molecule
        ndx2 : int
            Index of the second molecule
        align : bool, optional
            Align the molecules before calculating the RMSD (def: True)
        reorder_method : {None, 'hungarian', 'intertia-hungarian', 'brute', 'distance'}, optional
            Method to reorder the atoms (def: None)

        Returns
        -------
        float
            RMSD between the two molecules
        """
        zatoms = self.zatoms
        coords = self.coords
        molec1 = (zatoms[ndx1], coords[ndx1])
        molec2 = (zatoms[ndx2], coords[ndx2])
        if reorder_method is not None:
            molec2 = self.reorder(molec1, molec2, reorder_method)
        rmsd_func = rmsd_fort.rmsd_kabsch if align else rmsd_fort.rmsd
        return rmsd_func(molec1[1], molec2[1])

    def rmsd_matrix(self, align:bool=True, reorder_method=None) -> np.ndarray:
        """
        Calculate the RMSD matrix between all molecules in the molecule group

        Parameters
        ----------
        align : bool, optional
            Align the molecules before calculating the RMSD (def: True)
        reorder_method : {None, 'hungarian', 'intertia-hungarian', 'brute', 'distance'}, optional
            Method to reorder the atoms (def: None)

        Returns
        -------
        np.ndarray
            RMSD matrix
        """
        if reorder_method is not None:
            pairs = list(itertools.combinations(range(self.n_molecs), 2))
            with mp.Pool(processes=mp.cpu_count()) as pool:
                rmsds = pool.starmap(self.rmsd_ndx, [(i, j, align, reorder_method) for i, j in pairs])
            matrix = np.zeros((self.n_molecs, self.n_molecs))
            for pair, rmsd in zip(pairs, rmsds):
                matrix[pair] = rmsd
            matrix += matrix.T
        else:
            matrix = rmsd_fort.rmsd_matrix(self.coords, align)
        return matrix

    @staticmethod
    def reorder(molec1:tuple, molec2:tuple, reorder_method:str) -> tuple:
        """
        Reorder the atoms in molecule 2 to match molecule 1

        Parameters
        ----------
        molec1 : tuple[ndarray, ndarray]
            Reference molecule, tuple of (zatoms, coordinates)
        molec2 : tuple[ndarray, ndarray]
            Molecule to reorder, tuple of (zatoms, coordinates)
        reorder_method : {'hungarian', 'intertia-hungarian', 'brute', 'distance'}, optional
            Method to reorder the atoms

        Returns
        -------
        tuple[ndarray, ndarray]
            The reordered molecule, tuple of (atoms, coordinates)
        """
        #TODO: Faster implementation in Fortran
        if reorder_method not in MolecGroup.reorder_methods:
            raise ValueError('Unknown reorder method')
        molec2_reorder = MolecGroup.reorder_methods[reorder_method](molec1[0], molec2[0], molec1[1], molec2[1])
        return molec2[0][molec2_reorder], molec2[1][molec2_reorder]
