# RMSD Clustering

*Clustering based on RMSD*

## Usage
```
  rmsd_clustering [-h] [options] INPUT [INPUT ...]
```

## Installation
```
  pip install git+https://github.com/boneta/RMSD-Clustering.git
```
### Requirements
- Fortran compiler (tested with *gfortran*)
- LAPACK-devel libraries
- *Python 3.8+* and the following packages:
  - *numpy*
  - *scipy*
  - *ase*
  - *rmsd*
  - *matplotlib*
  - *seaborn*

## Detailed usage

### Input
Input files can be of any kind supported by [ASE](https://wiki.fysik.dtu.dk/ase/) (e.g. *xyz*, *pdb*, *vasp*). The input can be either a single file or several. File formats can be mixed and they will be automatically inferred individually.  
A format can be forced by using the *--in-type* option, in which case all input files must be of the same format.

### Reordering
Atom order is assumed to be the same for all structures. If this is not the case, the *--reorder* option can be used to reorder the atoms using different methods: 'hungarian', 'inertia-hungarian', 'brute' or 'distance'.

### Clustering
The clustering is performed based on the calculation of an RMSD matrix between all pairs of structures. Then they are grouped if the RMSD is below a given threshold set with the *--cutoff* option. If *--n-clusters* is used instead, the structures will be grouped into *n* clusters setting the cutoff accordingly.

### Examples
```
  rmsd_clustering molecs.xyz
  rmsd_clustering POSCAR_01 POSCAR_02 POSCAR_03
  rmsd_clustering --in-type vasp --reorder brute CONTCAR_*
  rmsd_clustering --n-clusters 20 molecs.xyz
```
