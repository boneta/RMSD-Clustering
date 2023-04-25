
F2PY = python3 -m numpy.f2py --f90flags="-O3 -fopenmp" -lgomp -lblas -llapack 

all:
	$(MAKE) clean
	$(F2PY) -c rmsd_clustering/rmsd_fort.f90 -m rmsd_clustering.rmsd_fort

clean:
	rm -f *.so