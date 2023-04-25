
import os

try:
    from numpy.distutils.core import Extension, setup
except ImportError:
    raise ImportError("NumPy is not installed. Please install NumPy manually before.")


# get text of README.md
current_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(current_path, "README.md")) as f:
    readme_text = f.read()

setup(
    name="rmsd-clustering",
    version="0.1.0",
    description="Clustering based on RMSD",
    long_description=readme_text,
    long_description_content_type="text/markdown",
    url="https://github.com/boneta/RMSD-Clustering",
    author="Sergio Boneta",
    author_email="boneta@unizar.es",
    license="GPLv3",
    python_requires='>=3.8',
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics"
        ],
    setup_requires=["numpy"],
    install_requires=[
        "numpy",
        "scipy",
        "ase",
        "rmsd",
        "matplotlib",
        "seaborn"
        ],
    ext_modules=[
        Extension(
            name="rmsd_clustering.rmsd_fort",
            sources=["rmsd_clustering/rmsd_fort.f90"],
            extra_f90_compile_args=["-O3", "-fopenmp"],
            extra_link_args=["-lgomp", "-lblas", "-llapack"]
            )
        ],
    packages=["rmsd_clustering"],
    include_package_data=False,
    entry_points={
        "console_scripts": ["rmsd_clustering=rmsd_clustering.__main__:main"]
        }
)
