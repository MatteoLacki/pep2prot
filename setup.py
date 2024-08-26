# This Python file uses the following encoding: utf-8
from setuptools import find_packages, setup

setup(
    name="pep2prot",
    packages=find_packages(),
    version="0.1.0",
    description="Turn a peptide report into a protein report.",
    long_description="Turn a peptide report into a protein report.",
    author="Mateusz Krzysztof Łącki",
    author_email="matteo.lacki@gmail.com",
    url="https://github.com/MatteoLacki/pep2prot.git",
    keywords=["mass spec", "peptide protein inference" "proteomics"],
    classifiers=[
        "Development Status :: 1 - Planning",
        "License :: OSI Approved :: BSD License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    install_requires=[
        "numpy",
        "networkx",  # TODO: try to get rid off
        "pandas",
        "furious_fastas",  # TODO: only need a fasta sequence reader really.
        "numba",
        "numba-progress",
        "tqdm",
        "duckdb",
    ],
    # include_package_data=True,
    # package_data={
    #     'data':
    #          ['data/contaminants_uniprot_format.fasta']
    # },
    # scripts = [
    #     "bin/pep2prot",
    #     "bin/_pep2prot.bat"
    # ]
)
