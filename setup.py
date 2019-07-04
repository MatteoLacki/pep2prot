# This Python file uses the following encoding: utf-8
from setuptools import setup, find_packages

setup(
    name='pep2prot',
    packages=find_packages(),
    version='0.0.1',
    description='Turn a peptide report into a protein report.',
    long_description='Turn a peptide report into a protein report.',
    author='Mateusz Krzysztof Łącki',
    author_email='matteo.lacki@gmail.com',
    url='https://github.com/MatteoLacki/pep2prot.git',
    keywords=[
        'Mass Spectrometry',
        'peptide protein inference'
        'proteomics'],
    classifiers=[
        'Development Status :: 1 - Planning',
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'],
    install_requires=[
        'networkx',
    ],
    # include_package_data=True,
    # package_data={
    #     'data':
    #          ['data/contaminants_uniprot_format.fasta']
    # },
    scripts = [
        "bin/pep2prot"
    ]
)
