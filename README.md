[![Tests](https://github.com/fowler-lab/sbmlcore/actions/workflows/tests.yaml/badge.svg)](https://github.com/fowler-lab/sbmlcore/actions/workflows/tests.yaml)
[![codecov](https://codecov.io/gh/fowler-lab/sbmlcore/branch/main/graph/badge.svg?token=P44BPYQBFS)](https://codecov.io/gh/fowler-lab/sbmlcore)
[![PyPI version](https://badge.fury.io/py/sbmlcore.svg)](https://badge.fury.io/py/sbmlcore)

# sbmlcore
Collection of core classes to help with building structure- and chemistry-based feature datasets to train machine learning models to predict antimicrobial resistance.

This is a pre-release alpha version - it may not be fully functional for your requirements and it is also subject to change with no notice!

We will be making a series of jupyter-notebooks demonstrating how to use the classes available [here](https://github.com/fowler-lab/sbmlcore-tutorials).

## Included features

### Changes in Amino Acid Properties

- Volume
- Hydropathy scales: Kyte-Doolittle ([paper](https://www.sciencedirect.com/science/article/abs/pii/0022283682905150?via%3Dihub)) and WimleyWhite ([paper](https://www.nature.com/articles/nsb1096-842))
- Molecular weight
- Isoelectric point

### Secondary structure

- STRIDE ([website](http://webclu.bio.wzw.tum.de/stride/) and [paper](https://onlinelibrary.wiley.com/doi/10.1002/prot.340230412))

### Solvent accessible surface areas

- FreeSASA ([paper](https://f1000research.com/articles/5-189/v1))
- STRIDE ([website](http://webclu.bio.wzw.tum.de/stride/) and [paper](https://onlinelibrary.wiley.com/doi/10.1002/prot.340230412))

### Likelihood of changes in protein function

- SNAP2 ([server](https://www.rostlab.org/services/snap/) and [paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-16-S8-S1))

### Effect of mutation on protein stability

- DeepDDG: a more recent neural network that claims to outperform DUET, PopMusic etc. ([paper](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00697) and  [server](http://protein.org.cn/ddg.html)). Can do all possible mutations in one job.

### Structural distances

- Distances between mutated residues and any atom/group of atoms of interest. Uses MDAnalysis ([paper1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3144279/) and [paper2](https://conference.scipy.org/proceedings/scipy2016/oliver_beckstein.html)).

## To potentially include at a later stage
 - Secondary structure: DSSP (do not anticipate much difference to STRIDE)
 - Protein stability:
    1. StabilityPredict. Online metapredictor, single amino acid at a time. Josh used in the pncA paper but had to contact them directly to run the entirity of PncA. ([paper](https://www.jbc.org/article/S0021-9258(20)34176-4/fulltext))
    2. DynaMUT. Also claims to outperform DUET etc. ([paper](https://academic.oup.com/nar/article/46/W1/W350/4990022)). Can process a list of specified mutations in one job. ([server](http://biosig.unimelb.edu.au/dynamut/))

PWF, 9 May 2023