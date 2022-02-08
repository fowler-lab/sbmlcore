# sbml-core
Collection of core classes and functions for structure-based machine learning to predict antimicrobial resistance


## Included features

### Secondary structure 

STRIDE and DSSP. Do not anticipate much difference.

### Effect of mutation on protein stability

There are numerous options here:

1. StabilityPredict. Online metapredictor, single amino acid at a time. Josh used in the pncA paper but had to contact them directly to run the entirity of PncA. (paper)[https://www.jbc.org/article/S0021-9258(20)34176-4/fulltext]
2. DeepDDG. More recent neural network that claims to outperform DUET, PopMusic etc. (paper)[https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.8b00697] (server)[http://protein.org.cn/ddg.html]. Can do all possible mutations in one job.
3. DynaMUT. Also claims to outperform DUET etc. (paper)[https://academic.oup.com/nar/article/46/W1/W350/4990022]. Can process a list of specified mutations in one job. (server)[http://biosig.unimelb.edu.au/dynamut/]
