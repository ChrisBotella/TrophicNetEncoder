# TrophicNetEncoder

These R scripts reproduce the simulated data and analysis of Botella et al. - 2021 - An appraisal of graph embeddings for analysing architectural variation across trophic networks.

## Required R packages

igraph, network, intergraph, ggplot2, sna, GGally, reticulate, umap (R installation only), moments, randomForest, energy

## Graph2Vec python installation

Graph2Vec is implemented in the python package karateclub (https://github.com/benedekrozemberczki/karateclub) and is called from R through the reticulate R package. We used karateclub version 1.0.14 with python 3.7.6, tested in Windows 10 and Ubuntu 18.04.

To install Graph2Vec, one must:
- Install python and packages numpy, networkx, karateclub (easy using pip, see https://github.com/benedekrozemberczki/karateclub).
- A small fix must be applied to the Graph2Vec implementation. In Lib/site-packages/karateclub/, open **estimator.py** and add "#" in the beginning of lines 74 and 75, as follow:
```
#self._check_connectivity(graph)
#self._check_directedness(graph)
```
- In the R script **functions_to_source.R**, specify the location of the python installation for R-reticulate, example:
```
use_python("C:/Users/JeanMichMich/miniconda3/python.exe",required = T)
```
and also specify the location of the python script **py_get_ebd.py** (contains the functions to load for R-reticulate), example:
```
source_python('C:/Users/JeanMichMich/py_get_ebd.py')
```

## How to reproduce results

All files must be put in one directory which must be specified in the beginning of each R script. Then simply run **simulate_nets.R** to generate the dataset of networks, then run **compute_embeddings.R** to compute the embedding of the networks with each method, then run **run_analysis_pipeline.R** to reproduce the tables and Figures of the manuscript and run **wood2015_Embedding.R** to reproduce analysis on the empirical marine foodwebs.
