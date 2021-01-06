# TrophicNetEncoder

These R scripts reproduce the simulated data and analysis of Botella et al. - 2021 - An appraisal of graph embeddings for analysing architectural variation across trophic networks.

## Required packages

igraph
network
intergraph
ggplot2
sna
GGally
reticulate
umap (R installation only)
moments
randomForest
energy



pcName =as.character(Sys.info()["nodename"]) 
if(pcName=="DESKTOP-RUARS8N"){
  use_python("C:/Users/user/miniconda3/python.exe",required = T)
  source_python('C:/Users/user/pCloud local/boulot/Github/EcoGraph Encoder/embed_reduce_analyse_wrappers/py_get_ebd.py')
}else if(pcName=="PORTHOS"){
  use_python("/home/christophe/miniconda3/bin/python",required = T)
  source_python('/home/christophe/pCloud local/boulot/Github/EcoGraph Encoder/embed_reduce_analyse_wrappers/py_get_ebd.py')
}

## Graph2Vec installation

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

## Procedure

All files of this repository should be contained in a same directory. First, run simulate_nets.R to generate the dataset of networks, then compute_embeddings.R to compute the embedding of the networks with each method, then run_analysis_pipeline.R to reproduce the tables and Figures of the manuscript.
