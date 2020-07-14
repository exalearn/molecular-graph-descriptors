# Molecular Graph Descriptors

This codebase accompanies the paper A Look Inside the Black Box: Using Graph-Theoretical Descriptors to Interpret a Continuous-Filter Convolutional Neural Network (CF-CNN) trained on the Global and Local Minimum Energy Structures of Neutral Water Clusters.

Molecular graphs and structures of over 5 million water clusters can be found [here](https://sites.uw.edu/wdbase/).

## Contents

### data

.xyz files of the lowest-energy structure of each cluster size *N*=3-30 and a set of 15 clusters and their predited potential energy.

### graphdescriptors

Code for generating graphs from molecular structures and computing various descriptors from .xyz files output from SchNetPack. The script get_descriptors.py collects a set of descriptors for all molecules in the .xyz file and outputs a csv file.

```
python get_descriptors.py --data_path ../data/test_predictions.xyz --output test_df.csv --min_dir ../data/
```

### schnetpack

Code amended from [SchNetPack](https://github.com/atomistic-machine-learning/schnetpack) to train and test the CF-CNN. See schnetpack/README.md for use.

## Requirements

* python 3
* pytorch (>= 0.4.1)
* h5py
* ASE
* networkx
* pandas
* numpy

## References
Please cite this paper if you find the code and datasets useful.

@article{bilbrey2020descriptors,
author = {Bilbrey,Jenna A.  and Heindel,Joseph P.  and Schram,Malachi  and Bandyopadhyay,Pradipta  and Xantheas,Sotiris S.  and Choudhury,Sutanay },
title = {A look inside the black box: Using graph-theoretical descriptors to interpret a Continuous-Filter Convolutional Neural Network (CF-CNN) trained on the global and local minimum energy structures of neutral water clusters},
journal = {The Journal of Chemical Physics},
volume = {153},
number = {2},
year = {2020}
}
