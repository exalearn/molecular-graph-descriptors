# Molecular Graph Descriptors

This codebase accompanies the paper A Look Inside the Black Box: Using Graph-Theoretical Descriptors to Interpret a Continuous-Filter Convolutional Neural Network (CF-CNN) trained on the Global and Local Minimum Energy Structures of Neutral Water Clusters.

Molecular graphs and structures of over 5 million water clusters can be found [here](https://sites.uw.edu/wdbase/).

## Contents

### graphdescriptors

Code for generating graphs from molecular structures and computing various descriptors.

### schnetpack

Code amended from [SchNetPack](https://github.com/atomistic-machine-learning/schnetpack) to train and test the CF-CNN.

## Requirements

* python 3
* pytorch (>= 0.4.1)
* h5py
* ASE
* networkx
* pandas
* numpy

## References

* Jenna A. Bilbrey, Joseph P. Heindel, Malachi Schram, Pradipta Bandyopadhyay, Soritis S. Xantheas, and Sutanay Choudhury. A Look Inside the Black Box: Using Graph-Theoretical Descriptors to Interpret a Continuous-Filter Convolutional Neural Network (CF-CNN) trained on the Global and Local Minimum Energy Structures of Neutral Water Clusters. *Journal of Chemical Physics*, **2020**.
