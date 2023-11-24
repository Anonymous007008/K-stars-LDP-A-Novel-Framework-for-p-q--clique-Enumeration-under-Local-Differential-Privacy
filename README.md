# K-stars-LDP-A-Novel-Framework-for-p-q--clique-Enumeration-under-Local-Differential-Privacy
## Overview
This is the official code of k-stars LDP. Datasets in the experiment should be put in the directory: ./data/<dataset>/ . The source code of k-stars LDP is put in the directory: ./cpp/

## Requirements
The gcc version we use in the experiment is 11.4.0.

## Environments
The experiments are conducted on the machine with Intel(R) Xeon(R) Gold 6230R CPU @ 2.10GHz, and NVIDIA GeForce RTX 3090 with 24GB memory and CUDA 11.8. The operating system is Ubuntu 18.04.6 with 216GB memory.

## Datasets
The datasets we use in the experiments are Gplus, IMDB, GitHub and Facebook. Please put then into the directory: ./data/<dataset>/.

## Run the k-stars LDP algorithms
Please run the following commands (<dataset> is Gplus, IMDB, GitHub or Facebook):
```
$ cd cpp/
$ make
$ chmod +x run_k-stars_LDP.sh
$ ./run_k-stars_LDP.sh [Dataset (Gplus/IMDB/GitHub/Facebook)]
$ cd ../
```
