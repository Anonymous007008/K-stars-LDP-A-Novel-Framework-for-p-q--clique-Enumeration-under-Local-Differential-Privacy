#!/bin/bash -x

if [ $# -ne 1 ]; then
    echo "USAGE: run_EvalLossLocal.sh [Dataset]"
    exit 1
fi

./k-stars_LDP ../data/${1}/edges.csv 1000 1-0.001 2 6-150 3-1 2
./k-stars_LDP ../data/${1}/edges.csv 2000 1-0.001 2 6-150 3-1 2
./k-stars_LDP ../data/${1}/edges.csv 3000 1-0.001 2 6-150 3-1 2
./k-stars_LDP ../data/${1}/edges.csv 4000 1-0.001 2 6-150 3-1 2
./k-stars_LDP ../data/${1}/edges.csv 5000 1-0.001 2 6-150 3-1 2
./k-stars_LDP ../data/${1}/edges.csv 6000 1-0.001 2 6-150 3-1 2
./k-stars_LDP ../data/${1}/edges.csv 7000 1-0.001 2 6-150 3-1 2
./k-stars_LDP ../data/${1}/edges.csv 8000 1-0.001 2 6-150 3-1 2
./k-stars_LDP ../data/${1}/edges.csv 9000 1-0.001 2 6-150 3-1 2
./k-stars_LDP ../data/${1}/edges.csv 10000 1-0.001 2 6-150 3-1 2

