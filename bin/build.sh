#!/usr/bin/env bash

BIN=$(dirname "$0")

sh ${BIN}/get_alignments.sh

python3 ${BIN}/build_gfedb.py -o $1 -k -v

sh ${BIN}/load_graph.sh $1

