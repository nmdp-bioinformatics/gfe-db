#!/usr/bin/env bash

BIN=$(dirname "$0")

# Check RELEASES 
if [ "$RELEASES" ]; then
	echo "Number of IMGT releases being loaded = " ${RELEASES}
else
	echo "Exiting. RELEASES env variable is not set."
	exit
fi

# Loading KIR
KIRFLAG=""
if [ "$KIR" == "True" ]; then
	echo "Loading KIR = " ${KIR}
	KIRFLAG="-k"
fi

sh ${BIN}/get_alignments.sh

python3 ${BIN}/build_gfedb.py -o $1 -r ${RELEASES} ${KIRFLAG} -v

