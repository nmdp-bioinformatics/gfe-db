#!/usr/bin/env bash

BIN=$(dirname "$0")

RELEASES=`echo ${RELEASES} | sed s'/"//g'`

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

ALIGNFLAG=""
if [ "$ALIGN" == "True" ]; then
	echo "Loading ALIGNMENTS = " ${ALIGN}
	ALIGNFLAG="-a"
	sh ${BIN}/get_alignments.sh
fi

#modified the offending lines of C 05:206, 05:208N that don't conform or are missing some information
cp ./mod-imgt/C_gen.sth ./data/3360


python3 ${BIN}/build_gfedb.py -o $1 -r ${RELEASES} ${KIRFLAG} ${ALIGNFLAG} -v
