#!/bin/bash

START_EXECUTION=$SECONDS

ROOT=$(dirname $(dirname "$0"))
BIN_DIR=$ROOT/bin
SRC_DIR=$ROOT/src
DATA_DIR=$ROOT/data

if [ ! -d "$DATA_DIR" ]; then
	echo "Creating new data directory in root..."
	mkdir -p $DATA_DIR/csv/
else
	# Remove previously created csv files
	rm -r $DATA_DIR/csv/*.csv
fi

# # For development
# export RELEASES="3420, 3430"
# export ALIGN=True
# export KIR=False
# export MEM_PROFILE=False

# Check if RELEASES is set
if [ -z ${RELEASES+x} ]; then 
    echo "RELEASES is not set. Please specify the release versions to load."; 
else 
    echo "Loading IMGT/HLA release versions: $RELEASES";
fi

# #RELEASES=$(echo "$RELEASES" | sed s'/"//g')
# echo "IMGT versions: $RELEASES"

# Load KIR data
echo "Check KIR..."
KIRFLAG=""
if [ "$KIR" == "True" ]; then
	echo "Loading KIR = $KIR"
	KIRFLAG="-k"
fi

# Load alignments data
echo "Check ALIGN..."
echo $ALIGN
ALIGNFLAG=""
if [ "$ALIGN" == "True" ]; then
	echo "Loading ALIGNMENTS = $ALIGN"
	ALIGNFLAG="-a"
	sh $BIN_DIR/get_alignments.sh
fi

# Memory profiling
MEM_PROFILE_FLAG=""
if [ "$MEM_PROFILE" == "True" ]; then
	echo "Memory profiling is set to $MEM_PROFILE."
	MEM_PROFILE_FLAG="-p"
	echo "" > summary_agg.txt
	echo "" > summary_diff.txt
fi

# # TO DO: Handle downloading dat files outside the python script
# for release in $RELEASES; do
# 	if [ ! -f "$DATA_DIR/hla.$release.dat" ]; then
# 		echo "Downloading DAT file for IPD-IMGT/HLA version $release..."
# 		curl -o hla.$release.dat -k https://media.githubusercontent.com/media/ANHIG/IMGTHLA/3420/hla.dat
# 	else
# 		echo "DAT file for IPD-IMGT/HLA version $release is already present..."
# 	fi
# done

# Build csv files
# $RELEASES=${echo "$RELEASES" | sed s'/,//g'}
for release in $RELEASES; do

	# # TO DO: handle downloading DAT files outside python script
	# NUM_ALLELES=$(cat $DATA_DIR/hla.$release.dat | grep -c "ID ")
	# echo "ALLELES: $NUM_ALLELES"
	
	release=$(echo "$release" | sed s'/,//g')

	echo -e "\n"
	python3 "$SRC_DIR"/gfedb.py \
		-o "$DATA_DIR/csv" \
		-r "$release" \
		$KIRFLAG \
		$ALIGNFLAG \
		$MEM_PROFILE_FLAG \
		-v \
		-l $1
		# -c "$NUM_ALLELES" \
	echo -e "\n"
done

END_EXECUTION=$(( SECONDS - start_release ))
echo "Finished in $END_EXECUTION seconds"