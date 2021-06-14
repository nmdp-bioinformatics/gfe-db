#!/bin/bash

START_EXECUTION=$SECONDS

ROOT=$(dirname $(dirname "$0"))
BIN_DIR=$ROOT/scripts
SRC_DIR=$ROOT/src
DATA_DIR=$ROOT/data
LOGS_DIR=$ROOT/logs
# GFE_BUCKET=gfe-db-4498

# aws stepfunctions get-activity-task ...

# For development
export RELEASES="3410"
export ALIGN=True
export KIR=False
export MEM_PROFILE=True

# #RELEASES=$(echo "$RELEASES" | sed s'/"//g')
# echo "IMGT versions: $RELEASES"

# Check if RELEASES is set
if [ -z ${RELEASES+x} ]; then 
    echo "RELEASES is not set. Please specify the release versions to load."; 
else 
    echo "Loading IMGT/HLA release versions: $RELEASES";
fi

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
	echo "Creating new directory in root: $DATA_DIR"
	mkdir -p $DATA_DIR
else
	echo "Data directory: $DATA_DIR"
fi

# Check if data directory exists
if [ ! -d "$LOGS_DIR" ]; then
	echo "Creating new directory in root: $LOGS_DIR"
	mkdir -p $LOGS_DIR
	touch $LOGS_DIR/logs.txt
else
	rm -f $LOGS_DIR/logs.txt
fi

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

# Build csv files
# $RELEASES=${echo "$RELEASES" | sed s'/,//g'}

# rm -rf $DATA_DIR/csv/*.csv

for release in $RELEASES; do

	release=$(echo "$release" | sed s'/,//g')

	# Check if data directory exists
	if [ ! -d "$DATA_DIR/$release/csv" ]; then
		echo "Creating new directory in root: $DATA_DIR/$release/csv..."
		mkdir -p $DATA_DIR/$release/csv
	else
		echo "CSV directory: $DATA_DIR/$release/csv"
	fi

	echo
	echo "Building release: $release..."

	# NUM_ALLELES=$(cat $DATA_DIR/hla.$release.dat | grep -c "ID ")
	# echo "ALLELES: $NUM_ALLELES"

	if [ -f $DATA_DIR/$release/hla.$release.dat ]; then
		echo "DAT file for release $release already exists"
	else
		echo "Downloading DAT file for release $release..."
		if [ "$(echo "$release" | bc -l)" -le 3350  ]; then
			imgt_hla_raw_url='https://raw.githubusercontent.com/ANHIG/IMGTHLA'
			echo "Downloaded: $imgt_hla_raw_url/$release/hla.dat to ..."
			curl -L $imgt_hla_raw_url/$release/hla.dat > $DATA_DIR/$release/hla.$release.dat
		else
			imgt_hla_media_url='https://media.githubusercontent.com/media/ANHIG/IMGTHLA'
			echo "Downloaded: $imgt_hla_media_url/$release/hla.dat to ..."
			curl -L $imgt_hla_media_url/$release/hla.dat > $DATA_DIR/$release/hla.$release.dat
		fi
	fi
	
	# echo -e "\n"
	python3 "$SRC_DIR"/gfedb.py \
		-o "$DATA_DIR/$release/csv" \
		-r "$release" \
		$KIRFLAG \
		$ALIGNFLAG \
		$MEM_PROFILE_FLAG \
		-v \
		-l $1
		# -c "$NUM_ALLELES" \
	# echo -e "\n"

	# # Copy CSVs to S3
	echo "Copying CSVs to S3..."
	aws s3 --recursive cp $DATA_DIR/$release/csv/ s3://$GFE_BUCKET/data/$release/csv/
done

END_EXECUTION=$(( SECONDS - $START_EXECUTION ))
echo "Finished in $END_EXECUTION seconds"

# Publish to SNS
