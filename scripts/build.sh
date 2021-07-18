#!/bin/bash

START_EXECUTION=$SECONDS

ROOT=$(dirname $(dirname "$0"))
BIN_DIR=$ROOT/scripts
SRC_DIR=$ROOT/src
DATA_DIR=$ROOT/data
LOGS_DIR=$ROOT/logs

export CYPHER_PATH=neo4j/cypher
# SCRIPT=load.cyp
export SCRIPT=load.cyp

# For development
export GFE_BUCKET=gfe-db-4498
export RELEASES="3420"
export ALIGN=True
export KIR=False
export MEM_PROFILE=True

# Check for environment variables
if [[ -z "${GFE_BUCKET}" ]]; then
	echo "GFE_BUCKET not set"
	exit 1
elif [[ -z "${RELEASES}" ]]; then
	echo "RELEASES not set. Please specify the release versions to load."
	exit 1
elif [[ -z "${ALIGN}" ]]; then
	echo "ALIGN not set"
	exit 1
elif [[ -z "${KIR}" ]]; then
	echo "KIR not set"
	exit 1
elif [[ -z "${MEM_PROFILE}" ]]; then
	echo "MEM_PROFILE not set"
	exit 1
else
	echo "Found environment variables:"
	echo -e "GFE_BUCKET: $GFE_BUCKET\nRELEASES: $RELEASES\nALIGN: $ALIGN\nKIR: $KIR\nMEM_PROFILE: $MEM_PROFILE"
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
if [ "$KIR" == "True" ]; then
	echo "Loading KIR = $KIR"
	KIRFLAG="-k"
else
	KIRFLAG=""
fi

# Load alignments data
if [ "$ALIGN" == "True" ]; then
	echo "Loading alignments..."
	ALIGNFLAG="-a"
	sh $BIN_DIR/get_alignments.sh
else
	ALIGNFLAG=""
fi

# Memory profiling
if [ "$MEM_PROFILE" == "True" ]; then
	echo "Memory profiling is set to $MEM_PROFILE."
	MEM_PROFILE_FLAG="-p"
	echo "" > summary_agg.txt
	echo "" > summary_diff.txt
else
	MEM_PROFILE_FLAG=""
fi

# Build csv files
RELEASES=`echo "${RELEASES}" | sed s'/"//'g | sed s'/,/ /g'`

for release in ${RELEASES}; do

	release=$(echo "$release" | sed s'/,//g')
	echo "Processing release: $release"

	# Check if data directory exists
	if [ ! -d "$DATA_DIR/$release/csv" ]; then
		echo "Creating new directory in root: $DATA_DIR/$release/csv..."
		mkdir -p $DATA_DIR/$release/csv
	else
		echo "CSV directory: $DATA_DIR/$release/csv"
	fi

	# Check if DAT file exists
	if [ -f $DATA_DIR/$release/hla.$release.dat ]; then
		echo "DAT file for release $release already exists"
	else
		echo "Downloading DAT file for release $release..."
		if [ "$(echo "$release" | bc -l)" -le 3350  ]; then
			imgt_hla_raw_url='https://raw.githubusercontent.com/ANHIG/IMGTHLA'
			echo "Downloading $imgt_hla_raw_url/$release/hla.dat to $DATA_DIR/$release/hla.$release.dat"
			curl -SL $imgt_hla_raw_url/$release/hla.dat > $DATA_DIR/$release/hla.$release.dat
		else
			imgt_hla_media_url='https://media.githubusercontent.com/media/ANHIG/IMGTHLA'
			echo "Downloading $imgt_hla_media_url/$release/hla.dat to $DATA_DIR/$release/hla.$release.dat"
			curl -SL $imgt_hla_media_url/$release/hla.dat > $DATA_DIR/$release/hla.$release.dat
		fi
	fi
	
	# python3 "$SRC_DIR"/gfedb.py \
	# 	-o "$DATA_DIR/$release/csv" \
	# 	-r "$release" \
	# 	$KIRFLAG \
	# 	$ALIGNFLAG \
	# 	$MEM_PROFILE_FLAG \
	# 	-v \
	# 	-l $1

	# Copy CSVs to S3
	echo -e "Uploading CSVs to s3://$GFE_BUCKET/data/$release/csv/:\n$(ls $DATA_DIR/$release/csv/)"
	aws s3 --recursive --quiet cp $DATA_DIR/$release/csv/ s3://$GFE_BUCKET/data/$release/csv/

	# # ls data/$release/csv
	# echo "Creating pre-signed URLs..."
	# export urls=()
	# for filename in ./data/$release/csv/*.csv; do
	# 	filename=`basename $filename`
	# 	# echo $filename
	# 	url=$(aws s3 presign --expires-in 3600 s3://$GFE_BUCKET/data/$release/csv/$filename)
	# 	urls+=($(echo $url | sed 's/\&/\\&/'))
	# done

	# echo "URLS:"
	# printf '%s\n' "${urls[@]}"

	# # Update URLs in Cypher script
	# sed -i.bak "s+file:///all_alignments.RELEASE.csv+"${urls[0]}"+g" $CYPHER_PATH/$SCRIPT
	# sed -i.bak "s+file:///all_cds.RELEASE.csv+"${urls[1]}"+g" $CYPHER_PATH/$SCRIPT
	# sed -i.bak "s+file:///all_features.RELEASE.csv+"${urls[2]}"+g" $CYPHER_PATH/$SCRIPT
	# sed -i.bak "s+file:///all_groups.RELEASE.csv+"${urls[3]}"+g" $CYPHER_PATH/$SCRIPT
	# sed -i.bak "s+file:///gfe_sequences.RELEASE.csv+"${urls[4]}"+g" $CYPHER_PATH/$SCRIPT

	sh $BIN_DIR/load_db.sh


done

END_EXECUTION=$(( SECONDS - $START_EXECUTION ))
echo "Finished in $END_EXECUTION seconds"
