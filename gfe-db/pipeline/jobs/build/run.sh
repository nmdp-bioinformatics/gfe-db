#!/bin/bash

START_EXECUTION=$SECONDS

export ROOT="$(dirname "$(dirname "$0")")"
export BIN_DIR=$ROOT/scripts
export SRC_DIR=$ROOT/src
export DATA_DIR=$ROOT/../data
export LOGS_DIR=$ROOT/logs

# Check for environment variables
if [[ -z "${GFE_BUCKET}" ]]; then
	echo "GFE_BUCKET not set. Please specify an S3 bucket."
	exit 1
fi

if [[ -z "${RELEASES}" ]]; then
	echo "RELEASES not set. Please specify the release versions to load."
	exit 1
fi

if [[ -z "${ALIGN}" ]]; then
	echo "ALIGN not set"
	ALIGN=False
fi

if [[ -z "${KIR}" ]]; then
	echo "KIR not set"
	KIR=False
fi

if [[ -z "${MEM_PROFILE}" ]]; then
	echo "MEM_PROFILE not set"
	MEM_PROFILE=False
fi

echo "Found environment variables:"
echo -e "GFE_BUCKET: $GFE_BUCKET\nRELEASES: $RELEASES\nALIGN: $ALIGN\nKIR: $KIR\nMEM_PROFILE: $MEM_PROFILE\nLIMIT: $LIMIT"

# Check limit
if [[ -z "${LIMIT}" ]]; then
	echo "No limit set, building GFEs for all alleles"
else
	echo "Build is limited to $LIMIT alleles"
fi

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
	echo "Creating new directory in root: $DATA_DIR"
	mkdir -p "$DATA_DIR"
else
	# TODO: get full path
	echo "Data directory: $DATA_DIR"
fi

# Check if logs directory exists
if [ ! -d "$LOGS_DIR" ]; then
	# TODO: get full path
	echo "Creating logs directory: $LOGS_DIR"
	mkdir -p "$LOGS_DIR"
else
	# TODO: get full path
	echo "Logs directory: $LOGS_DIR"
fi

# Memory profiling
if [ "$MEM_PROFILE" == "True" ]; then
	echo "Memory profiling is set to $MEM_PROFILE."
	MEM_PROFILE_FLAG="-p"
	touch "$LOGS_DIR/mem_profile_agg.txt"
	touch "$LOGS_DIR/mem_profile_diff.txt"
else
	MEM_PROFILE_FLAG=""
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
	sh "$BIN_DIR/get_alignments.sh"
else
	ALIGNFLAG=""
fi

# Build csv files
RELEASES=$(echo "${RELEASES}" | sed s'/"//'g | sed s'/,/ /g')
# exit 1 # TODO test state machine error handling
for release in ${RELEASES}; do

	release=$(echo "$release" | sed s'/,//g')
	echo "Processing release: $release"

	# Check if data directory exists
	if [ ! -d "$DATA_DIR/$release/csv" ]; then
		# TODO: get full path
		echo "Creating new directory in root: $DATA_DIR/$release/csv..."
		mkdir -p "$DATA_DIR/$release/csv"
	else
		# TODO: get full path
		echo "CSV directory: $DATA_DIR/$release/csv"
	fi

	# Check if DAT file exists
	if [ -f "$DATA_DIR/$release/hla.$release.dat" ]; then
		echo "DAT file for release $release already exists"
	else
		echo "Downloading DAT file for release $release..."
		if [ "$(echo "$release" | bc -l)" -lt 3350  ]; then

			# Should be environment variable
			imgt_hla_raw_url='https://raw.githubusercontent.com/ANHIG/IMGTHLA'
			echo "Downloading $imgt_hla_raw_url/$release/hla.dat to $DATA_DIR/$release/hla.$release.dat"
			curl -SL "$imgt_hla_raw_url/$release/hla.dat" > "$DATA_DIR/$release/hla.$release.dat"
		else
			imgt_hla_media_url='https://media.githubusercontent.com/media/ANHIG/IMGTHLA'
			echo "Downloading $imgt_hla_media_url/$release/hla.dat to $DATA_DIR/$release/hla.$release.dat"
			curl -SL "$imgt_hla_media_url/$release/hla.dat" > "$DATA_DIR/$release/hla.$release.dat"
		fi
	fi
	
	# Builds CSV files
	python3 "$SRC_DIR"/app.py \
		-o "$DATA_DIR/$release/csv" \
		-r "$release" \
		$KIRFLAG \
		$ALIGNFLAG \
		$MEM_PROFILE_FLAG \
		-v \
		-l $LIMIT
    build_exit_status=$?
    echo "Build exit status (1:CRITICAL, 2:WARNING): $build_exit_status"
    
    # Notify missing alleles
    if [ $build_exit_status -eq 2 ]; then
    echo "WARNING: Some alleles failed to build, please see logs for error messages"
    fi

    # fail for any exit code other than 0 or 2. 2 is a warning for missing data but not fatal.
    if [ $build_exit_status -ne 0 ] && [ $build_exit_status -ne 2 ]; then
    echo "CRITICAL: Build failed, please see logs for error messages"
    exit 1
    fi

	# TODO: Use this S3 hierarchy: root/release/csv | logs
	echo -e "Uploading CSVs to s3://$GFE_BUCKET/data/$release/csv/:\n$(ls $DATA_DIR/$release/csv/)"
	aws s3 --recursive cp "$DATA_DIR/$release/csv/" s3://$GFE_BUCKET/data/$release/csv/ > "$LOGS_DIR/s3Copy$$LOG_FILE"
	mv "$LOGS_DIR/gfeBuildLogs.txt" "$LOGS_DIR/gfeBuildLogs.$release.txt"
	mv "$LOGS_DIR/s3Copy$$LOG_FILE" "$LOGS_DIR/s3CopyLog.$release.txt"

	if [ "$MEM_PROFILE" == "True" ]; then
		mv "$LOGS_DIR/mem_profile_agg.txt" "$LOGS_DIR/mem_profile_agg.$release.txt"
		mv "$LOGS_DIR/mem_profile_diff.txt" "$LOGS_DIR/mem_profile_diff.$release.txt"
	fi

	echo -e "Uploading logs to s3://$GFE_BUCKET/logs/$release/:\n$(ls $LOGS_DIR/)"
	aws s3 --recursive cp "$LOGS_DIR/" s3://$GFE_BUCKET/logs/pipeline/build/$release/ > "$LOGS_DIR/s3CopyLog.Local.txt"

done

END_EXECUTION=$(( SECONDS - $START_EXECUTION ))
echo "Finished in $END_EXECUTION seconds"
exit 0

# For debugging to keep the build server running
# sleep 1h
