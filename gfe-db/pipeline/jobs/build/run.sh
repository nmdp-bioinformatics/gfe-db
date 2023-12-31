#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

START_EXECUTION=$SECONDS

# export ROOT="$(dirname "$(dirname "$0")")"
export ROOT="$(dirname "$0")"
export BIN_DIR=$ROOT/scripts
export SRC_DIR=$ROOT/src
export DATA_DIR=$ROOT/data
export LOGS_DIR=$ROOT/logs

get_download_url() {
    owner="$1"
    repo="$2"
    asset_path="$3"
    commit_sha="$4"

    base_url="https://api.github.com"
    endpoint="/repos/${owner}/${repo}/contents/${asset_path}"
    url="${base_url}${endpoint}"

    # # Authorization header
    # auth_header="Authorization: token ${GITHUB_PERSONAL_ACCESS_TOKEN}"

    # Content-Type header
    content_type_header="Content-Type: application/json"

    # Accept header
    accept_header="Accept: application/vnd.github.v3+json"

    # X-GitHub-Api-Version header
    x_github_api_version_header="X-GitHub-Api-Version: 2022-11-28"

    # GET request with headers and ref parameter
    # response=$(curl -s -H "${auth_header}" -H "${content_type_header}" -H "${accept_header}" -H "${x_github_api_version_header}" "${url}?ref=${commit_sha}")
    response=$(curl -s -H "${content_type_header}" -H "${accept_header}" -H "${x_github_api_version_header}" "${url}?ref=${commit_sha}")

    # Print the response
    echo "${response}" | jq -r '.download_url'

}

# Takes the download url and downloads the asset
get_asset() {
    download_url="$1"
    asset_path="$2"

    # Download the asset
    response=$(curl -sSL -o "${asset_path}" -w "%{http_code}" "${download_url}")

    if [ "$response" -ne 200 ]; then
        echo "Failed to download asset. HTTP response code: $response"
        exit 1
    else
        echo "Successfully downloaded asset. HTTP response code: $response"
    fi
}

if [[ -z "${EVENT}" ]]; then
	echo "No event found. Exiting..."
    exit 1
else
	echo "Found event"
    echo "$EVENT"
fi

# parse event
version=$(echo "$EVENT" | jq -r '.state.execution.version')
commit_sha=$(echo "$EVENT" | jq -r '.state.commit.sha')
repository_owner=$(echo "$EVENT" | jq -r '.state.repository.owner')
repository_name=$(echo "$EVENT" | jq -r '.state.repository.name')
align=$(echo "$EVENT" | jq -r '.state.execution.input_parameters.align')
kir=$(echo "$EVENT" | jq -r '.state.execution.input_parameters.kir')
mem_profile=$(echo "$EVENT" | jq -r '.state.execution.input_parameters.mem_profile')
limit=$(echo "$EVENT" | jq -r '.state.execution.input_parameters.limit')
s3_path=$(echo "$EVENT" | jq -r '.input.s3_path')

# Refactor the above variable validations into a for loop
for var in version commit_sha align kir mem_profile limit repository_owner repository_name s3_path; do
    if [[ -z "${!var}" ]] || [[ "${!var}" == "null" ]]; then
        echo "\`$var\` not set. Please specify a value."
        exit 1
    fi
    echo "$var: ${!var}"
done

if [[ "${limit}" == "-1" ]] || [[ -z "${limit}" ]]; then
    echo "No limit set, building GFEs for all alleles"
elif [[ "${limit}" =~ ^[0-9]+$ ]] && [[ "${limit}" -gt 0 ]]; then
    echo "Build is limited to $limit alleles"
else
    echo "Invalid limit specified. Please specify either a positive integer or -1 for no limit."
    exit 1
fi

echo "Found environment variables"

# Check if data directory exists
# TODO: get full path for each
if [ ! -d "$DATA_DIR" ]; then
	echo "Creating new directory in root: $DATA_DIR"
	mkdir -p "$DATA_DIR"
else
	echo "Data directory: $DATA_DIR"
fi

# Check if logs directory exists
if [ ! -d "$LOGS_DIR" ]; then
	echo "Creating logs directory: $LOGS_DIR"
	mkdir -p "$LOGS_DIR"
else
	echo "Logs directory: $LOGS_DIR"
fi

# TODO test memory profiling for build job
# Memory profiling
if [ "$mem_profile" == "true" ]; then
	echo "Memory profiling is set to $mem_profile."
	MEM_PROFILE_FLAG="-p"
	touch "$LOGS_DIR/mem_profile_agg.txt"
	touch "$LOGS_DIR/mem_profile_diff.txt"
else
	MEM_PROFILE_FLAG=""
fi

# Load kir data
if [ "$kir" == "true" ]; then
	echo "Loading kir = $kir"
	KIRFLAG="-k"
else
	KIRFLAG=""
fi

# Load alignments data
if [ "$align" == "true" ]; then
	echo "Loading alignments..."
	ALIGNFLAG="-a"
	sh "$BIN_DIR/get_alignments.sh"
else
	ALIGNFLAG=""
fi

# Check for FEATURE_SERVICE_URL
if [[ -z "${FEATURE_SERVICE_URL}" ]]; then
	echo "No FEATURE_SERVICE_URL set, building GFEs with default feature service."
else
	echo "Using Feature Service: ${FEATURE_SERVICE_URL}"
fi

# Build csv files
RELEASES=$(echo "${RELEASES}" | sed s'/"//'g | sed s'/,/ /g')
# exit 1 # TODO test state machine error handling

echo "Processing release version: $version"

# Check if data directory exists
# TODO: get full path for each
if [ ! -d "$DATA_DIR/$version/csv" ]; then
    echo "Creating new directory in root: $DATA_DIR/$version/csv"
    mkdir -p "$DATA_DIR/$version/csv"
else
    echo "CSV directory: $DATA_DIR/$version/csv"
fi

# Check if DAT file exists
if [ -f "$DATA_DIR/$version/hla.$version.dat" ]; then
    echo "DAT file for release $version already exists"
else

    # download_url works for all releases including 3440 and earlier
    echo "Fetching DAT file for release $version..."
    download_url="$(get_download_url "$repository_owner" "$repository_name" "hla.dat" "$commit_sha")"
    get_asset "$download_url" "$DATA_DIR/$version/hla.$version.dat"
fi

# Builds CSV files
# TODO booleans for kir, align, mem_profile are lower case, limit is now -1 instead of none
# TODO implement s3_path
python3 "$SRC_DIR"/app.py \
	-o "$DATA_DIR/$version/csv" \
	-r "$version" \
	$KIRFLAG \
	$ALIGNFLAG \
	$MEM_PROFILE_FLAG \
	-v \
	-l $limit \
    -u $FEATURE_SERVICE_URL
    
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

# TODO: Use this S3 hierarchy: root/release/data/csv | logs
echo -e "Uploading data to s3://$GFE_BUCKET/data/$version"
res=$(aws s3 cp --recursive "$DATA_DIR/$version/" "s3://$GFE_BUCKET/data/$version/")

echo $res

if [ "$mem_profile" == "true" ]; then
	mv "$LOGS_DIR/mem_profile_agg.txt" "$LOGS_DIR/mem_profile_agg.$version.txt"
	mv "$LOGS_DIR/mem_profile_diff.txt" "$LOGS_DIR/mem_profile_diff.$version.txt"
fi

echo -e "Uploading logs to s3://$GFE_BUCKET/logs/$version"
aws s3 --recursive cp "$LOGS_DIR/" s3://$GFE_BUCKET/logs/pipeline/build/$version/logs/


END_EXECUTION=$(( SECONDS - $START_EXECUTION ))
echo "Finished in $END_EXECUTION seconds"

# For debugging to keep the build server running
if [ "$DEBUG" == "true" ]; then
    echo "DEBUG mode is set to $DEBUG. Sleeping..."
    while true; do sleep 1000; done
fi

exit 0