BIN_DIR=$(dirname "$0")
CSV_DATA_DIR="data/csv"
#LIMIT=$1

# For development
export IMGT="3420"
export RELEASES="3420"  # this value should be either 3360 or 3370 
export ALIGN=True
export KIR=False

RELEASES=$(echo "${RELEASES}" | sed s'/"//g')
NUM_ALLELES=$(cat data/hla.$RELEASES.dat | grep -c "ID ")
echo "RELEASES: ""$RELEASES"
echo "NUM_ALLELES: $NUM_ALLELES"

# Check RELEASES 
echo "Check releases..."
if [ "$RELEASES" ]; then
	echo "Number of IMGT releases being loaded = ${RELEASES}"
else
	echo "Exiting. RELEASES env variable is not set."
	exit 1
fi

# Loading KIR
echo "Check KIR..."
KIRFLAG=""
if [ "$KIR" == "True" ]; then
	echo "Loading KIR = " "${KIR}"
	KIRFLAG="-k"
fi

echo "Check ALIGN..."
echo $ALIGN
ALIGNFLAG=""
if [ "$ALIGN" == "True" ]; then
	echo "Loading ALIGNMENTS = ${ALIGN}"
	ALIGNFLAG="-a"
	sh "${BIN_DIR}"/get_alignments.sh
fi

echo "Creating new data directory in root..."
mkdir -p "${CSV_DATA_DIR}"

# Profile load script
echo "Building GFE data..."

# Run load script
echo "" > summary_agg.txt
echo "" > summary_diff.txt

# Remove previously created csv files
rm $CSV_DATA_DIR/*.csv

# Build csv files
python3 "${BIN_DIR}"/build_gfedb.py \
	-o "${CSV_DATA_DIR}" \
	-r "${RELEASES}" \
	${KIRFLAG} \
	${ALIGNFLAG} \
	-v \
	-c "${NUM_ALLELES}" \
	-l $1

# # coverage run "${BIN_DIR}"/build_gfedb_optimized.py \
# python3 -m filprofiler run "${BIN_DIR}"/build_gfedb_optimized.py \
# 	-o "${CSV_DATA_DIR}" \
# 	-r "${RELEASES}" \
# 	${KIRFLAG} \
# 	${ALIGNFLAG} \
# 	-v \
# 	-l $1 # > memory_profile_logging.log