BIN_DIR=$(dirname "$0")
CSV_DATA_DIR="../data/csv"

RELEASES=$(echo "${RELEASES}" | sed s'/"//g')

# Check RELEASES 
if [ "$RELEASES" ]; then
	echo "Number of IMGT releases being loaded = " "${RELEASES}"
else
	echo "Exiting. RELEASES env variable is not set."
	exit 1
fi

# Loading KIR
KIRFLAG=""
if [ "$KIR" == "True" ]; then
	echo "Loading KIR = " "${KIR}"
	KIRFLAG="-k"
fi

ALIGNFLAG=""
if [ "$ALIGN" == "True" ]; then
	echo "Loading ALIGNMENTS = " "${ALIGN}"
	ALIGNFLAG="-a"
	sh "${BIN_DIR}"/get_alignments.sh
fi

mkdir -p "${CSV_DATA_DIR}"
python3 "${BIN_DIR}"/build_gfedb.py -o "${CSV_DATA_DIR}" -r "${RELEASES}" ${KIRFLAG} ${ALIGNFLAG} -v
