BIN_DIR=$(dirname "$0")
SRC_DIR=$(dirname $(dirname "$0"))/src
DATA_DIR=$(dirname $(dirname "$0"))/data

# For development
export RELEASES="3420 3430"  # this value should be either 3360 or 3370 
export ALIGN=True
export KIR=False

RELEASES=$(echo "$RELEASES" | sed s'/"//g')
echo "IMGT versions: $RELEASES"

# Loading KIR
echo "Check KIR..."
KIRFLAG=""
if [ "$KIR" == "True" ]; then
	echo "Loading KIR = $KIR"
	KIRFLAG="-k"
fi

echo "Check ALIGN..."
echo $ALIGN
ALIGNFLAG=""
if [ "$ALIGN" == "True" ]; then
	echo "Loading ALIGNMENTS = $ALIGN"
	ALIGNFLAG="-a"
	sh $BIN_DIR/get_alignments.sh
fi

if [ ! -d "$DIRECTORY" ]; then
	echo "Creating new data directory in root..."
	mkdir -p $DATA_DIR
else
	# Remove previously created csv files
	rm $DATA_DIR/*.csv
fi

# Run load script
echo "" > summary_agg.txt
echo "" > summary_diff.txt

# Build csv files
for release in $RELEASES; do

	echo "Building GFE data for version $release..."
	NUM_ALLELES=$(cat $DATA_DIR/hla.$release.dat | grep -c "ID ")
	echo "Total alleles: $NUM_ALLELES"

	python3 "$SRC_DIR"/gfedb.py \
		-o "$CSV_DATA_DIR" \
		-r "$release" \
		$KIRFLAG \
		$ALIGNFLAG \
		-v \
		-c "$NUM_ALLELES" \
		-l $1
done
