#!/usr/bin/env bash

set -x

BIN=$(dirname "$0")

DATA_DIR=${BIN}/../data
mkdir -p "${DATA_DIR}"

#base_url="https://media.githubusercontent.com/media/ANHIG/IMGTHLA"
base_url="https://raw.githubusercontent.com/ANHIG/IMGTHLA"

RELEASES=`echo ${RELEASES} | sed s'/"//'g | sed s'/,/ /g'`
echo "ALIGN RELEASES = ${RELEASES}"

loci="A B C DRB1 DQB1 DPB1 DPA1 DQA1"

for dbversion in ${RELEASES};do
  dbversion_trimmed=$(echo ${dbversion} | sed 's/\.//g')
  mkdir -p "${DATA_DIR}/${dbversion_trimmed}"
	for loc in ${loci}
	do
		msf_url="${base_url}/${dbversion_trimmed}/msf/${loc}_gen.msf"
		curl -L "${msf_url}" -o "${DATA_DIR}/${dbversion_trimmed}/${loc}_gen.msf"

		msf_url="${base_url}/${dbversion_trimmed}/msf/${loc}_nuc.msf"
		curl -L "${msf_url}" -o "${DATA_DIR}/${dbversion_trimmed}/${loc}_nuc.msf"

		msf_url="${base_url}/${dbversion_trimmed}/msf/${loc}_prot.msf"
		curl -L "${msf_url}" -o "${DATA_DIR}/${dbversion_trimmed}/${loc}_prot.msf"
	done
done

if [ "$KIR" == "True" ]; then
	kir_loci="'KIR3DS1 KIR3DP1 KIR3DL3 KIR3DL2 KIR3DL1 KIR2DS5 KIR2DS4 KIR2DS3 KIR2DS2 KIR2DS1 KIR2DP1 KIR2DL4"
	kir_base="ftp://ftp.ebi.ac.uk/pub/databases/ipd/kir/msf"

	mkdir -p "${DATA_DIR}/kir/${dbversion_trimmed}"

	for kir_locus in ${kir_loci}
	do
		kir_url="${kir_base}/${kir_locus}_gen.msf"
		curl -L "${kir_url}" -o "${DATA_DIR}/kir/${kir_locus}_gen.msf"
	done
fi

