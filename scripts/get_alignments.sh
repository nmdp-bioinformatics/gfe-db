#!/bin/bash

BIN_DIR=$(dirname "$0")

DATA_DIR=${BIN_DIR}/../data
mkdir -p "${DATA_DIR}"

base_url="https://raw.githubusercontent.com/ANHIG/IMGTHLA"

RELEASES=`echo ${RELEASES} | sed s'/"//'g | sed s'/,/ /g'`
echo "ALIGN RELEASES = ${RELEASES}"

loci="A B C DRB1 DQB1 DPB1 DPA1 DQA1"

for dbversion in ${RELEASES}; do
  dbversion_trimmed=$(echo ${dbversion} | sed 's/\.//g')
  mkdir -p "${DATA_DIR}/${dbversion_trimmed}/alignments"
	for loc in ${loci}
	do

		gen_msf="${DATA_DIR}/${dbversion_trimmed}/alignments/${loc}_gen.msf"
		nuc_msf="${DATA_DIR}/${dbversion_trimmed}/alignments/${loc}_nuc.msf"
		prot_msf="${DATA_DIR}/${dbversion_trimmed}/alignments/${loc}_prot.msf"
	
		if [ ! -f $gen_msf ]; then
			msf_url="${base_url}/${dbversion_trimmed}/msf/${loc}_gen.msf"
			echo "Downloading ${loc}_gen.msf..."
			curl -L "${msf_url}" -o "${gen_msf}"
		else
   			echo "${loc}_gen.msf already exists, skipping download..."
		fi

		if [ ! -f $nuc_msf ]; then
			msf_url="${base_url}/${dbversion_trimmed}/msf/${loc}_nuc.msf"
			echo "Downloading ${loc}_nuc.msf..."
			curl -L "${msf_url}" -o "${nuc_msf}"
		else
   			echo "${loc}_nuc.msf already exists, skipping download..."
		fi

		if [ ! -f $prot_msf ]; then
			msf_url="${base_url}/${dbversion_trimmed}/msf/${loc}_prot.msf"
			echo "Downloading ${loc}_prot.msf..."
			curl -L "${msf_url}" -o "${prot_msf}"
		else
   			echo "${loc}_prot.msf already exists, skipping download..."
		fi
	done
done

if [ "$KIR" = "True" ]; then
	kir_loci="'KIR3DS1 KIR3DP1 KIR3DL3 KIR3DL2 KIR3DL1 KIR2DS5 KIR2DS4 KIR2DS3 KIR2DS2 KIR2DS1 KIR2DP1 KIR2DL4"
	kir_base="ftp://ftp.ebi.ac.uk/pub/databases/ipd/kir/msf"

	mkdir -p "${DATA_DIR}/${dbversion_trimmed}/kir"

	for kir_locus in ${kir_loci}
	do
		kir_url="${kir_base}/${kir_locus}_gen.msf"
		curl -L "${kir_url}" -o "${DATA_DIR}/${dbversion_trimmed}/kir/${kir_locus}_gen.msf"
	done
fi

