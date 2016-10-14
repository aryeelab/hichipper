#!/bin/bash

# Parse parameters
WK_DIR=$1
OUT_NAME=$2
SAMPLE=$3
BAM_ONE=$4
BAM_TWO=$5
PEAK_PAD=$6
MERGE_GAP=$7
MIN_QUAL=30
READ_LEN=75

# Make useful shortcuts
LOG_FILE="${OUT_NAME}/hichipper.log"

echo "`date`: Processing ${SAMPLE}" | tee -a $LOG_FILE

# Process .bam into bedpe
samtools view -F2304 "${WK_DIR}/${BAM_ONE}" | awk -v READ_LEN="$READ_LEN" -v MIN_QUAL="$MIN_QUAL" '{if ($5>=MIN_QUAL) print $3,$4,$4+READ_LEN,$1; else print "*","*","*",$1}' OFS='\t' > "${OUT_NAME}/${SAMPLE}_pos_r1.bed.tmp"
samtools view -F2304 "${WK_DIR}/${BAM_TWO}" | awk -v READ_LEN="$READ_LEN" -v MIN_QUAL="$MIN_QUAL" '{if ($5>=MIN_QUAL) print $3,$4,$4+READ_LEN,$1; else print "*","*","*",$1}' OFS='\t' > "${OUT_NAME}/${SAMPLE}_pos_r2.bed.tmp"
Total_PETs1=`samtools flagstat "${WK_DIR}/${BAM_ONE}" | head -1 | awk '{print $1}'`
Total_PETs2=`samtools flagstat "${WK_DIR}/${BAM_TWO}" | head -1 | awk '{print $1}'`

# The example files that I generated aren't actually paired, so I'm putting this condition in there
if [ "$Total_PETs1" -lt "$Total_PETs2" ]; then
	echo "WARNING-- .bam 1 has more reads than .bam 2. Check to see if they actually match!" | tee -a $LOG_FILE
	head "-${Total_PETs1}" "${OUT_NAME}/${SAMPLE}_pos_r2.bed.tmp" > temppets
	mv temppets "${OUT_NAME}/${SAMPLE}_pos_r2.bed.tmp"
elif [ "$Total_PETs2" -lt "$Total_PETs1" ]; then
	echo "WARNING-- .bam 2 has more reads than .bam 1. Check to see if they actually match!" | tee -a $LOG_FILE
	head "-${Total_PETs2}" "${OUT_NAME}/${SAMPLE}_pos_r1.bed.tmp" > temppets
	mv temppets "${OUT_NAME}/${SAMPLE}_pos_r1.bed.tmp"
fi

paste "${OUT_NAME}/${SAMPLE}_pos_r1.bed.tmp" "${OUT_NAME}/${SAMPLE}_pos_r2.bed.tmp" | cut -f 1-3,5-8 > "${OUT_NAME}/${SAMPLE}_rawinteractions.tmp" 

# Keep only interactions where both reads are mapped
awk '{if ($1 != "*" && $4 != "*") print}' "${OUT_NAME}/${SAMPLE}_rawinteractions.tmp" > "${OUT_NAME}/${SAMPLE}_interactions.unsorted.bedpe.tmp"

# Swap anchors if necessary
awk '{if ($1<$4 || ($1==$4 && $2<=$5)) print $0; else print $4,$5,$6,$1,$2,$3,$7}' 'OFS=\t' "${OUT_NAME}/${SAMPLE}_interactions.unsorted.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}_interactions.unsorted.orderedanchors.bedpe.tmp"
sort -k1,1n -k2,2n -k4,4n -k5,5n "${OUT_NAME}/${SAMPLE}_interactions.unsorted.orderedanchors.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}_interactions.bedpe.tmp"
Mapped_PETs_q30=`wc -l "${OUT_NAME}/${SAMPLE}_interactions.bedpe.tmp" | awk '{print $1}'`

# Remove duplicates
echo "`date`: Writing unique (i.e. deduplicated) interactions (based on chr1, pos1, chr2, pos2)" | tee -a $LOG_FILE
sort -k1,1n -k2,2n -k4,4n -k5,5n --unique "${OUT_NAME}/${SAMPLE}_interactions.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp"
Mapped_unique_PETs_q30=`wc -l "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" | awk '{print $1}'`
Mapped_unique_intrachromosomal_q30=`awk '$1 == $4 {print $0}' "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" | wc -l | awk '{print $1}'`
Mapped_unique_intrachromosomal_q30_5kb=`awk '$1 == $4 && (($5+$6) - (($2+$3)/2)> 5000) {print $0}' "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" | wc -l | awk '{print $1}'`

echo "Mapped_unique_PETs_q30=${Mapped_unique_PETs_q30}" 
echo "Mapped_unique_intrachromosomal_q30=${Mapped_unique_intrachromosomal_q30}"
echo "Mapped_unique_intrachromosomal_q30_5kb=${Mapped_unique_intrachromosomal_q30_5kb}" 


# Call peaks
cut -f 1-3,7 "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" >  "${OUT_NAME}/${SAMPLE}_left.dedup.bed.tmp"
cut -f 4-6,7 "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" >  "${OUT_NAME}/${SAMPLE}_right.dedup.bed.tmp"
cat "${OUT_NAME}/${SAMPLE}_left.dedup.bed.tmp" "${OUT_NAME}/${SAMPLE}_right.dedup.bed.tmp" >> "${OUT_NAME}/${SAMPLE}_reads.bed.tmp"	
macs2 callpeak -t "${OUT_NAME}/${SAMPLE}_reads.bed.tmp" -g hs -f BED -n "${OUT_NAME}/${SAMPLE}_temporary" -p 0.01 --nomodel

# Pad peaks
awk -v PEAK_PAD="$PEAK_PAD" '{$2-=PEAK_PAD; $3+=PEAK_PAD}1' OFS="\t" "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak" > "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded"
echo "`date`: Padded peaks by ${PEAK_PAD}bp" | tee -a $LOG_FILE
awk '$2 > 0 && $5 > 0 {print $0}' "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded" > "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded.clean"

# Merge gaps
bedtools merge -d $MERGE_GAP -i "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded.clean" > "${OUT_NAME}/${SAMPLE}_temporary_peaks.merged.bed.tmp"
bedtools intersect -loj -a "${OUT_NAME}/${SAMPLE}_left.dedup.bed.tmp" -b "${OUT_NAME}/${SAMPLE}_temporary_peaks.merged.bed.tmp" | awk '{print $5,$6,$7,$4}' OFS='\t' > "${OUT_NAME}/${SAMPLE}_anchor1.bed.tmp"
bedtools intersect -loj -a "${OUT_NAME}/${SAMPLE}_right.dedup.bed.tmp" -b "${OUT_NAME}/${SAMPLE}_temporary_peaks.merged.bed.tmp" | awk '{print $5,$6,$7,$4}' OFS='\t' > "${OUT_NAME}/${SAMPLE}_anchor2.bed.tmp"

# Make final output
paste "${OUT_NAME}/${SAMPLE}_anchor1.bed.tmp" "${OUT_NAME}/${SAMPLE}_anchor2.bed.tmp" | cut -f 1-3,5-8 > "${OUT_NAME}/${SAMPLE}_anchor.interactions.tmp"
awk '{if ($1 != "." && $4 != ".") print}' "${OUT_NAME}/${SAMPLE}_anchor.interactions.tmp" > "${OUT_NAME}/${SAMPLE}_anchor.interactions.bedpe.tmp"
cut -f1-6  "${OUT_NAME}/${SAMPLE}_anchor.interactions.bedpe.tmp" | sort | uniq -c | awk '{print $2,$3,$4,$5,$6,$7,".",$1}' >  "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp"
awk '$1 != $4 {print $0}' "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}.inter.loop_counts.bedpe"
awk '$1 == $4 && $2 != $5 {print $0}' "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}.intra.loop_counts.bedpe"

# Write out summary statistics
echo "Total_PETs=${Total_PETs1}" > "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_PETs_q30=${Mapped_PETs_q30}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_PETs_q30=${Mapped_unique_PETs_q30}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intrachromosomal_q30=${Mapped_unique_intrachromosomal_q30}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intrachromosomal_q30_5kb=${Mapped_unique_intrachromosomal_q30_5kb}" >> "${OUT_NAME}/${SAMPLE}.stat"
