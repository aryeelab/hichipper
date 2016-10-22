#!/bin/bash

# Parse parameters
WK_DIR=$1
OUT_NAME=$2
SAMPLE=$3
BAM_ONE=$4
BAM_TWO=$5
PEAK_PAD=$6
MERGE_GAP=$7
MIN_QUAL=$8
READ_LEN=$9
MIN_DIST=${10}
MAX_DIST=${11}
MACS2_STRING=${12}

# Make useful shortcuts
LOG_FILE="${OUT_NAME}/${OUT_NAME}.hichipper.log"

hichipper --version | tee -a $LOG_FILE
echo "`date`: Processing ${SAMPLE}" | tee -a $LOG_FILE

# Process .bam into bedpe
if [ ! -f "${BAM_ONE}" ] ; then
    echo "File ${BAM_ONE} is not in ${OUT_NAME}, aborting." | tee -a $LOG_FILE
    exit
fi

if [ ! -f "${BAM_TWO}" ] ; then
    echo "File ${BAM_TWO} is not in ${OUT_NAME}, aborting." | tee -a $LOG_FILE
    exit
fi

# Create bed files of the reads, but we remove those with poor quality. 
bedtools bamtobed -i "${BAM_ONE}" | paste -d'\t' /dev/stdin <(samtools view -F2304 "${BAM_ONE}" | awk -v MIN_QUAL="$MIN_QUAL" '{if ($5>=MIN_QUAL) print $3,$4,$1; else print "*","*",$1}' OFS='\t') |  awk '{print $7,$8,$3,$4}' OFS='\t' > "${OUT_NAME}/${SAMPLE}_pos_r1.bed.tmp"
bedtools bamtobed -i "${BAM_TWO}" | paste -d'\t' /dev/stdin <(samtools view -F2304 "${BAM_TWO}" | awk -v MIN_QUAL="$MIN_QUAL" '{if ($5>=MIN_QUAL) print $3,$4,$1; else print "*","*",$1}' OFS='\t') |  awk '{print $7,$8,$3,$4}' OFS='\t' > "${OUT_NAME}/${SAMPLE}_pos_r2.bed.tmp"

if [ ! -f "${OUT_NAME}/${SAMPLE}_pos_r2.bed.tmp" ] ; then
    echo "File specified doesn't seem like a .bam and/or samtools is not accessible, aborting." | tee -a $LOG_FILE
    exit
fi

# The example files that I generated aren't actually paired, so I'm putting this condition as a warning
Total_PETs1=`samtools flagstat "${WK_DIR}/${BAM_ONE}" | head -1 | awk '{print $1}'`
Total_PETs2=`samtools flagstat "${WK_DIR}/${BAM_TWO}" | head -1 | awk '{print $1}'`
if [ "$Total_PETs1" -lt "$Total_PETs2" ]; then
	echo "`date`: WARNING-- .bam 1 has more reads than .bam 2. Check to see if they actually match!" | tee -a $LOG_FILE
	head "-${Total_PETs1}" "${OUT_NAME}/${SAMPLE}_pos_r2.bed.tmp" > temppets
	mv temppets "${OUT_NAME}/${SAMPLE}_pos_r2.bed.tmp"
elif [ "$Total_PETs2" -lt "$Total_PETs1" ]; then
	echo "`date`: WARNING-- .bam 2 has more reads than .bam 1. Check to see if they actually match!" | tee -a $LOG_FILE
	head "-${Total_PETs2}" "${OUT_NAME}/${SAMPLE}_pos_r1.bed.tmp" > temppets
	mv temppets "${OUT_NAME}/${SAMPLE}_pos_r1.bed.tmp"
fi

# Pair; keep reads that are both mapped
paste "${OUT_NAME}/${SAMPLE}_pos_r1.bed.tmp" "${OUT_NAME}/${SAMPLE}_pos_r2.bed.tmp" | cut -f 1-3,5-8 > "${OUT_NAME}/${SAMPLE}_rawinteractions.tmp" 
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
Mapped_unique_intrachromosomal_q30_Min_Dist=`awk -v MIN_DIST="$MIN_DIST" '$1 == $4 && (($5+$6) - (($2+$3)/2)>=MIN_DIST) {print $0}' "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" | wc -l | awk '{print $1}'`

echo "`date`: Mapped_unique_PETs_q30=${Mapped_unique_PETs_q30}" | tee -a $LOG_FILE
echo "`date`: Mapped_unique_intrachromosomal_q30=${Mapped_unique_intrachromosomal_q30}"| tee -a $LOG_FILE
echo "`date`: Mapped_unique_intrachromosomal_q30_Min_Dist=${Mapped_unique_intrachromosomal_q30_Min_Dist}" | tee -a $LOG_FILE

# Call peaks
cut -f 1-3,7 "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" >  "${OUT_NAME}/${SAMPLE}_left.dedup.bed.tmp"
cut -f 4-6,7 "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" >  "${OUT_NAME}/${SAMPLE}_right.dedup.bed.tmp"
cat "${OUT_NAME}/${SAMPLE}_left.dedup.bed.tmp" "${OUT_NAME}/${SAMPLE}_right.dedup.bed.tmp" >> "${OUT_NAME}/${SAMPLE}_reads.bed.tmp"	
macs2 callpeak -t "${OUT_NAME}/${SAMPLE}_reads.bed.tmp" -g hs -f BED -n "${OUT_NAME}/${SAMPLE}_temporary" $MACS2_STRING | tee -a $LOG_FILE

# Check to see if MACS2 worked
if [ ! -f "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak" ] ; then
    echo "MACS2 execution was not successful--modify command if possible" | tee -a $LOG_FILE
    exit
fi

# Pad peaks
awk -v PEAK_PAD="$PEAK_PAD" '{$2-=PEAK_PAD; $3+=PEAK_PAD}1' OFS="\t" "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak" > "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded"
echo "`date`: Padded peaks by ${PEAK_PAD}bp" | tee -a $LOG_FILE
awk '$2 > 0 && $5 > 0 {print $0}' "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded" > "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded.clean"
mv "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded.clean" "${OUT_NAME}/${SAMPLE}.peaks" 
NUM_PEAKS=`wc -l "${OUT_NAME}/${SAMPLE}.peaks" | awk '{print $1}'`

# Merge gaps
echo "`date`: Intersecting PETs with anchors" | tee -a $LOG_FILE
bedtools merge -d $MERGE_GAP -i "${OUT_NAME}/${SAMPLE}.peaks" > "${OUT_NAME}/${SAMPLE}_temporary_peaks.merged.bed.tmp"
bedtools intersect -loj -a "${OUT_NAME}/${SAMPLE}_left.dedup.bed.tmp" -b "${OUT_NAME}/${SAMPLE}_temporary_peaks.merged.bed.tmp" | awk '{print $5,$6,$7,$4}' OFS='\t' > "${OUT_NAME}/${SAMPLE}_anchor1.bed.tmp"
bedtools intersect -loj -a "${OUT_NAME}/${SAMPLE}_right.dedup.bed.tmp" -b "${OUT_NAME}/${SAMPLE}_temporary_peaks.merged.bed.tmp" | awk '{print $5,$6,$7,$4}' OFS='\t' > "${OUT_NAME}/${SAMPLE}_anchor2.bed.tmp"

# Make final output
paste "${OUT_NAME}/${SAMPLE}_anchor1.bed.tmp" "${OUT_NAME}/${SAMPLE}_anchor2.bed.tmp" | cut -f 1-3,5-8 > "${OUT_NAME}/${SAMPLE}_anchor.interactions.tmp"
awk '{if ($1 != "." && $4 != ".") print}' "${OUT_NAME}/${SAMPLE}_anchor.interactions.tmp" > "${OUT_NAME}/${SAMPLE}_anchor.interactions.bedpe.tmp"
cut -f1-6  "${OUT_NAME}/${SAMPLE}_anchor.interactions.bedpe.tmp" | sort | uniq -c | awk '{print $2,$3,$4,$5,$6,$7,".",$1}' >  "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp"
Mapped_unique_intrachromosomal_q30_5kb=`awk -v DIST=5000 '$1 == $4 && (($5+$6)/2 - (($2+$3)/2)>=DIST) {print $0}' "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" | wc -l | awk '{print $1}'`
echo "`date`: Mapped_unique_intrachromosomal_q30_5kb=${Mapped_unique_intrachromosomal_q30_5kb}"| tee -a $LOG_FILE

# Produce final output
awk '$1 != $4 {print $0}' "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}.inter.loop_counts.bedpe"
awk '$1 == $4 {print $0}' "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}.intra.loop_counts.bedpe"
awk -v MIN_DIST="$MIN_DIST" -v MAX_DIST="$MAX_DIST" '$1 == $4 && $2 != $5 && (($5+$6)/2 - ($2+$3)/2)>=MIN_DIST && (($5+$6)/2 - ($2+$3)/2)<=MAX_DIST {print $0}' "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}.filt.intra.loop_counts.bedpe"
Loop_PETs=`awk '{sum += $8} END {print sum}' "${OUT_NAME}/${SAMPLE}.filt.intra.loop_counts.bedpe"`
NUM_SHORT_PETs=`awk -v MIN_DIST="$MIN_DIST" '$1 == $4 && (($5+$6)/2 - ($2+$3)/2)<MIN_DIST {print $0}' "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp" | wc -l | awk '{print $1}'`
NUM_LONG_PETs=`awk  -v MAX_DIST="$MAX_DIST" '$1 == $4 && $2 != $5 && (($5+$6)/2 - ($2+$3)/2)>MAX_DIST {print $0}' "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp" | wc -l | awk '{print $1}'`
echo "`date`: Loop_PETs=${Loop_PETs}" | tee -a $LOG_FILE

# Write out summary statistics
echo "Total_PETs=${Total_PETs1}" > "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_PETs_q30=${Mapped_PETs_q30}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_PETs_q30=${Mapped_unique_PETs_q30}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intra_q30=${Mapped_unique_intrachromosomal_q30}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intra_q30_5kb=${Mapped_unique_intrachromosomal_q30_5kb}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intra_q30_loops=${Loop_PETs}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "NUM_SHORT_PETs=${NUM_SHORT_PETs}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "NUM_LONG_PETs=${NUM_LONG_PETs}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Number_of_Peaks=${NUM_PEAKS}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "MIN_LENGTH=${MIN_DIST}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "MAX_LENGTH=${MAX_DIST}" >> "${OUT_NAME}/${SAMPLE}.stat"