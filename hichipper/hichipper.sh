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
MIN_DIST=${9}
MAX_DIST=${10}
MACS2_STRING=${11}

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
samtools view "${BAM_ONE}" | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | paste -d'\t' /dev/stdin <(samtools view -F2304 "${BAM_ONE}" | awk -v MIN_QUAL="$MIN_QUAL" '{if ($5>=MIN_QUAL) print $3,$4,$1; else print "*","*",$1}' OFS='\t') | awk '{print $2"\t"$3"\t"$3+$1"\t."}' > "${OUT_NAME}/${SAMPLE}_pos_r1.bed.tmp"
samtools view "${BAM_TWO}" | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | paste -d'\t' /dev/stdin <(samtools view -F2304 "${BAM_TWO}" | awk -v MIN_QUAL="$MIN_QUAL" '{if ($5>=MIN_QUAL) print $3,$4,$1; else print "*","*",$1}' OFS='\t') | awk '{print $2"\t"$3"\t"$3+$1"\t."}' > "${OUT_NAME}/${SAMPLE}_pos_r2.bed.tmp"

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

# Pair; keep reads that are both mapped then swap anchors if necessary
paste "${OUT_NAME}/${SAMPLE}_pos_r1.bed.tmp" "${OUT_NAME}/${SAMPLE}_pos_r2.bed.tmp" | cut -f 1-3,5-8 | awk '{if ($1 != "*" && $4 != "*") print}' > "${OUT_NAME}/${SAMPLE}_interactions.unsorted.bedpe.tmp"
awk '{if ($1<$4 || ($1==$4 && $2<=$5)) print $0; else print $4,$5,$6,$1,$2,$3,$7}' 'OFS=\t' "${OUT_NAME}/${SAMPLE}_interactions.unsorted.bedpe.tmp" | sort -k1,1n -k2,2n -k4,4n -k5,5n  >  "${OUT_NAME}/${SAMPLE}_interactions.bedpe.tmp"

echo "`date`: Writing unique (i.e. deduplicated) interactions (based on chr1, pos1, chr2, pos2)" | tee -a $LOG_FILE
sort -k1,1n -k2,2n -k4,4n -k5,5n --unique "${OUT_NAME}/${SAMPLE}_interactions.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp"

# Summary stats
Mapped_PETs_q30=`wc -l "${OUT_NAME}/${SAMPLE}_interactions.bedpe.tmp" | awk '{print $1}'`
Mapped_unique_PETs_q30=`wc -l "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" | awk '{print $1}'`
echo "`date`: Mapped_unique_PETs_q30=${Mapped_unique_PETs_q30}" | tee -a $LOG_FILE
Mapped_unique_intrachromosomal_q30=`awk '$1 == $4 {print $0}' "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" | wc -l | awk '{print $1}'`
echo "`date`: Mapped_unique_intrachromosomal_q30=${Mapped_unique_intrachromosomal_q30}"| tee -a $LOG_FILE

Mapped_unique_intra_q30_small=`awk -v MIN_DIST="$MIN_DIST" '$1 == $4 && (($5+$6)/2 - ($2+$3)/2)<=MIN_DIST {print $0}' "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" | wc -l | awk '{print $1}'`
Mapped_unique_intra_q30_med=`awk -v MIN_DIST="$MIN_DIST" -v MAX_DIST="$MAX_DIST" '$1 == $4 && (($5+$6)/2 - ($2+$3)/2)>=MIN_DIST && (($5+$6)/2 - ($2+$3)/2)<=MAX_DIST {print $0}' "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" | wc -l | awk '{print $1}'`
Mapped_unique_intra_q30_large=`awk -v MAX_DIST="$MAX_DIST" '$1 == $4 && (($5+$6)/2 - ($2+$3)/2)>=MAX_DIST {print $0}' "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" | wc -l | awk '{print $1}'`
echo "`date`: Mapped_unique_intra_q30_small=${Mapped_unique_intra_q30_small}" | tee -a $LOG_FILE
echo "`date`: Mapped_unique_intra_q30_med=${Mapped_unique_intra_q30_med}" | tee -a $LOG_FILE
echo "`date`: Mapped_unique_intra_q30_large=${Mapped_unique_intra_q30_large}" | tee -a $LOG_FILE

# Split
cut -f 1-3,7 "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" >  "${OUT_NAME}/${SAMPLE}_left.dedup.bed.tmp"
cut -f 4-6,7 "${OUT_NAME}/${SAMPLE}_interactions.dedup.bedpe.tmp" >  "${OUT_NAME}/${SAMPLE}_right.dedup.bed.tmp"

# Call peaks
cat "${OUT_NAME}/${SAMPLE}_left.dedup.bed.tmp" "${OUT_NAME}/${SAMPLE}_right.dedup.bed.tmp" >> "${OUT_NAME}/${SAMPLE}_reads.bed.tmp"	
macs2 callpeak -t "${OUT_NAME}/${SAMPLE}_reads.bed.tmp" -g hs -f BED -n "${OUT_NAME}/${SAMPLE}_temporary" $MACS2_STRING | tee -a $LOG_FILE

# Check to see if MACS2 worked
if [ ! -f "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak" ] ; then
    echo "MACS2 execution was not successful--modify command if possible" | tee -a $LOG_FILE
    exit
fi

# Pad peaks
awk -v PEAK_PAD="$PEAK_PAD" '{$2-=PEAK_PAD; $3+=PEAK_PAD}1' OFS="\t" "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak" > "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded"
echo "`date`: Padding peaks by ${PEAK_PAD}bp" | tee -a $LOG_FILE
awk '$2 > 0 && $5 > 0 {print $0}' "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded" > "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded.clean"
mv "${OUT_NAME}/${SAMPLE}_temporary_peaks.narrowPeak.padded.clean" "${OUT_NAME}/${SAMPLE}.peaks" 
NUM_PEAKS=`wc -l "${OUT_NAME}/${SAMPLE}.peaks" | awk '{print $1}'`
echo "`date`: Total number of peaks called: ${NUM_PEAKS}" | tee -a $LOG_FILE

# Merge gaps; check bedtools
echo "`date`: Intersecting PETs with anchors" | tee -a $LOG_FILE
bedtools merge -d $MERGE_GAP -i "${OUT_NAME}/${SAMPLE}.peaks" > "${OUT_NAME}/${SAMPLE}_temporary_peaks.merged.bed.tmp"
if [ ! -f  "${OUT_NAME}/${SAMPLE}_temporary_peaks.merged.bed.tmp" ] ; then
    echo "Bedtools does not appear" | tee -a $LOG_FILE
    exit
fi

bedtools intersect -loj -a "${OUT_NAME}/${SAMPLE}_left.dedup.bed.tmp" -b "${OUT_NAME}/${SAMPLE}_temporary_peaks.merged.bed.tmp" | awk '{print $5,$6,$7,$4}' OFS='\t' > "${OUT_NAME}/${SAMPLE}_anchor1.bed.tmp"
bedtools intersect -loj -a "${OUT_NAME}/${SAMPLE}_right.dedup.bed.tmp" -b "${OUT_NAME}/${SAMPLE}_temporary_peaks.merged.bed.tmp" | awk '{print $5,$6,$7,$4}' OFS='\t' > "${OUT_NAME}/${SAMPLE}_anchor2.bed.tmp"

# Getting close
paste "${OUT_NAME}/${SAMPLE}_anchor1.bed.tmp" "${OUT_NAME}/${SAMPLE}_anchor2.bed.tmp" | cut -f 1-3,5-8 | awk '{if ($1 != "." && $4 != ".") print}'  > "${OUT_NAME}/${SAMPLE}_anchor.interactions.bedpe.tmp"
cut -f1-6  "${OUT_NAME}/${SAMPLE}_anchor.interactions.bedpe.tmp" | sort | uniq -c | awk '{print $2,$3,$4,$5,$6,$7,".",$1}' >  "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp"
Mapped_unique_intra_q30_anchor=`awk '$1 == $4 {print $0}' "${OUT_NAME}/${SAMPLE}_anchor.interactions.bedpe.tmp" | wc -l | awk '{print $1}'`
Mapped_unique_intra_q30_anchor_small=`awk -v MIN_DIST="$MIN_DIST" '$1 == $4 && (($5+$6)/2 - ($2+$3)/2)<=MIN_DIST {print $0}' "${OUT_NAME}/${SAMPLE}_anchor.interactions.bedpe.tmp" | wc -l | awk '{print $1}'`
Mapped_unique_intra_q30_anchor_med=`awk -v MIN_DIST="$MIN_DIST" -v MAX_DIST="$MAX_DIST" '$1 == $4 && (($5+$6)/2 - ($2+$3)/2)>=MIN_DIST && (($5+$6)/2 - ($2+$3)/2)<=MAX_DIST {print $0}' "${OUT_NAME}/${SAMPLE}_anchor.interactions.bedpe.tmp" | wc -l | awk '{print $1}'`
Mapped_unique_intra_q30_anchor_large=`awk -v MAX_DIST="$MAX_DIST" '$1 == $4 && (($5+$6)/2 - ($2+$3)/2)>=MAX_DIST {print $0}' "${OUT_NAME}/${SAMPLE}_anchor.interactions.bedpe.tmp" | wc -l | awk '{print $1}'`
echo "`date`: Mapped_unique_intra_q30_anchor=${Mapped_unique_intra_q30_anchor}" | tee -a $LOG_FILE
echo "`date`: Mapped_unique_intra_q30_anchor_small=${Mapped_unique_intra_q30_anchor_small}" | tee -a $LOG_FILE
echo "`date`: Mapped_unique_intra_q30_anchor_med=${Mapped_unique_intra_q30_anchor_med}" | tee -a $LOG_FILE
echo "`date`: Mapped_unique_intra_q30_anchor_large=${Mapped_unique_intra_q30_anchor_large}" | tee -a $LOG_FILE

# Produce final output
awk '$1 != $4 {print $0}' "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}.inter.loop_counts.bedpe"
awk '$1 == $4 {print $0}' "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}.intra.loop_counts.bedpe"
awk -v MIN_DIST="$MIN_DIST" -v MAX_DIST="$MAX_DIST" '$1 == $4 && $2 != $5 && (($5+$6)/2 - ($2+$3)/2)>=MIN_DIST && (($5+$6)/2 - ($2+$3)/2)<=MAX_DIST {print $0}' "${OUT_NAME}/${SAMPLE}.loop_counts.bedpe.tmp" > "${OUT_NAME}/${SAMPLE}.filt.intra.loop_counts.bedpe"

# Finalize 
Loop_PETs=`awk '{sum += $8} END {print sum}' "${OUT_NAME}/${SAMPLE}.filt.intra.loop_counts.bedpe"`
echo "`date`: Loop_PETs=${Loop_PETs}" | tee -a $LOG_FILE

# Write out summary statistics
echo "Total_PETs=${Total_PETs1}" > "${OUT_NAME}/${SAMPLE}.stat"

echo "Mapped_PETs_q30=${Mapped_PETs_q30}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_PETs_q30=${Mapped_unique_PETs_q30}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intra_q30=${Mapped_unique_intrachromosomal_q30}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intra_q30_small=${Mapped_unique_intra_q30_small}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intra_q30_med=${Mapped_unique_intra_q30_med}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intra_q30_large=${Mapped_unique_intra_q30_large}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intra_q30_anchor=${Mapped_unique_intra_q30_anchor}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intra_q30_anchor_small=${Mapped_unique_intra_q30_anchor_small}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intra_q30_anchor_med=${Mapped_unique_intra_q30_anchor_med}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Mapped_unique_intra_q30_anchor_large=${Mapped_unique_intra_q30_anchor_large}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "Number_of_Peaks=${NUM_PEAKS}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "MIN_LENGTH=${MIN_DIST}" >> "${OUT_NAME}/${SAMPLE}.stat"
echo "MAX_LENGTH=${MAX_DIST}" >> "${OUT_NAME}/${SAMPLE}.stat"