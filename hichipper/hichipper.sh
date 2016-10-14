#!/bin/bash

MIN_QUAL=30
READ_LEN=75
LOG_FILE=HiChIPPer.log
PEAK_PAD=1500
MERGE_GAP=1500

# Parse command line samples
if [ "$#" -eq 0 ]; then
  { echo "ERROR-- Must specify at least one sample"; exit 1; }
else
  myArray=("$@")
fi

# Define function
processEverything () {
    file=$1
    echo "`date`: Processing ${file}"
    
    # Process .bam into bedpe
    samtools view -F2304 "${file}_1_hg19.bwt2merged.bam" | awk -v READ_LEN="$READ_LEN" -v MIN_QUAL="$MIN_QUAL" '{if ($5>=MIN_QUAL) print $3,$4,$4+READ_LEN,$1; else print "*","*","*",$1}' OFS='\t' > "${file}_pos_r1.bed"
	samtools view -F2304 "${file}_2_hg19.bwt2merged.bam" | awk -v READ_LEN="$READ_LEN" -v MIN_QUAL="$MIN_QUAL" '{if ($5>=MIN_QUAL) print $3,$4,$4+READ_LEN,$1; else print "*","*","*",$1}' OFS='\t' > "${file}_pos_r2.bed"
	paste "${file}_pos_r1.bed" "${file}_pos_r2.bed" | cut -f 1-3,5-8 > "${file}_interactions.tmp" 
	Total_PETs=`samtools flagstat "${file}_1_hg19.bwt2merged.bam" | head -1 | awk '{print $1}'`

	# Keep only interactions where both reads are mapped
	awk '{if ($1 != "*" && $4 != "*") print}' "${file}_interactions.tmp" > "${file}_interactions.unsorted.bedpe"
	
	# Swap anchors if necessary
	awk '{if ($1<$4 || ($1==$4 && $2<=$5)) print $0; else print $4,$5,$6,$1,$2,$3,$7}' 'OFS=\t' "${file}_interactions.unsorted.bedpe" > "${file}_interactions.unsorted.orderedanchors.bedpe"
	sort -k1,1V -k2,2n -k4,4V -k5,5n "${file}_interactions.unsorted.orderedanchors.bedpe" > "${file}_interactions.bedpe"
	Mapped_PETs_q30=`wc -l "${file}_interactions.bedpe" | awk '{print $1}'`
	
	# Remove duplicates
	echo "`date`: Writing unique (i.e. deduplicated) interactions (based on chr1, pos1, chr2, pos2)" | tee -a $LOG_FILE
	sort -k1,1V -k2,2n -k4,4V -k5,5n --unique "${file}_interactions.bedpe" > "${file}_interactions.dedup.bedpe"
	Mapped_unique_PETs_q30=`wc -l "${file}_interactions.dedup.bedpe" | awk '{print $1}'`
	Mapped_unique_intrachromosomal_q30=`awk '$1 == $4 {print $0}' "${file}_interactions.dedup.bedpe" | wc -l | awk '{print $1}'`
	Mapped_unique_intrachromosomal_q30_5kb=`awk '$1 == $4 && (($5+$6) - (($2+$3)/2)> 5000) {print $0}' "${file}_interactions.dedup.bedpe" | wc -l | awk '{print $1}'`
	
	# Call peaks
	cut -f 1-3,7 "${file}_interactions.dedup.bedpe" >  "${file}_left.dedup.bed"
	cut -f 4-6,7 "${file}_interactions.dedup.bedpe" >  "${file}_right.dedup.bed"
	cat "${file}_left.dedup.bed" "${file}_right.dedup.bed" >> "${file}_reads.bed"	
	macs2 callpeak -t "${file}_reads.bed" -g hs -f BED -n $file -p 0.01 --nomodel
	
	# Pad peaks
	awk -v PEAK_PAD="$PEAK_PAD" '{$2-=PEAK_PAD; $3+=PEAK_PAD}1' OFS="\t" "${file}_peaks.narrowPeak" > "${file}_peaks.narrowPeak.padded"
	echo "`date`: Padded peaks by ${PEAK_PAD}bp" | tee -a $LOG_FILE
	awk '$2 > 0 && $5 > 0 {print $0}' "${file}_peaks.narrowPeak.padded" > "${file}_peaks.narrowPeak.padded.clean"
	
	# Merge gaps
	bedtools merge -d $MERGE_GAP -i "${file}_peaks.narrowPeak.padded.clean" > "${file}_peaks.merged.bed"
	bedtools intersect -loj -a "${file}_left.dedup.bed" -b "${file}_peaks.merged.bed" | awk '{print $5,$6,$7,$4}' OFS='\t' > "${file}_anchor1.bed"
	bedtools intersect -loj -a "${file}_right.dedup.bed" -b "${file}_peaks.merged.bed" | awk '{print $5,$6,$7,$4}' OFS='\t' > "${file}_anchor2.bed"
	
	# Make final output
	paste "${file}_anchor1.bed" "${file}_anchor2.bed" | cut -f 1-3,5-8 > "${file}_anchor.interactions.tmp"
	awk '{if ($1 != "." && $4 != ".") print}' "${file}_anchor.interactions.tmp" > "${file}_anchor.interactions.bedpe"
	cut -f1-6  "${file}_anchor.interactions.bedpe" | sort | uniq -c | awk '{print $2,$3,$4,$5,$6,$7,".",$1}' >  "${file}.loop_counts.bedpe"
	awk '$1 != $4 {print $0}' "${file}.loop_counts.bedpe" > "${file}.inter.loop_counts.bedpe"
	awk '$1 == $4 && $2 != $5 {print $0}' "${file}.loop_counts.bedpe" > "${file}.intra.loop_counts.bedpe"
	
	# Write out summary statistics
	echo "Total_PETs=${Total_PETs}" > "${file}.stat"
	echo "Mapped_PETs_q30=${Mapped_PETs_q30}" >> "${file}.stat"
	echo "Mapped_unique_PETs_q30=${Mapped_unique_PETs_q30}" >> "${file}.stat"
	echo "Mapped_unique_intrachromosomal_q30=${Mapped_unique_intrachromosomal_q30}" >> "${file}.stat"
	echo "Mapped_unique_intrachromosomal_q30_5kb=${Mapped_unique_intrachromosomal_q30_5kb}" >> "${file}.stat"

}

# Parallel Execution?
for sample in "${myArray[@]}"; do
	processEverything "$sample"
done
