#!/bin/bash

############################################
# do this before running dada2.R script!
############################################

# adjust default fastq names to sid_R[12].fastq.gz
rename 's/_S.*\.fastq\.gz/_R1\.fastq\.gz/' *_R1*
rename 's/_S.*\.fastq\.gz/_R2\.fastq\.gz/' *_R2*

# print command strings to text file
outpref="fastp-r-W5-M20-l150-Q-p100-D1__"
for r1 in $(ls --color=never *_R1*); do
	r2=`echo $r1 | sed 's/R1/R2/'`
	sid=`echo $r1 | sed 's/_R1.*gz//'`
	echo fastp -i $r1 -I $r2 \
		-o $outpref$r1 -O $outpref$r2 \
		--unpaired1 ${outpref}UnP_$r1 --unpaired2 ${outpref}UnP_$r2 \
		-R "${sid} fastp report" -j ${sid}.fastp.json -h ${sid}.fastp.html \
		--dont_overwrite \
		-Q -p -P 100 -D 1 \
		-r -W 5 -M 20 \
		-l 150 >> cmdToParallel.txt
done

# run all commands in parallel
cat cmdToParallel.txt | parallel

# cleanup directory
DIRECTORY="../fastp/fastp.fq"
if [ ! -d "$DIRECTORY" ]; then
  mkdir -p "$DIRECTORY"
fi
mv fastp*.fastq.gz "$DIRECTORY"

DIRECTORY="../fastp/fastp.fq/UnP"
if [ ! -d "$DIRECTORY" ]; then
  mkdir "$DIRECTORY"
fi
mv ../fastp/fastp.fq/*__UnP_*.fastq.gz "$DIRECTORY"

DIRECTORY="../fastp/fastp.reports"
if [ ! -d "$DIRECTORY" ]; then
  mkdir "$DIRECTORY"
fi
mv *.fastp* "$DIRECTORY"