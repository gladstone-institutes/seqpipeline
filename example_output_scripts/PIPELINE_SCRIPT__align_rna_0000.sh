#!/bin/bash -u
set -e
set -o pipefail
if [[ -e "./Pipeline_01_Spliced_Tophat_Alignment/accepted_hits.bam" ]]
then
    echo '[OK] All the specified output files already exist--so we are exiting the script now.' 
    exit 0
fi
if [[ ! -e "/data/a1.fq.gz" ]]; then
      echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/a1.fq.gz>
      exit 49
fi
mkdir -p "./Pipeline_01_Spliced_Tophat_Alignment"
tophat -o ./Pipeline_01_Spliced_Tophat_Alignment --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search  --GTF=___PLACEHOLDER_FOR_ANNOTATION_FILE___ --num-threads=4   ___PLACEHOLDER_FOR_ANNOTATION_FILE___ /data/a1.fq.gz 
if [[ ! -e "./Pipeline_01_Spliced_Tophat_Alignment/accepted_hits.bam" ]]; then
      echo '[ERROR] Cannot find the following required file, which we expect to exist: <./Pipeline_01_Spliced_Tophat_Alignment/accepted_hits.bam>
      exit 49
fi
