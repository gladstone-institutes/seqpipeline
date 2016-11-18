#!/bin/bash -u
set -e
set -o pipefail

### You can run this script with a command like one of these:
###                                           * bash ./this_script_name.sh
###             Or submit it to a cluster:    * qsub --some-options ./this_script_name.sh

echo 'First, we are checking to make sure all the FASTQ input files exist...'

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a2.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a2.mm9.chr19.fq.gz>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a3.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a3.mm9.chr19.fq.gz>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a1.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a1.mm9.chr19.fq.gz>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b1.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b1.mm9.chr19.fq.gz>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b3.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b3.mm9.chr19.fq.gz>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b2.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b2.mm9.chr19.fq.gz>'
    exit 49
fi


###### Handling RNA-seq here 
mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a1.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a1.mm9.chr19.fq.gz>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
tophat -o /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1 --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search  --GTF=/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf --num-threads=4   /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a1.mm9.chr19.fq.gz 
test -e /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam || mv -f /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a2.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a2.mm9.chr19.fq.gz>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
tophat -o /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2 --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search  --GTF=/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf --num-threads=4   /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a2.mm9.chr19.fq.gz 
test -e /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam || mv -f /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b1.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b1.mm9.chr19.fq.gz>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
tophat -o /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1 --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search  --GTF=/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf --num-threads=4   /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b1.mm9.chr19.fq.gz 
test -e /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam || mv -f /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b2.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b2.mm9.chr19.fq.gz>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
tophat -o /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2 --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search  --GTF=/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf --num-threads=4   /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b2.mm9.chr19.fq.gz 
test -e /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam || mv -f /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a3.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a3.mm9.chr19.fq.gz>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
tophat -o /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search  --GTF=/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf --num-threads=4   /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a3.mm9.chr19.fq.gz 
test -e /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam || mv -f /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b3.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b3.mm9.chr19.fq.gz>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
tophat -o /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search  --GTF=/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf --num-threads=4   /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b3.mm9.chr19.fq.gz 
test -e /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam || mv -f /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X1"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam.sorted_by_name.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
samtools sort -n /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam.sorted_by_name.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X1/htseq_count.Test_3_Compare.X1.matrix.txt" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
htseq-count  --format=bam --order=name --mode=intersection-nonempty --stranded=no --minaqual=10 --type=exon --idattr=gene_id  /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X1/accepted_hits.bam.sorted_by_name.bam /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X1/htseq_count.Test_3_Compare.X1.matrix.txt
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X1/htseq_count.Test_3_Compare.X1.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X1/htseq_count.Test_3_Compare.X1.matrix.txt>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X2"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam.sorted_by_name.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
samtools sort -n /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam.sorted_by_name.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X2/htseq_count.Test_3_Compare.X2.matrix.txt" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
htseq-count  --format=bam --order=name --mode=intersection-nonempty --stranded=no --minaqual=10 --type=exon --idattr=gene_id  /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/X2/accepted_hits.bam.sorted_by_name.bam /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X2/htseq_count.Test_3_Compare.X2.matrix.txt
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X2/htseq_count.Test_3_Compare.X2.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X2/htseq_count.Test_3_Compare.X2.matrix.txt>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y1"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam.sorted_by_name.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
samtools sort -n /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam.sorted_by_name.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y1/htseq_count.Test_3_Compare.Y1.matrix.txt" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
htseq-count  --format=bam --order=name --mode=intersection-nonempty --stranded=no --minaqual=10 --type=exon --idattr=gene_id  /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y1/accepted_hits.bam.sorted_by_name.bam /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y1/htseq_count.Test_3_Compare.Y1.matrix.txt
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y1/htseq_count.Test_3_Compare.Y1.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y1/htseq_count.Test_3_Compare.Y1.matrix.txt>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y2"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam.sorted_by_name.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
samtools sort -n /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam.sorted_by_name.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y2/htseq_count.Test_3_Compare.Y2.matrix.txt" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
htseq-count  --format=bam --order=name --mode=intersection-nonempty --stranded=no --minaqual=10 --type=exon --idattr=gene_id  /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Y2/accepted_hits.bam.sorted_by_name.bam /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y2/htseq_count.Test_3_Compare.Y2.matrix.txt
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y2/htseq_count.Test_3_Compare.Y2.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y2/htseq_count.Test_3_Compare.Y2.matrix.txt>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3A"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam.sorted_by_name.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
samtools sort -n /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam.sorted_by_name.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3A/htseq_count.Test_3_Compare.Z3A.matrix.txt" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
htseq-count  --format=bam --order=name --mode=intersection-nonempty --stranded=no --minaqual=10 --type=exon --idattr=gene_id  /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3A/accepted_hits.bam.sorted_by_name.bam /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3A/htseq_count.Test_3_Compare.Z3A.matrix.txt
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3A/htseq_count.Test_3_Compare.Z3A.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3A/htseq_count.Test_3_Compare.Z3A.matrix.txt>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3B"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam.sorted_by_name.bam" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
samtools sort -n /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam.sorted_by_name.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam.sorted_by_name.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam.sorted_by_name.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3B/htseq_count.Test_3_Compare.Z3B.matrix.txt" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
htseq-count  --format=bam --order=name --mode=intersection-nonempty --stranded=no --minaqual=10 --type=exon --idattr=gene_id  /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01t_Align_Tophat_Dir/Test_3_Compare/Z3B/accepted_hits.bam.sorted_by_name.bam /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.gtf > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3B/htseq_count.Test_3_Compare.Z3B.matrix.txt
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3B/htseq_count.Test_3_Compare.Z3B.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3B/htseq_count.Test_3_Compare.Z3B.matrix.txt>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_03_EdgeR_Diff_Expr_Dir"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X1/htseq_count.Test_3_Compare.X1.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X1/htseq_count.Test_3_Compare.X1.matrix.txt>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X2/htseq_count.Test_3_Compare.X2.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X2/htseq_count.Test_3_Compare.X2.matrix.txt>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y1/htseq_count.Test_3_Compare.Y1.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y1/htseq_count.Test_3_Compare.Y1.matrix.txt>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y2/htseq_count.Test_3_Compare.Y2.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y2/htseq_count.Test_3_Compare.Y2.matrix.txt>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3A/htseq_count.Test_3_Compare.Z3A.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3A/htseq_count.Test_3_Compare.Z3A.matrix.txt>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3B/htseq_count.Test_3_Compare.Z3B.matrix.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3B/htseq_count.Test_3_Compare.Z3B.matrix.txt>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_03_EdgeR_Diff_Expr_Dir/edgeR.Test_3_Compare.out.txt" ]]
then
    echo '[SKIPPING] All the specified output files already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
Rscript /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/run_edgeR_diff_expr.R --verbose  --groups="1,1,2,2,3,3" --countfiles="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X1/htseq_count.Test_3_Compare.X1.matrix.txt,/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/X2/htseq_count.Test_3_Compare.X2.matrix.txt,/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y1/htseq_count.Test_3_Compare.Y1.matrix.txt,/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Y2/htseq_count.Test_3_Compare.Y2.matrix.txt,/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3A/htseq_count.Test_3_Compare.Z3A.matrix.txt,/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_02_HTSeq_Count_Dir/Test_3_Compare/Z3B/htseq_count.Test_3_Compare.Z3B.matrix.txt" --outfile="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_03_EdgeR_Diff_Expr_Dir/edgeR.Test_3_Compare.out.txt"
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_03_EdgeR_Diff_Expr_Dir/edgeR.Test_3_Compare.out.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_03_EdgeR_Diff_Expr_Dir/edgeR.Test_3_Compare.out.txt>'
    exit 49
fi

