#!/bin/bash -u
set -e
set -o pipefail

### You can run this script with a command like one of these:
###                                           * bash ./this_script_name.sh
###             Or submit it to a cluster:    * qsub --some-options ./this_script_name.sh

echo 'First, we are checking to make sure all the FASTQ input files exist...'

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a1.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a1.mm9.chr19.fq.gz>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a2.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a2.mm9.chr19.fq.gz>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b1.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b1.mm9.chr19.fq.gz>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b2.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b2.mm9.chr19.fq.gz>'
    exit 49
fi


###### Handling ChIP-seq here 
mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a1.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a1.mm9.chr19.fq.gz>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
bowtie2 --threads=4 -x /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr -U "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a1.mm9.chr19.fq.gz"  | samtools view -@ 4 -b -  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/all.bam && samtools view -@ 4 -b -F 0x104 /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/all.bam  | samtools sort -  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam &&  /bin/rm /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/all.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a2.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a2.mm9.chr19.fq.gz>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
bowtie2 --threads=4 -x /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr -U "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/a2.mm9.chr19.fq.gz"  | samtools view -@ 4 -b -  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/all.bam && samtools view -@ 4 -b -F 0x104 /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/all.bam  | samtools sort -  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam &&  /bin/rm /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/all.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b1.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b1.mm9.chr19.fq.gz>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
bowtie2 --threads=4 -x /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr -U "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b1.mm9.chr19.fq.gz"  | samtools view -@ 4 -b -  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/all.bam && samtools view -@ 4 -b -F 0x104 /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/all.bam  | samtools sort -  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam &&  /bin/rm /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/all.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b2.mm9.chr19.fq.gz" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b2.mm9.chr19.fq.gz>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
bowtie2 --threads=4 -x /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr -U "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/b2.mm9.chr19.fq.gz"  | samtools view -@ 4 -b -  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/all.bam && samtools view -@ 4 -b -F 0x104 /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/all.bam  | samtools sort -  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam &&  /bin/rm /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/all.bam
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/ALPHA"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam.converted_to.bed" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam.converted_to.bed) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
bam2bed --do-not-sort  < /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam | cut -f 1-3  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam.converted_to.bed
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam.converted_to.bed" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam.converted_to.bed>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
bam2bed --do-not-sort  < /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam | cut -f 1-3  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/ALPHA/ALPHA_BCP_peaks.bed" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/ALPHA/ALPHA_BCP_peaks.bed) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
BCP_HM -1 /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam.converted_to.bed -2 /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed -f 200 -w 200 -p 0.001 -3 /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/ALPHA/ALPHA_BCP_peaks.bed
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/ALPHA/ALPHA_BCP_peaks.bed" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/ALPHA/ALPHA_BCP_peaks.bed>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/BETA"

# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam.converted_to.bed" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam.converted_to.bed) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
bam2bed --do-not-sort  < /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam | cut -f 1-3  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam.converted_to.bed
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam.converted_to.bed" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam.converted_to.bed>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
bam2bed --do-not-sort  < /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam | cut -f 1-3  > /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/BETA/BETA_BCP_peaks.bed" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/BETA/BETA_BCP_peaks.bed) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
BCP_HM -1 /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam.converted_to.bed -2 /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed -f 200 -w 200 -p 0.001 -3 /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/BETA/BETA_BCP_peaks.bed
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/BETA/BETA_BCP_peaks.bed" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05b_BCP_Peak_Dir/Test_Experiment_ChIP/BETA/BETA_BCP_peaks.bed>'
    exit 49
fi

