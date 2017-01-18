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

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/ALPHA"

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


# ===== Checking for a required file
if [[ ! -e "/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.chrom.sizes.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.chrom.sizes.txt>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr_chromosome_fastas" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr_chromosome_fastas>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/A-2016-08-August/github_seqpipeline/resources/GEM_Read_Distribution_default.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/A-2016-08-August/github_seqpipeline/resources/GEM_Read_Distribution_default.txt>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/ALPHA/ALPHA_result.htm" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/ALPHA/ALPHA_result.htm) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
java -Xmx25G -jar /data/applications/2015_06/bin/gem.jar --t 1 --q 2 --k_min 6 --k_max 13 --d /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/A-2016-08-August/github_seqpipeline/resources/GEM_Read_Distribution_default.txt --g /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.chrom.sizes.txt --genome /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr_chromosome_fastas --expt /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA/accepted_hits.bam --ctrl /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/ALPHA_CHIP_CONTROL/accepted_hits.bam --f SAM  --sl --outBED  --out /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/ALPHA
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/ALPHA/ALPHA_result.htm" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/ALPHA/ALPHA_result.htm>'
    exit 49
fi

mkdir -p "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/BETA"

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


# ===== Checking for a required file
if [[ ! -e "/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.chrom.sizes.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.chrom.sizes.txt>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr_chromosome_fastas" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr_chromosome_fastas>'
    exit 49
fi


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/A-2016-08-August/github_seqpipeline/resources/GEM_Read_Distribution_default.txt" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/A-2016-08-August/github_seqpipeline/resources/GEM_Read_Distribution_default.txt>'
    exit 49
fi


#-----------------------
if [[ -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/BETA/BETA_result.htm" ]]
then
    echo '[SKIPPING] All the specified output files (/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/BETA/BETA_result.htm) already exist--so we are skipping this step' 
else
    echo '[RUNNING]...'
java -Xmx25G -jar /data/applications/2015_06/bin/gem.jar --t 1 --q 2 --k_min 6 --k_max 13 --d /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/A-2016-08-August/github_seqpipeline/resources/GEM_Read_Distribution_default.txt --g /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr.chrom.sizes.txt --genome /data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr_chromosome_fastas --expt /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam --ctrl /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam --f SAM  --sl --outBED  --out /data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/BETA
fi
#-----------------------


# ===== Checking for a required file
if [[ ! -e "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/BETA/BETA_result.htm" ]]; then
    echo '[ERROR] Cannot find the following required file, which we expect to exist: </data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/Pipeline_05g_GEM_Peak_Dir/Test_Experiment_ChIP/BETA/BETA_result.htm>'
    exit 49
fi

