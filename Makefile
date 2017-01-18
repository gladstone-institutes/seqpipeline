

all: test_3groups


test_2groups:
	echo "Running a test with 2 groups and 4 samples (2 replicates in A, 2 replicates in B)"
	python2  ./pipeline.py --basedir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/" \
                --outdir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/" \
		--experiment-id="Test_Experiment" \
		--sample-ids="A1,A2,B1,B2" \
		--rna-samples=a1.mm9.chr19.fq.gz,a2.mm9.chr19.fq.gz,b1.mm9.chr19.fq.gz,b2.mm9.chr19.fq.gz \
		--groups=1,1,2,2 --species=mm9  --script="script_test.sh"


test_3groups:
	echo "Running a test with 3 groups and 6 samples (2 replicates in 'X', 2 replicates in 'Y', 2 replicates in 'Z')"
	python2  ./pipeline.py --basedir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/" \
                --outdir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/" \
		--experiment-id="Test_3_Compare" \
		--sample-ids="X1,X2,Y1,Y2,Z3A,Z3B" \
		--groups="1,1,2,2,3,3" \
		--rna-samples=a1.mm9.chr19.fq.gz,a2.mm9.chr19.fq.gz,b1.mm9.chr19.fq.gz,b2.mm9.chr19.fq.gz,a3.mm9.chr19.fq.gz,b3.mm9.chr19.fq.gz \
		--species=mm9  --script="script_3_compare_test.sh"

test_gem:
	python2  ./pipeline.py \
--script="$@.sh" \
 --outdir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/" \
--basedir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/" \
--experiment-id="Test_Experiment_ChIP" \
--sample-ids="ALPHA,BETA" \
--groups=1,2 \
--peak-caller=GEM \
--chip-samples=a1.mm9.chr19.fq.gz,b1.mm9.chr19.fq.gz \
 --chip-inputs=a2.mm9.chr19.fq.gz,b2.mm9.chr19.fq.gz \
--species=mm9


test_bcp:
	python2  ./pipeline.py \
--script="$@.sh" \
 --outdir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/" \
--basedir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/" \
--experiment-id="Test_Experiment_ChIP" \
--sample-ids="ALPHA,BETA" \
--groups=1,2 \
--peak-caller=BCP \
--chip-samples=a1.mm9.chr19.fq.gz,b1.mm9.chr19.fq.gz \
 --chip-inputs=a2.mm9.chr19.fq.gz,b2.mm9.chr19.fq.gz \
--species=mm9


# Test BCP:
#  BCP_HM -1 Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA/accepted_hits.bam.converted_to.bed -2 Pipeline_01b_Align_Bowtie_Dir/Test_Experiment_ChIP/BETA_CHIP_CONTROL/accepted_hits.bam.converted_to.bed -f 2 -w 1 -p 1.1 -3 out.bed
