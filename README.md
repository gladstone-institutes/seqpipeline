# seqpipeline (pipeline.py)

This is bioinformatics sequence-analysis pipelining software.

Currently, it is written in python, with a single R script to handle 'edgeR' differential expression.

The code is located at https://github.com/gladstone-institutes/seqpipeline

# How to run it:

1. Make a new directory where you're going to run everything.
2. Download the test data (FASTQ reads) from http://gb.ucsf.edu/bio/public/kp-600/test_data/
      * There are six of these files.
3. Put that test FASTQ data into a new folder named 'test_data'
4. You can now run one of the commands in the 'Makefile' for this project, which at the moment has two options:
     * test_2groups
     * test_3groups
5. If you run "make test_2groups", pipeline.py will be run and will generate an output file named 'script_test.sh'
6. You can then invoke that script by running "bash script_test.sh". That is how you actually run all the bioinformatics tools.

# Example command:

    python2  ./pipeline.py --basedir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/" \
                --outdir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/" \
		--experiment-id="Test_3_Compare" \
		--sample-ids="X1,X2,Y1,Y2,Z3A,Z3B" \
		--groups="1,1,2,2,3,3" \
		--rna-samples=a1.mm9.chr19.fq.gz,a2.mm9.chr19.fq.gz,b1.mm9.chr19.fq.gz,b2.mm9.chr19.fq.gz,a3.mm9.chr19.fq.gz,b3.mm9.chr19.fq.gz \
		--species=mm9  --script="script_3_compare_test.sh"


# To-do:

Currently, only RNA-seq has been properly debugged in the updated 'pipeline.py' program.

=========================================

# Required software:
Tophat (splice-aware aligner)
*     Executable name: tophat (version 2.1.1)
     *     To install: See details at: http://ccb.jhu.edu/software/tophat/index.shtml

Bowtie (non-spliced aligner)
*     Executable name: bowtie2 (version 2.2.4)
     *     To install: See details at: http://bowtie-bio.sourceforge.net/index.shtml

BCP (ChIP-seq peak caller for everything except TF binding)
*     Executable name: BCP_HM (version 1.1)
     *     To install: Available at https://cb.utdallas.edu/BCP/
                      (Paper: http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1002613 )

bam2bed (Part of the 'bedops' suite. Converts BAM regions to BED files.)
*     Executable name: bam2bed (version 2.4.20)
     *     Note: This program is required only for running BCP in the CHIPseq pipeline--nothing else uses it.
     *     To install: Available here: http://bedops.readthedocs.io/en/latest/content/installation.html
                      (Paper: http://bioinformatics.oxfordjournals.org/content/28/14/1919.abstract )

GEM (motif-aware ChIP-seq peak caller for TF binding sites)
*     Executable name: gem.jar (version 2.5)
     *     To install: See instructions at: http://groups.csail.mit.edu/cgs/gem/
     *                (We were using version 2.5, but version 2.7+ is available now.)

htseq-count (reads -> gene-level counts)
*     Executable name: htseq-count (version 0.6.0)
     *     To install: pip2.7 install HTseq
     *                (Note: must be installed via pip (or other package manager). Do not just copy the binaries--it will not work.)

samtools
*     Executable name: samtools (version 1.3)
     *     To install: yum install samtools (on Redhat/CentOS)
     *                (Available through your package manager. Other examples: brew install samtools (Mac Homebrew), apt-get install samtools (Ubuntu))

edgeR (R library for differential expression)
*     Executable name: NA (version 3.14.0)
     *     To install: source('https://bioconductor.org/biocLite.R'); biocLite('edgeR')
                      (Available through Bioconductor.)

java
*     Executable name: java (version 1.8.0)
     *     To install: Install via package manager or on Oracle's web site: https://java.com/en/
                      (May also be possible to install using your package manager (e.g. 'yum install java'). Java is required to run 'gem.jar' and 'MarkDuplicates.jar')
