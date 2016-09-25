# seqpipeline

This is bioinformatics sequence-analysis pipelining software.

Currently, it is written entirely in R.

This is located at https://github.com/gladstone-institutes/seqpipeline


Required software:
Tophat (splice-aware aligner)
     Executable name: tophat (version 2.1.1)
          To install: See details at: http://ccb.jhu.edu/software/tophat/index.shtml

Bowtie (non-spliced aligner)
     Executable name: bowtie2 (version 2.2.4)
          To install: See details at: http://bowtie-bio.sourceforge.net/index.shtml

BCP (ChIP-seq peak caller for everything except TF binding)
     Executable name: BCP_HM (version 1.1)
          To install: Available at https://cb.utdallas.edu/BCP/
                      (Paper: http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1002613)

GEM (motif-aware ChIP-seq peak caller for TF binding sites)
     Executable name: gem.jar (version 2.5)
          To install: See instructions at: http://groups.csail.mit.edu/cgs/gem/
                      (We were using version 2.5, but version 2.7+ is available now.)

htseq-count (reads -> gene-level counts)
     Executable name: htseq-count (version 0.6.0)
          To install: pip2.7 install HTseq
                      (Note: must be installed via pip (or other package manager). Do not just copy the binaries--it will not work.)

samtools
     Executable name: samtools (version 1.3)
          To install: yum install samtools (on Redhat/CentOS)
                      (Available through your package manager. Other examples: brew install samtools (Mac Homebrew), apt-get install samtools (Ubuntu))

edgeR (R library for differential expression)
     Executable name: NA (version 3.14.0)
          To install: source('https://bioconductor.org/biocLite.R'); biocLite('edgeR')
                      (Available through Bioconductor.)

java
     Executable name: java (version 1.8.0)
          To install: Install via package manager or on Oracle's web site: https://java.com/en/
                      (May also be possible to install using your package manager (e.g. 'yum install java'). Java is required to run 'gem.jar' and 'MarkDuplicates.jar')
