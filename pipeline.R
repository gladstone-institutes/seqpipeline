#!/usr/bin/env Rscript
# This is an R script that manages the re-running of scripts for B2B

options(stringsAsFactors=FALSE);
options(error=recover)

VERBOSE       <- TRUE
GLOBAL_ERRORS <- c("") # Problems are logged to this *GLOBAL* variable

file.nonzero.exists <- function(f) { return (file.exists(f)&&file.info(f)$size>0) }
is.nothing          <- function(x) { return (is.null(x) || is.na(x) || toupper(x) %in% c("NULL","N/A","NA","NONE","")) }
print0              <- function(...) { print(paste0(...)) }
statusPrint         <- function(...) { print0("[Status Update] ", ...); }
verboseStatusPrint  <- function(...) { if (VERBOSE) { statusPrint(...); } }
system0             <- function(...) { print0("[SYSTEM CALL]: ", ...); if (GLOBAL_DRY_RUN) { print0("(Not run--this is a DRY RUN") } else { system(paste0(...))} }
errlog              <- function(..., fatal=FALSE) { msg<-paste0(...); print0(msg); warning(msg); GLOBAL_ERRORS <<- append(GLOBAL_ERRORS, msg); if (fatal) { stop(msg); } }
warnlog             <- function(...) { msg<-paste0(...); print0(msg); warning(msg); GLOBAL_ERRORS <<- append(GLOBAL_ERRORS, msg); }
die_if_file_missing <- function(f, ...) {
     if(!file.exists(f)) {
          errlog("[Fatal error: missing file}: ", f, " ", ..., fatal=T)
          return(FALSE);
     }
     return(TRUE);
}
     

get_exe_path <- function(exe_name, how_to_install=NULL) {
     if (is.null(exe_name) || is.na(exe_name)) { return(NA); }
     result = tryCatch({
          exe_path <- system2("which", args=c(exe_name), stdout=T);
     }, error = function(e) {
          exe_path <- NULL # failure
     })
     if (length(exe_path) == 0 || is.na(exe_path) || is.null(exe_path) || !file.exists(exe_path)) {
          errlog("Failure to find the required executable named: ", exe_name)
          if (!is.null(how_to_install)) { errlog("You may be able to install it as follows:\n", how_to_install) }
     }
     return(exe_path)
}

write_global_errors_to_file <- function(filename) {
     fileConn <- file(ERR_LOG_FILE, open="w") # write to new file
     writeLines(GLOBAL_ERRORS, fileConn)
     close(fileConn)
}


record_exid_failure <- function(exid_that_failed, comment="") {
     FAILED_IDS[[exid]] <<- comment # GLOBAL VARIABLE
     errlog("[", exid, " ERROR]: ", comment)
}


get_java_and_jar_path <- function(jar_full_path, gigabytes_ram) {
     stopifnot(!missing(gigabytes_ram))
     return( paste0(get_exe_path("java"), " -Xmx", gigabytes_ram, "G -jar ", jar_full_path) )
}

# How to test:
#      ./0c_Analyze_kp_600_august_2016.R  

# Things you will have to set:
#        Set ASSAY_TYPE to either "RNASEQ" or "CHIPSEQ"

# This script does the following:
#      * Generates a file (ALIGN_CMD_FILE) with the qsub commands for running bowtie/tophat alignments on all relevant input files
#      * Generates a file (COMMANDS_FILE) with the downstream a qsub commands for running bowtie/tophat alignments on all relevant input files
#              * Also adds to this file the HTSEQ commands for counting reads.
#              * For CHIP-seq, also adds the 'gem' commands

# ************ CURRENTLY: also has the option to run edgeR analysis on RNAseq data, if "SHOULD_RUN_EDGER" is specified.
#         * note that this is run in this script itself, and NOT via qsub! That's not ideal

# Sample input file:
# X.out.RNA.4a_with_manual_fixes.txt

#INFILE <- "ZZZZ_TEST.txt" ; warning("USING A FAKE INPUT FILE")

software.list <- list("Tophat (splice-aware aligner)"                              =list(exe="tophat",      version="2.1.1" ,install="See details at: http://ccb.jhu.edu/software/tophat/index.shtml")
                    , "Bowtie (non-spliced aligner)"                               =list(exe="bowtie2",     version="2.2.4" ,install="See details at: http://bowtie-bio.sourceforge.net/index.shtml")
                    , "BCP (ChIP-seq peak caller for everything except TF binding)"=list(exe="BCP_HM",      version="1.1"   ,install="Available at https://cb.utdallas.edu/BCP/", install_comment="Paper: http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1002613")  # There are two versions of BCP - BCP_TF and BCP_HM. The latter is for histone marks - what you are interested in. BCP takes bed files as input.
                    , "GEM (motif-aware ChIP-seq peak caller for TF binding sites)"=list(exe="gem.jar",     version="2.5"   ,install="See instructions at: http://groups.csail.mit.edu/cgs/gem/", install_comment="We were using version 2.5, but version 2.7+ is available now.")
                    , "htseq-count (reads -> gene-level counts)"                   =list(exe="htseq-count", version="0.6.0" ,install="pip2.7 install HTseq", install_comment="Note: must be installed via pip (or other package manager). Do not just copy the binaries--it will not work.")
                    , "samtools"                                                   =list(exe="samtools",    version="1.3"   ,install="yum install samtools (on Redhat/CentOS)", install_comment="Available through your package manager. Other examples: brew install samtools (Mac Homebrew), apt-get install samtools (Ubuntu)")
                    , "edgeR (R library for differential expression)"              =list(exe=NA,            version="3.14.0",install="source('https://bioconductor.org/biocLite.R'); biocLite('edgeR')", install_comment="Available through Bioconductor.")
                    , "java"                                                       =list(exe="java"       , version="1.8.0" ,install="Install via package manager or on Oracle's web site: https://java.com/en/", install_comment="May also be possible to install using your package manager (e.g. 'yum install java'). Java is required to run 'gem.jar' and 'MarkDuplicates.jar'")
                    #, "test"=list(exe="testfdsafsda")
                      )

print0("------------------------------")
print0("Required software:")
for (desc in names(software.list)) {
     details <- software.list[[desc]]
     print0(desc)
     SPACER <- "     "
     print0(SPACER, "Executable name: ", details$exe, " (version ", details$version, ")")
     print0(SPACER, SPACER, "To install: ", details$install)
     if (!is.null(details$install_comment)) { print0(SPACER, SPACER, "            (", details$install_comment, ")") }
}
print0("------------------------------")

for (desc in names(software.list)) {
     details <- software.list[[desc]]
     get_exe_path(details$exe, paste0("How to install: ", details$install))
}

require("edgeR")


# ============= Mark duplicates =============
#MARK_DUP_JAR = "/data/applications/picard/picard-tools-1.114/MarkDuplicates.jar";
#dedup_cmd = paste(JAVA_EXE, " -Xmx15g -jar ", MARK_DUP_JAR, " INPUT=", finalbam, " REMOVE_DUPLICATES
# Note: although the above command RUNS the mark-duplicates command, it doesn't REMOVE duplicates---it just marks them

ASSAY_TYPE <- "RNASEQ"
#ASSAY_TYPE <- "CHIPSEQ"

BASEDIR <- getwd()

COMMANDS_FILE <- file.path(BASEDIR, "D1b_followup_commands_autogenerated.txt")
ALIGN_CMD_FILE = "D1a_alignment_commands_autogenerated.txt"

SHOULD_USE_FAKE_SIMULATED_DATA <- FALSE   # Should we use fake simulated data? Set to TRUE to not use real files
SIM_STATUS                     <- ifelse(SHOULD_USE_FAKE_SIMULATED_DATA, "--FAKE_SIMULATED_DATA--", "")

ERR_LOG_FILE  <- file.path(BASEDIR, "Z.errors.txt") # <-- print the GLOBAL_ERRORS variable to this file afterward

FAILED_IDS    <- list()
OK_IDS        <- list()


GEM_RAM_GB <- 25  # in gigabytes. Sometimes crsahes if it's < 10

GLOBAL_DRY_RUN     <- TRUE
SAMTOOLS_EXE       <- get_exe_path("samtools")
TOPHAT_N_THREADS   <- 4
BOWTIE_N_THREADS   <- 4
SAMTOOLS_N_THREADS <- 4
COMMENT_OUT_SHELL_COMMAND_STR <- "##### "

GLOBAL_QSUB_COUNT <- 1 # 
GLOBAL_QSUB_PREFIX <- format(Sys.time(), "%d_%H:%M:%S") # should be something unique to each run of this program

B2FRAME            <- "This_is_a_B2B_Data_Frame_Type" # just some text to use for sanity checking inputs later
GEM_READ_DIST_FILE <- file.path(BASEDIR, "DataFiles/Read_Distribution_default.txt")
COUNT_DIR          <- file.path(BASEDIR, "D5_htseq_count")
OUT_EDGER_DIR      <- file.path(BASEDIR, "E1_rnaseq_edger_diff_expr")
OUT_GEM_DIR        <- file.path(BASEDIR, "F1a_chipseq_gem")
OUT_BCP_DIR        <- file.path(BASEDIR, "F1b_chipseq_bcp")

NA_STRINGS         <- c("NA","na","NULL","null","N/A","n/a")

clear_out_file <- function(filename) { connection <- file(COMMANDS_FILE, open="w"); close(connection); stopifnot(file.exists(filename)); }
append_commands <- function(filename, lines, commented_out=FALSE) {
     connection <- file(COMMANDS_FILE, open="a");
     if (commented_out) { lines <- paste("### ", lines) } # should we actually print these as COMMENTED OUT lines?
     writeLines(lines, connection);
     close(connection);
}

read_sample_info_from_file <- function(filename) {
     REQUIRED_COLUMN_NAMES <- c("species", "file1", "file2", "proj_type", "expr_group", "chip_input", "pubtype") # These must all be the headers in 'filename'. Additional columns are ALLOWED.
     statusPrint("Reading the 'all' data matrix from input file ", INFILE)
     dall     <- read.table(filename, sep="\t", header=T, row.names=1, na.strings=NA_STRINGS)
     exid.vec <- gsub("X.*", "R", rownames(dall), perl=T, ignore.case=T) # full experiment ID. So for example, an sample ID might be "399X2". But the actual EXPERIMENT id there is "399R." In other words: take off the "X" and everything beyond, and add R
     missing_headers.vec <- REQUIRED_COLUMN_NAMES[!(REQUIRED_COLUMN_NAMES %in% colnames(dall))]
     if (length(missing_headers.vec) > 0) {
          errlog("Fatal error--in the input file '", filename, "', we need to find a set of required column names. However, the following names were missing: ", paste(missing_headers.vec, collapse=", "))
          stop("Improperly formatted input file! Check the header names. Note: the case must also match! Usually this means all lower-case.")
     }
     final_frame <- data.frame(dall, "exid"=exid.vec)
     attr(final_frame, B2FRAME) <- B2FRAME # <-- dubious way to indicate that this is our SPECIFIC kind of data frame
     return(final_frame)
}

assert_type_is_b2frame <- function(x) {
     if (B2FRAME != attr(x, B2FRAME)) {
          errlog("Warning! It looks like an incorrect data frame type was passed into this function. Exiting now!...");
          stop(paste0("The input variable must have an 'attr' that is set to the value '", B2FRAME, "'"))
     }
}

get_unique_experiment_ids_from_b2frame <- function(df) { # returns a vector of unique experiment ids
     assert_type_is_b2frame(df)
     exid.vec <- unique(df[,"exid",drop=TRUE])
     return(exid.vec)
}

get_species_from_b2frame <- function(df) { # returns a single string with the species in it. ERROR raised if there are multiple species.
     assert_type_is_b2frame(df)
     species <- unique(df[["species"]])
     if (1!=length(species)) {
          errlog("Somehow you appear to have more than one (different) species indicated for this experiment!")
          stop("Multiple species in one experiment is not currently supported!")
     }
     return(species)
}

get_b2frame_for_exid_from_b2frame <- function(df, experiment_id) {
     assert_type_is_b2frame(df)
     finalframe <- df[ df$exid == experiment_id, , drop=F] # only the records for THIS SPECIFIC experiment
     return(finalframe)
}

get_sids_for_exid_from_b2frame <- function(df, experiment_id) { # gets a vector of sample IDs associated with this experiment ID
     assert_type_is_b2frame(df)
     interesting_rownames.vec <- sort(unique(rownames(df[ df$exid == experiment_id, , drop=F]))) # only the records for THIS SPECIFIC experiment
     stopifnot(length(interesting_rownames.vec) > 0)
     return(interesting_rownames.vec) # rownames are the sample IDs
}

clusterize_command <- function(cmd) { # Returns a qsub wrapper for the command 'cmd'
     # Here's a weird 'feature' of R: we are keeping track of the NUMBER of qsub jobs here by adding an attribute to this function.
     # (This is basically a static variable).
     # This is TERRIBLE practice.
     QSUB_EXE <- system2("which", args=c("qsub"),  stdout=T); stopifnot(file.exists(QSUB_EXE))
     GLOBAL_QSUB_COUNT <<- (GLOBAL_QSUB_COUNT+1) # increment this GLOBAL variable each time a job is submitted
     username <- Sys.info()[["user"]]
     jobname  <- paste0("B2B_", GLOBAL_QSUB_PREFIX, "_", sprintf("%05d", GLOBAL_QSUB_COUNT), "_", username, "_", randint)
     cmd      <- paste0("echo \"", cmd, "\" | ", QSUB_EXE, " -V -N ", "'", jobname, "'")
     return(cmd)
}

agw_construct_cmd <- function(...) {
     return(paste0(...)) # may have to add qsub stuff here
}

agw_get_annotation <- function(assembly, file_type, must_exist=TRUE) {
     stopifnot(file_type %in% c("bowtie_index", "gtf", "fasta", "fa_by_chrom", "chr_sizes")) # sanity check the input 'file_type' string variable
     mapping <- list("hg19"     ="/data/info/genome/hg19_ensembl_igenome_with_chr/hg19.chr"
                     , "mm9"    ="/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr"
                     , "danRer7"="/data/info/genome/danRer7_ensembl_igenome_ucsc/danRer7_ucsc"
                     , "galGal4"="/data/info/genome/galGal4_ensembl_igenome_ucsc/galGal4_ucsc"
                     )
     pre = mapping[[assembly]] # file prefix
     if (is.null(pre)) { stop(paste0("error, <", assembly, "> is an unrecognized assembly -- note that you must use the UCSC id for the assembly (like 'hg19' or 'mm9'), not the common species name (so not 'mouse')")) }
     
     file_list <- list("bowtie_index"=paste0(pre), "gtf"=paste0(pre,".gtf"), "fasta"=paste0(pre,".fa")
                     , "fa_by_chrom"=paste0(pre,"_chromosome_fastas"), "chr_sizes"=paste0(pre,"chrom.sizes.txt"))
     the_file <- file_list[[file_type]]
     if (must_exist && !file.nonzero.exists(the_file) && !file.nonzero.exists(paste0(the_file, ".1.bt2"))) {
          errlog("Missing the REQUIRED input file of type '", file_type, "', which we expected (but failed) to find at the following location: ", the_file)
          stop("Missing a REQUIRED input file.")
     }
     return(the_file)
}

agw_align <- function(exid, fq1, fq2=NULL, outdir, aligner="none") {
     stopifnot(aligner %in% c("bowtie","tophat"))	
     die_if_file_missing(fq1, "[ERR]: Input FASTQ file ", fq1, " (forward (pair 1) and/or unpaired) seems to be missing!")
     (is.nothing(fq2) || die_if_file_missing(fq2, "[ERR]: Input FASTQ file ", fq2, " (#2, reverse pair file) seems to be missing!"))

     isPaired <- !is.nothing(fq2)
     innerDistStr <- ifelse(isPaired, " --mate-inner-dist=150 ", " ")

     btie   <- agw_get_annotation(species, "bowtie_index")
     ###fa   <- agw_get_annotation(species, "fasta")

     all_reads_bam         <- file.path(outdir, "all.bam")
     filtered_final_bam    <- file.path(outdir, "accepted_hits.bam")

     cmd <- paste0("mkdir -p ", outdir)
     if (aligner=="tophat") {
          gtf      <- agw_get_annotation(species, "gtf")    # <-- only actually used by tophat
          pair2Str <- ifelse(isPaired, yes=paste0(" ", fq2), no="") # might be blank, if unpaired input
          cmd <- paste0(cmd, "\n"
                      , get_exe_path('tophat'), " -o ", outdir
                      , " --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search "
                      , " --GTF=", gtf, " --num-threads=", TOPHAT_N_THREADS, " ", innerDistStr, " ", btie, " ", fq1, " ", pair2Str)
          
     } else if (aligner == "bowtie") {
          inputFqString <- ifelse(isPaired, yes=paste0(" -1 ", fq1," -2 ", fq2), no=paste0(" -U ", fq1))
          cmd <- paste0(cmd, "\n"
                      , get_exe_path('bowtie2'), " --threads=", BOWTIE_N_THREADS, " -x ", btie, " ", inputFqString
                      , " | " , SAMTOOLS_EXE, " view -@ ", SAMTOOLS_N_THREADS, " -b -                            > ", all_reads_bam
                      , " && ", SAMTOOLS_EXE, " view -@ ", SAMTOOLS_N_THREADS, " -b -F 0x104 ", everythingbam, " > ", filtered_final_bam)
          warning("are these position sorted? I think they are not")
          # samtools view -b -F 0x104 in.bam > mapped_primary.bam This is how you get ONLY primary mapped reads
     } else {
          stop("unknown aligner: only bowtie and tophat are allowed for now!")
     }
     return(agw_construct_cmd(cmd));
}

agw_get_chip_input_id <- function(exid, sid, dframe) {
     thisdat <- dframe[ dframe[,"exid"]==exid, , drop=F]  # only the records for THIS SPECIFIC experiment
     sampdat <- thisdat[rownames(thisdat)==sid, , drop=F][1,,drop=F]
     inid <- sampdat$chip_input; stopifnot(length(inid) == 1)
     return(inid)
}

agw_get_chip_input_bam_for_sample <- function(exid, sid, dframe) {
     inid <- agw_get_chip_input_id(exid, sid, dframe)
     return (ifelse(is.nothing(inid), yes=NA, no=agw_get_bam_path(experiment_id=exid, sample_id=inid)))
}

agw_sample_is_an_input <- function(exid, sid, df) {
     # return whether the sample is a chipseq 'input' file (i.e., a control).
     # Must start with the literal text "INPUT"
     return (grepl("^input", agw_get_chip_input_id(exid, sid, df), perl=T, ignore.case=T))
}

agw_get_gem_chipseq_cmd <- function(species, inputBam=NA, expBam, outPrefix) {
     # inputbam means; the INPUT file for chipseq. There MIGHT possibly not be one.
     # expbam is the (required) experimental bam file. 
     threads <- 1 # gem uses a ton of threads for some reason, so 1 is even conservative
     qval    <- 2	      # sets the q-value threshold (-log10, so 2 = 0.01)
     kmin    <- 6	      # minimum kmer length. From Monkey.
     kmax    <- 13	      # maximum kmer length. From Monkey.
     genomeByChromFastaDir <- agw_get_annotation(species, "fa_by_chrom")
     chrSizesFile          <- agw_get_annotation(species, "chr_sizes")

     if (!is.nothing(inputBam) && !file.nonzero.exists(inputBam)) {
          errlog("The missing BAM file we were looking for was named <<", inputBam, ">>. It was part of experiment with output prefix <", outPrefix, ">.")
          return(paste0("echo 'failure for", inputBam, "'"))
          #stopifnot(is.na(inputBam) || file.exists(inputBam))
     }
     stopifnot(file.exists(expBam))

     #GEM's help page is at http://groups.csail.mit.edu/cgs/gem/
     #Below is an example command I used to run GEM on Tbx5, a TF, Chip-seq
     #         java -Xmx10G -jar /home/rthomas/GEM/gem.jar --q 2 --d /home/rthomas/GEM/Read_Distribution_default.txt /home/rthomas/Annotation/mouse/mm9/mm9.chrom.sizes --genome /data/info/genome/mm9_ensembl_igenome_with_chr/mm9_per_chr_fasta_dir   --expt /home/rthomas/ChipSeqChallenge/RealData/04a.mapping/Tbx5_CM_WT_A.1_mm9.chr_q30.bam  --ctrl /home/rthomas/ChipSeqChallenge/RealData/04a.mapping/Input_CM.1_mm9.chr_q30.bam  --f BAM  --out Tbx5_Chipseq_WT_CorrectMM9 --k_min 6 --k_max 13 --t 8
     gemCmd <- paste0(get_java_and_jar_path(jar_full_path=get_exe_path("gem.jar"), gigabytes_ram=GEM_RAM_GB)
                    , " --t ", threads
                    , " --q ", qval
                    , " --k_min ", kmin
                    , " --k_max ", kmax
                    , " --d ", GEM_READ_DIST_FILE  # readdist file / read distribution file
                    , " --g ", chrSizesFile
                    , " --genome ", genomeByChromFastaDir
                    , " --expt ", expBam
                    , ifelse(is.na(inputBam), " ", paste0(" --ctrl ", inputBam))  # can be omitted if there is NO input file
                    , " --f SAM --sl --outBED"
                    , " --out ", outPrefix) # outprefix creates a new folder with this name, respecting existing subdirectories
     return(agw_construct_cmd(gemCmd))
}

agw_get_bcp_chipseq_cmd <- function() {

     # Attached is the manual for running BCP. There are two versions of BCP - BCP_TF and BCP_HM. The latter is for histone marks - what you are interested in. BCP takes bed files as input.
     #This is an example command I ran on a sample of H3K4me3 data:
     #BCP_HM -1 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1.bed -2 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed -f 200 -w 200 -p 0.001 -3 GM12878_H3K4me3_Rep1_EncodeAlign_results001_HM.bed
     
     #Note, there are changes to this command if you have Chip-exo data. You will need to replace  "Read_Distribution_default.txt" by "Read_Distribution_ChIP-exo.txt" and will need to add --smooth 3 option to the commas.
     #-Reuben

     # convert input
     warning("we need to convert bam1 to bed1")
     warning("we need to convert bam2 to bed2")
     warning("we need to convert bam3 to bed3")
     
     #BCP_HM -1 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1.bed
     #       -2 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed
     #       -f 200 -w 200 -p 0.001 -3 GM12878_H3K4me3_Rep1_EncodeAlign_results001_HM.bed
     bed1 <- "/home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1.bed"
     bed2 <- "/home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed"
     bed3 <- "GM12878_H3K4me3_Rep1_EncodeAlign_results001_HM.bed"
     f    <- 200
     w    <- 200
     p    <-   0.001
     
     bcpCmd <- paste0(get_exe_path("BCP_HM")
                      , " -1 ", bed1
                      , " -2 ", bed2
                      , " -f ", f
                      , " -w ", w
                      , " -p ", p
                      , " -3 ", bed3)
     stop("Not implemented")
     return(agw_construct_cmd(bcpCmd))
}

agw_get_bam_path <- function(experiment_id, sample_id) {
     bam_path <- file.path(BAM_DIR, experiment_id, sample_id, paste0("accepted_hits.bam"))
     return(bam_path)
}

agw_get_count_path <- function(experiment_id, sample_id) {
     count_path <- file.path(COUNT_DIR, experiment_id, sample_id, paste0("htseq_count.",experiment_id,"-",sample_id, ".counts.txt"))
     return(list("path"=count_path, "dir"=file.path(COUNT_DIR, experiment_id, sample_id)))
}

agw_get_gem_path <- function(experiment_id, sample_id) {
     parentDir <- file.path(OUT_GEM_DIR, experiment_id, sample_id)
     prefixBase <- "gem"
     return(list("prefix"=file.path(parentDir, prefixBase)
               , "dir"=parentDir
               , "created_file"=file.path(parentDir, prefixBase, paste0(prefixBase,"_result.htm"))))
}

agw_de_two_groups <- function(theRG, exprDatTab, g1, g2, theGroupVec, writeOutputToThisFile=FALSE) {
     # g1 and g2 are both NUMERIC
     stopifnot(g1 != g2)
     stopifnot(g1 >= 1);
     stopifnot(g2 >= 1);
     D <- DGEList(counts=theRG, group=theGroupVec)
     D <- edgeR::estimateCommonDisp(D)       
     D <- edgeR::estimateTagwiseDisp(D) # exactTest for negative bionomial distribution
     groupCounts        <- table(exprDatTab$expr_group)
     bothHaveReplicates <- (groupCounts[[g1]] >= 2 && groupCounts[[g2]] >= 2)
     if (bothHaveReplicates) {
          et <- edgeR::exactTest(D, pair=c(g1, g2))
     } else {
          HARD_CODED_DISPERSION <- 0.1 # may be low for human (0.4 seems recommended?). 0.1 is recommended for mouse and other genetically-identical samples though
          # Note about HARD_CODED_DISPERSION: this number was used by Stacia during initial analysis, so we need to keep it the same here too!
          et <- edgeR::exactTest(D, pair=c(g1, g2), dispersion=HARD_CODED_DISPERSION)
     }
     tt <- data.frame(topTags(et, n=nrow(et), sort.by="none"), check.names=FALSE) # gets the FDR calculated, too
     t.df <- tt[, c("logFC", "PValue", "FDR")] # only the three columns of interest---omit the logCPM column
     colnames(t.df) <- paste( paste0(g1,":",g2), colnames(t.df), sep=" ")

     if (writeOutputToThisFile) {
          outTable <- data.frame("IDENTIFIER"=rownames(t.df), t.df, check.names=FALSE)
          print0("Writing to the output file <", writeOutputToThisFile, ">")
          write.table(outTable, file=writeOutputToThisFile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
     }
     return(t.df) # Return the table with logFC, PValue, FDR. It's a data frame now, instead of whatever topTags natively returns.
}

agw_handle_rna_diff_expression <- function(exid, sid, df, outdir, cmdfile) {
     assert_type_is_b2frame(df)
     groupcounts             <- table(df$expr_group)
     smallest_num_replicates <- min(groupcounts)
     stopifnot(!missing(outdir)); stopifnot(!grep(" ", outdir)) # no spaces in outdir!
     append_commands(cmdfile, agw_construct_cmd(paste0("mkdir -p ", outdir)))
     edgerOut <- file.path(outdir, paste0("edgeR_differential_expression_", exid, "", SIM_STATUS, ".tsv.txt"))
     if (smallest_num_replicates < 2) {
          record_exid_failure(exid, "No replicates, skip this comparison FOR NOW ONLY.")
     } else if (!all(sapply(count_filenames.vec, file.nonzero.exists))) {
          record_exid_failure(exid, "EDGER IS MISSING SOME HTSEQ INPUT FILES for that experiment ID so we are skipping this differential expression.")
     } else if (file.exists(edgerOut)) {
          errlog(exid, " Looks like the edgeR out ALREADY EXISTS for experiment ID ", exid, ", in the location <", edgerOut, "> so we will not be re-running it.")
          OK_IDS[[exid]] <- "OK" # we (probably) did this properly
     } else if (length(count_filenames.vec) != length(df$expr_group)) {
          record_exid_failure(exid, paste0("missing input files for experiment ", exid, ": The number of 'count filenames' is NOT equal to the length of the group vector. This means some experimental groups were not aligned and/or processed downstream."))
          FAILED_IDS[[exid]] <- 1
     } else {
          print0("Running EdgeR for experiment ", exid, "...")
          stopifnot(all(sapply(count_filenames.vec, file.nonzero.exists))) # all the input files must ALREADY exist, please!
          RG <- edgeR::readDGE(count_filenames.vec, header=FALSE) # htseq-count does NOT have a header!
          # """Each file is assumed to contain digital gene expression data for
          # one genomic sample or count library, with gene identifiers in the
          # first column and counts in the second column. Gene identifiers are
          # assumed to be unique and not repeated in any one file.  The
          # A count of zero will be entered for any
          # gene that was not found in any particular sample.
          #By default, the files are assumed to be tab-delimited and to
          # contain column headings. Other file formats can be handled by..."""
          
          # If there are replicates:
          # Figure out groups for pairwise analysis
          groups <- df$expr_group
          stopifnot(is.integer(groups))
          if (any(sort(unique(groups)) != seq_along(unique(groups)))) {
               errlog("ERROR 162F: the groups were NOT numbered 1-through-n with no gaps! This indicates something we did not expect in the input data.", fatal=TRUE)
               stop()
          }
          
          elist = list()
          for (g1 in sort(unique(groups))) { # g1 and g2 are both NUMERIC GROUP IDs
               for (g2 in sort(unique(groups))) {
                    if (g1 >= g2) { next; } # don't double-calculate, or calculate g1 vs itself
                    print0("Comparing group ", g1, " vs ", g2, " (comparison name is ", g1, ":", g2, ")...")
                    elist[[ 1+length(elist) ]] <- agw_de_two_groups(theRG=RG, exprDatTab=df, g1=g1, g2=g2, theGroupVec=groups)
               }
          }
          if (length(elist) > 0) {
               eframe <- do.call(cbind, elist)
               outTable <- data.frame("IDENTIFIER"=rownames(elist[[1]]), eframe, check.names=FALSE)
               write.table(outTable, file=edgerOut, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
               OK_IDS[[exid]] <- "OK" # we (probably) did this properly
          } else {
               edgerOut <- paste0(edgerOut, "--FAILED-TO-GENERATE-RESULTS.txt")
               write.table(c("no results generated"), file=edgerOut)
               record_exid_failure(exid, "Warning: did not find any group comparisons for this experiment! This may be a programming error (?) or it may be due to unusual input data.")
          }
          print0("Wrote results to the output file <", edgerOut, ">")
          #topTags(et)
          # If there are no replicates:
          # D <- DGEList(counts=RG,group=1:2:3:4)
          # et <- exactTest(D,dispersion=0.1)
          ## topTags(et)
     }
}

agw_handle_rna_counting_per_feature <- function(exid, sid, species, cmdfile) {
     bam <- agw_get_bam_path(experiment_id=exid, sample_id=sid)
     gtf <- agw_get_annotation(species, "gtf"); stopifnot(file.nonzero.exists(gtf))
     ccc    <- agw_get_count_path(experiment_id=exid, sample_id=sid)$"path"
     cccDir <- agw_get_count_path(experiment_id=exid, sample_id=sid)$"dir"
     print0("Validating that counts file <", ccc, "> exists for sample ", sid)
     #htCmd <- paste0("mkdir -p ", cccDir, "; ", HTSEQ_COUNT_EXE, " --format=bam --order=pos --mode=intersection-nonempty --stranded=no --minaqual=10 --type=exon --idattr=gene_id ", bam, " ", gtf, " >| ", ccc)
     bamNameSorted <- paste0(bam, ".sorted_by_name.bam")
     htCmd <- paste0("mkdir -p ", cccDir, "; ", SAMTOOLS_EXE, " sort -n ", bam, " > ", bamNameSorted, " && ", get_exe_path('htseq-count'), " --format=bam --order=name --mode=intersection-nonempty --stranded=no --minaqual=10 --type=exon --idattr=gene_id ", bamNameSorted, " ", gtf, " >| ", ccc)
     if (!file.nonzero.exists(ccc)) {
          append_commands(cmdfile, agw_construct_cmd(htCmd)) # the command for running htseq-count
     } else {
          verboseStatusPrint("Not re-running HTSeq, because the output file already exists.")
          append_commands(cmdfile, agw_construct_cmd(htCmd), comment_out=TRUE)
     }
     count_filenames.vec <- append(count_filenames.vec, ccc)
}

agw_handle_chip_peak_calls <- function(exid, sid, df_for_this_exid, cmdfile) {
     if (agw_sample_is_an_input(exid, sid, df_for_this_exid)) { # WE don't run GEM or BCP on an input sample by itself!
          verboseStatusPrint("No need to peak-call an INPUT file (", exid, " / ", sid, ")")
          return;
     }

     peak_type <- "TF" # or "HISTONE" or "NON-TF"
     
     if (peak_type == "TF") {
          ggg    <- agw_get_gem_path(exid, sid)
          input  <- agw_get_chip_input_bam_for_sample(exid, sid, df_for_this_exid)
          if (is.na(input)) {
               msg <- paste0("[", exid, " FAILURE]: Failure to find a matching INPUT sample for experiment ", exid, " / sample ", sid, " -- not running GEM")
               errlog(msg)
               gemCmd <- paste0("# ", msg)
          } else {
               print0("[OK] Generating a list of commands to ", cmdfile)
               gemCmd <- paste0("mkdir -p ", ggg$dir, " && ", agw_get_gem_chipseq_cmd(species, inputBam=input, expBam=bam, outPrefix=ggg$prefix))
          }
          do_outputs_already_exist <- file.nonzero.exists(ggg$created_file) # if the output file ALREADY EXISTS, then comment out the command
          append_commands(cmdfile, agw_construct_cmd(gemCmd), commented_out=do_outputs_already_exist) # append the command for running gem
     } else if (peak_type == "HISTONE" || peak_type == "NON-TF") {
          stop("can't handle BCP so far")
     } else {
          stop("unknown peak type")
     }
}

# =====================  DONE WITH FUNCTIONS ==========================

#MARK_DUP_JAR = "/data/applications/picard/picard-tools-1.114/MarkDuplicates.jar";
if (ASSAY_TYPE == "CHIPSEQ") {
     INFILE      <- file.path(BASEDIR, "X.out.CHIP.3d_chipseq_first_entry_only_for_dupes.txt")
     BAM_DIR     <- file.path(BASEDIR, "D2b_bowtie_output")
     ALIGNED_DIR <- file.path(getwd(), "D2b_bowtie_output")
     ALIGNER     <- "bowtie2"
} else if (ASSAY_TYPE == "RNASEQ") {
     INFILE      <- file.path(BASEDIR, "X.out.RNA.4a_with_manual_fixes.txt")
     BAM_DIR     <- file.path(BASEDIR, "D2a_tophat_output")
     ALIGNED_DIR <- file.path(getwd(), "D2a_tophat_output")
     ALIGNER     <- "tophat"
} else {
     stop("unrecognized assay type")
}

clear_out_file(COMMANDS_FILE)
append_commands(COMMANDS_FILE, "echo 'Starting the set of commands for this specific input file'")

statusPrint("Reading the 'all' data matrix from input file ", INFILE)
alldat      <- read_sample_info_from_file(INFILE)

append_commands(COMMANDS_FILE, "echo '[DONE] with the commands for this specific input file'")

#for (exid in "24R") { # unique(dat[,"exid"])) {
for (exid in get_unique_experiment_ids_from_b2frame(alldat)) {
     subdat                      <- get_b2frame_for_exid_from_b2frame(alldat, exid)
     sample_ids_for_exid.vec     <- get_sids_for_exid_from_b2frame(df=subdat, experiment_id=exid)
     species                     <- get_species_from_b2frame(subdat)
     print0("Handling experiment ", exid, " (", species, "). The full set of data for this experiment is:"); print(subdat)
     f1.vec <- subdat[["file1"]]
     f2.vec <- subdat[["file2"]]

     # handle the alignments
     for (ridx in seq_along(rownames(subdat))) {
          item      <- rownames(subdat)[ridx] # item name (unique)
          exid      <- subdat[ridx,    "exid"]  # experiment ID (shared)
          file1     <- subdat[ridx,   "file1"]  # pair 1 for paired-end
          file2     <- subdat[ridx,   "file2"]  # pair 2 for paired-end
          aligner   <- ALIGNER
          paired1   <- file.path(BASEDIR, file1)
          paired2   <- if (is.nothing(file2)) { NULL } else { file.path(BASEDIR, file2) }
          outdir    <- file.path(ALIGNED_DIR, exid, item)
          align_cmd <- agw_align(exid=exid, fq1=paired1, fq2=paired2, outdir=outdir, aligner=aligner)
          append_commands(COMMANDS_FILE, align_cmd) # align it!
     }
     
     exid_is_missing_a_sample <- FALSE
     
     if (SHOULD_USE_FAKE_SIMULATED_DATA) {
          warning("Note: using FAKE SIMULATED data!")
          count_filenames.vec <- paste("FAKE.", seq_along(f1.vec), ".counts.tmp", sep='')
          for (f in count_filenames.vec) {
               # Making a fake simulated file now!
               n <- 300 # <-- number of fake genes
               v <- as.integer(floor(rexp(n, 0.01) + 1))
               df <- data.frame(v)
               rownames(df) <- paste("ENSFAKE_", seq(from=1, to=n), sep="")
               write.table(df, file=f, col.names=FALSE, quote=FALSE, row.names=T, sep="\t")
          }
     } else {
          # Get real data
          count_filenames.vec <- c()
          for (sid in sample_ids_for_exid.vec) {
               bam <- agw_get_bam_path(experiment_id=exid, sample_id=sid)
               if (!file.nonzero.exists(bam)) {
                    errlog("[", exid, " FAILURE]: [MISSING BAM FILE in experiment ID ", exid, ", sample ID ", sid, ". Specifically, file <", bam, "> does not exist.")
                    exid_is_missing_a_sample <- TRUE
               } else {
                    if      (ASSAY_TYPE == "RNASEQ" ) { agw_handle_rna_counting_per_feature(exid=exid, sid=sid, species=species, cmdfile=COMMANDS_FILE) }
                    else if (ASSAY_TYPE == "CHIPSEQ") { agw_handle_chip_peak_calls(exid=exid, sid=sid, df_for_this_exid=subdat, cmdfile=COMMANDS_FILE) }
                    else { stop("incorrect assay type") }
               }
          }

          if (exid_is_missing_a_sample) {
               record_exid_failure(exid, paste0("[SKIPPING all remaining analysis steps for experiment ID ", exid, " due to a missing file issue earlier. See the logs for the exact sample ID that was missing."))
          } else {
               if (ASSAY_TYPE == "RNASEQ") {
                    agw_handle_rna_diff_expression(exid, sid, subdat, outdir=file.path(OUT_EDGER_DIR, exid), cmdfile=COMMANDS_FILE)
               }
          }
          #browser()
     }
}

write_global_errors_to_file(filename=ERR_LOG_FILE)

print("---------------")
print(GLOBAL_ERRORS, collapse="\n")
print("---------------")
if (length(FAILED_IDS) > 0) {
     print0("ERROR: Failure for the following ", length(FAILED_IDS), " experiment IDs:")
     print(paste(names(FAILED_IDS), collapse=", "))
} else {
     print0("[OK] Commands were successfully written for every ID.")
}
print0("Appeared to succeed with the following ", length(OK_IDS), " experiment IDs:")
print(paste(names(OK_IDS), collapse=", "))
print("---------------")
