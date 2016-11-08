#!/usr/bin/env python2

# ./1_pipeline.py --rna-samples=Z3_Test_Data/a1.mm9.chr19.fq.gz --groups='1'

from __future__ import print_function # sets print("...") to have parens like in python 3
from __future__ import division # defaults to floating point division. No more "1/2 == 0"
import sys
import re
import pdb #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques
import os.path
import argparse # requires python 2.7

ASSAY_RNA  = 1
ASSAY_CHIP = 2


global_file_paths_we_already_aligned = dict()
global_num_script_files_written = 0

# ======= Paths to the programs we need. Assumed to be on the user's $PATH ========
TOPHAT_PATH   = "tophat"
BOWTIE_PATH  = "bowtie2"
SAMTOOLS_PATH = "samtools"

ALIGNER_TOPHAT  = TOPHAT_PATH
ALIGNER_BOWTIE = BOWTIE_PATH

TOPHAT_N_THREADS   = 4
BOWTIE_N_THREADS  = 4
SAMTOOLS_N_THREADS = BOWTIE_N_THREADS

NUM_ZEROS_TO_PAD_IN_SCRIPT_NAMES = 4
EXIT_CODE_IF_MISSING_FILE = 49 # just some arbitrary non-zero number
LITERAL_BACKSLASH="\\" # reduces to just one backslash
LITERAL_DQUOTE="\""

DEFAULT_SCRIPT_PREFIX     = "./PIPELINE_SCRIPT_"
DEFAULT_SPLICED_ALIGN_DIR = "./Pipeline_01_Spliced_Tophat_Alignment"
DEFAULT_DNA_ALIGN_DIR     = "./Pipeline_01_DNA_Bowtie_Alignment"

opt = None # <-- "opt" is a global variable that stores the cmd line arguments. This should really be the ONLY non-constant global!

def enquote(s):
    # should we check to see if s has quotes in it already?
    return("\"" + str(s) + "\"")

def verbosePrint(m):
    if (opt.verbose): printStderr(m)
    return

def progressPrint(m):
    printStderr(m) # print even if we are NOT in verbose mode
    return

def printStderr(*args, **keyargs):
    print(*args, file=sys.stderr, **keyargs)

def withPrependedBasedir(list_of_files, basedir):
    if basedir is not None:
        list_of_files = [os.path.join(basedir, fff) for fff in list_of_files] # python list comprehension---add the basedir to the beginning of each file path
        pass
    return list_of_files

def dieBadArgs(errMsg):
    sys.exit("ERROR: Problem with the command line arguments! Specifically: " + errMsg + ". Try using '--help' to see all the possible arguments.")
    return

def handleCmdLineArgs():
    DESC = '''
A program for running a standard ChIP- or RNA-seq pipeline.
Written for Python 2.7+ by Alex Williams, 2016. (Python 2.7 is needed for the 'argparse' library)

What this program does:

1) You specify an experiment and the input FASTQ files, and perhaps a few other files.

2) You run a command line invocation like:
     python2  pipeline.py --rna-samples=/data/a1.fq.gz,/data/a2.fq.gz --groups=1,2 --species=mm9  --out="SCRIPTFILE_OUTPUT_PREFIX"

3) The output is a shell script that you can then submit on your cluster.
'''

    EPILOG = '''
Version history: (none yet)
'''

    parser = argparse.ArgumentParser(description=DESC, epilog=EPILOG, add_help=True, formatter_class=argparse.RawTextHelpFormatter)
    #, formatter_class=argparse.RawTextHelpFormatter)

    #parser.add_argument('--sup', '-s', dest='sup', default=None, metavar="\"UMI_LEFT,LINKER_LEFT,UMI_RIGHT,LINKER_RIGHT\"", type=str, help="Default: if unspecified (or set as --sup=GUESS), we will just guess the UMI lengths by looking at the consensus sequence of the first 1000 reads (and the linker wil be assumed to be ZERO bases--it will be included as part of the 'payload' sequence). 'Sup' stands for 'Supplementary' lengths (note: numbers, NOT actual sequence) for the UMI_LEFT, LINKER_LEFT, UMI_RIGHT, and LINKER_RIGHT.\nExample: --sup=\"8,2,2,5\"\n   * for a sample with an 8+5 UMI and linkers of length 2.")

    parser.add_argument('-v',  '--verbose', dest='verbose', action='store_true', help="Verbose debugging")
    parser.add_argument('--groups', '-g', type=str, default=None, dest='groups', help="Specify experiment groups (numeric). Must be the same length as the number of files.")

    parser.add_argument('--basedir', '-b', type=str, default=None, dest='basedir', help="Optional. Specify a base directory to prepend to all paths.!")

    parser.add_argument('--chip-samples', '-c', type=str, default=None, dest='chip_samples', help="ChIP-seq only: Specify ChIP-seq sample files, comma-delimited. Must match the order of the corresponding --chip_inputs files!")
    parser.add_argument('--chip-inputs' , '-i', type=str, default=None, dest='chip_inputs', help="ChIP-seq only: Specify ChIP-seq input files, comma-delimited. Must match the order of the corresponding --chip_samples files!")
    parser.add_argument('--rna-samples', '--rna', '-r', type=str, default=None, dest='rna_samples', help="RNA-seq only: Specify RNA-seq input files, comma-delimited.")

    parser.add_argument('--species', '-s', type=str, default=None, dest='species', help="Required: The short name for the species.")
    
    parser.add_argument('--delim', type=str, default=',', dest='delim', help="Default file/group delimiter. Normally should be a comma, unless you're doing something very unusual.")

    parser.add_argument('--align-dir', type=str, default=None, dest='align_dir', help="Default location to put the ALIGNED reads after running tophat or bowtie. Ideally a full path.")

    parser.add_argument('--out', '-o', type=str, default=DEFAULT_SCRIPT_PREFIX, dest='script_prefix', help="The output script files to write. For example, --out=PIPELINE would generate files with names like: PIPELINE_001.sh, PIPELINE_002.sh, etc. Can be a path with a folder component as well.")

    global opt # <-- critical, since we modify (global variable) "got" below. Do not remove this!
    opt = parser.parse_args()

    assay = None
    if (opt.species is None):
        dieBadArgs("Species must be specified! It must be one species for ALL samples currently.")
        pass

    if (opt.rna_samples is not None and (opt.chip_samples is not None or opt.chip_inputs is not None)):
        dieBadArgs("You cannot specify BOTH RNA-seq and ALSO ChIP-seq at the same time! Pick either RNA-seq or ChIP-seq")
        pass

    if (opt.groups is None):
        dieBadArgs("You must specify the groups that each sample belongs to, or 'NA' or a blank value.")
        pass

    if (opt.rna_samples is not None):
        assay = ASSAY_RNA
        if (opt.align_dir is None):
            opt.align_dir = DEFAULT_SPLICED_ALIGN_DIR
            pass
        samples = withPrependedBasedir(opt.rna_samples.split(opt.delim), opt.basedir)
        # Apparently we're running an RNA-seq project
        pass
    elif (opt.chip_samples is not None):
        assay = ASSAY_CHIP
        if (opt.align_dir is None):
            opt.align_dir = DEFAULT_DNA_ALIGN_DIR
            pass
        samples = withPrependedBasedir(opt.chip_samples.split(opt.delim), opt.basedir)
        if (opt.chip_inputs is None):
            dieBadArgs("If you specify --chip-samples, you have to also specify --chip-inputs, but it appears that was NOT specified.")
            pass
        inputs  = withPrependedBasedir(opt.chip_inputs.split(opt.delim), opt.basedir)
        # Apparently we're running a ChIP-seq project
        if len(inputs) != len(samples):
            dieBadArgs("Somehow your --chip-inputs was a different length from the --chip-samples! Fix this.")
            pass
        pass
    else:
        dieBadArgs("You have to AT LEAST have an RNA-seq sample to process or a ChIP-seq sample to process. Specify at least one!")
        pass



    groups = opt.groups.split(opt.delim)
    if len(groups) != len(samples):
        dieBadArgs("Somehow your '--groups' was not equal in length to your specified number of sample files!")
        pass
    
    verbosePrint("Groups: " + str(groups))

    # Finally, actually generate the commands...

    if assay == ASSAY_RNA:
        handleRNA(groups=groups,  species=opt.species, output_bam_dir=opt.align_dir, samples=samples,                script_prefix=opt.script_prefix)
    elif assay == ASSAY_CHIP:
        handleCHIP(groups=groups, species=opt.species, output_bam_dir=opt.align_dir, samples=samples, inputs=inputs, script_prefix=opt.script_prefix)
    else:
        sys.exit("Something went horribly wrong!")
        pass

    
    return  # End of command-line-reading function


def listify_remove_none(aaa):
    if not isinstance(aaa, (list, tuple)):
        aaa = [aaa] # wrap it in a list
    return [x for x in aaa if x is not None] # remove NONE elements
    
class OurScript(object):
    """Script text for writing to a file"""
    def __init__(self):
        self.lines = ["#!/bin/bash -u"
                      , "set -e"
                      , "set -o pipefail"]

    def append(self, text):
        self.lines.append(text) # add exactly one line
        return

    def appendCheckForRequiredFiles(self, file_list):
        file_list = listify_remove_none(file_list) # make sure it's a [ ] list, and not just like ONE file someone passed in
        if (len(file_list) == 0):
            return # no need to do anything if we have zero items
        for f in file_list:
            self.append("if [[ ! -e " + enquote(f) + " ]]; then")
            self.append("      echo '[ERROR] Cannot find the following required file, which we expect to exist: <" + f + ">")
            self.append("      exit " + str(EXIT_CODE_IF_MISSING_FILE)) # exit with a somewhat arbitrary code number if there is a missing file
            self.append("fi")
            pass
        return

    def appendScriptExitIfAllFilesExist(self, file_list):
        file_list = listify_remove_none(file_list) # make sure it's a [ ] list, and not just like ONE file someone passed in
        if (len(file_list) == 0):
            return # no need to do anything if we have zero items

        multiIfString = "if " + " && ".join(  [ "[[ -e " + enquote(g) + " ]]" for g in file_list ] ) # <-- python list comprehension that generates results like "if [[ -e myFile ]] && [[ -e file2 ]] && [[ -e file3 ]]"
        self.append(multiIfString)
        self.append("then")
        self.append("    echo '[OK] All the specified output files already exist--so we are exiting the script now.' ")
        self.append("    exit 0")
        self.append("fi")
        return
        

    def writeToDisk(self, script_prefix=None):
        if (script_prefix is None):
            print("\n".join(self.lines)) # just print to stdout
        else:
            global global_num_script_files_written
            scriptname = script_prefix + "_" + str(global_num_script_files_written).zfill(NUM_ZEROS_TO_PAD_IN_SCRIPT_NAMES) + ".sh"
            global_num_script_files_written += 1
            with open(scriptname, 'w') as fff:
                fff.write("\n".join(self.lines) + "\n")
                pass
            pass
            progressPrint("Wrote script data to the following filename: " + scriptname)
        return


def newScriptText():
    return 

def handleRNA(groups, species, output_bam_dir, samples, script_prefix):
    
    for samp in samples:
        generateAlignmentCmd(fastq1=samp, species=species, outdir=output_bam_dir, aligner=ALIGNER_TOPHAT, script_prefix=script_prefix+"_align_rna")
        pass

    pass

def handleCHIP(groups, species, output_bam_dir, samples, inputs, script_prefix):

    for samp,inpt in zip(samples, inputs): # get matching samples/inputs
        generateAlignmentCmd(fastq1=samp, species=species, outdir=output_bam_dir, aligner=ALIGNER_BOWTIE, script_prefix=script_prefix+"_align_chip_sample")
        generateAlignmentCmd(fastq1=inpt, species=species, outdir=output_bam_dir, aligner=ALIGNER_BOWTIE, script_prefix=script_prefix+"_align_chip_input")
        # now run gem or whatever
        # finally, do something else
        pass
    pass


def agw_get_annotation(species, thing):
    return "___PLACEHOLDER_FOR_ANNOTATION_FILE___"


def generateAlignmentCmd(fastq1, fastq2=None, species=None, outdir=None, aligner=None, script_prefix=None):
    '''Note that this is will only generate ONE alignment command per file, even if it's called multiple times.'''
    
    fqkey = " | ".join([str(fastq1), str(fastq2), str(species), str(aligner)]) # just a unique key for this pair of reads, species, and aligner
    if fqkey in global_file_paths_we_already_aligned:
        verbosePrint("[Skipping...]: Already generated an alignment command for the file/pair " + str(fastq1) + " and " + str(fastq2) + " with species " + str(species) + " and the aligner " + str(aligner))
        return # exit early
    global_file_paths_we_already_aligned[fqkey] = True # remember that we saw this file...
    
    isPaired = fastq2 is not None
    innerDistStr = " --mate-inner-dist=150 " if isPaired else " "

    bowtieIndex = agw_get_annotation(species, "bowtie_index")
    all_reads_bam = os.path.join(outdir, "all.bam")
    aligned_sorted_bam  = os.path.join(outdir, "accepted_hits.bam")

    if aligner == ALIGNER_TOPHAT:
        gtf = agw_get_annotation(species, "gtf")    # <-- only actually used by tophat
        cmdStr  = TOPHAT_PATH + " -o " + outdir 
        cmdStr += " " + "--min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search "
        cmdStr += " " + "--GTF=" + gtf
        cmdStr += " " + "--num-threads=" + str(TOPHAT_N_THREADS)
        cmdStr += " " + innerDistStr
        cmdStr += " " + bowtieIndex
        cmdStr += " " + fastq1 + " "
        cmdStr += " " + fastq2  if isPaired  else "" # might or might not have a second file in the pair!
        
    elif aligner == ALIGNER_BOWTIE:
        if isPaired:
            inputFqString = " -1 " + enquote(fastq1) + " -2 " + enquote(fastq2) + " "
        else:
            inputFqString = " -U " + enquote(fastq1) + " "
            pass
        cmdStr  = BOWTIE_PATH
        cmdStr += " " + "--threads=" + str(BOWTIE_N_THREADS)
        cmdStr += " " + "-x " + bowtieIndex
        cmdStr += " " + inputFqString
        cmdStr += " | "  + SAMTOOLS_PATH + " view -@ " + str(SAMTOOLS_N_THREADS) + " -b - " +                        " > " + all_reads_bam
        cmdStr += " && " + SAMTOOLS_PATH + " view -@ " + str(SAMTOOLS_N_THREADS) + " -b -F 0x104 " + all_reads_bam + " "           # samtools view -b -F 0x104 in.bam > mapped_primary.bam This is how you get ONLY primary mapped reads
        cmdStr += " | samtools sort - " + " > " + aligned_sorted_bam
    else:
        raise Exception("Unrecognized aligner! We currently only know about tophat and bowtie2. Note that this is case-sensitive. This is a CODING bug and should not be related to a user-specified file. Your specified aligner was: " + aligner)
        pass
    
    xcmd = OurScript()
    xcmd.appendScriptExitIfAllFilesExist(file_list=[aligned_sorted_bam])     # exit EARLY if all the output files already exist
    xcmd.appendCheckForRequiredFiles(file_list=[fastq1, fastq2]) # require that these files exist (or are 'None')

    xcmd.append("mkdir -p " + enquote(outdir))
    xcmd.append(cmdStr) # the actual alignment command!
    xcmd.appendCheckForRequiredFiles(file_list=[aligned_sorted_bam]) # make sure the file got generated!
    xcmd.writeToDisk(script_prefix=script_prefix)

    # if ALL the generated files already exist, then do not run this script!

        

    return

# Must come at the VERY END!
if __name__ == "__main__":
    handleCmdLineArgs()
    pass




xyz = '''


#!/usr/bin/env Rscript
# This is an R script that manages the re-running of scripts for B2B

# Required R packages:
#       * optparse (option parsing) library. Available in CRAN
#         Install it by typing this in R:   install.packages('optparse')



# by Alex Williams

# How to use it:

# Example usage:
#   Rscript   this_b2b_script_name.R   --type=RNA   --info=sample_info.rna.txt  --out=w.samples

# Or:
#   Rscript   this_b2b_script_name.R   --type=CHIP  --info=sample_info.chip.txt  --out=chip.samples


cerr <- function(...) cat(paste0(..., "\n"), file=base::stderr()) # Print to STDERR
cout <- function(...) cat(paste0(..., "\n"), file=base::stdout()) # Print to STDOUT

suppressPackageStartupMessages(require("optparse"))

if (!require("optparse")) {
     stop("This scripts requires the 'optparse' library in R. Install it with: install.packages('optparse'). (Optparse is used to read command line arguments)")
}

require("edgeR")


options(stringsAsFactors=FALSE)
if (interactive()) { options(error=recover) } # <-- only when running interactively, or else 'stop' and 'stopifnot' fail to work.

option_list <- list(
     make_option(c("-v", "--verbose"), action="store_true", default=FALSE
                 ,help="Print extra-verbose output")
   , make_option(c("-t", "--type"), type="character", default=NULL
                 ,help="REQUIRED. Specify the type of input: must be either --type=RNA or --type=CHIP. Case-insensitive.")
   , make_option(c("-i", "--info"), type="character", default=NULL
                 ,help="REQUIRED. Specify the 'info' filename with the matrix of metadata about each project.")

   , make_option(c("-o", "--out"), type="character", default=NULL
                ,help="REQUIRED. Specify an output file PREFIX for the generated commands.")
)
opt <- parse_args(OptionParser(option_list=option_list))

opt$type <- toupper(opt$type) # upper-case the analysis assay type ("rna" --> "RNA")
if (is.null(opt$out))  { stop("[MISSING REQUIRED '--out=' ARGUMENT]: You must specify a prefix for your output files. For example, --out=my_samples . Thne, the per-sample files would be written as: my_samples.1.cmd.txt, my_samples.2.cmd.txt, etc...") }
if (is.null(opt$type)) { stop("[MISSING REQUIRED '--type=' ARGUMENT]: Please specify an analysis type with the '--type=...' command line option. Valid options are: --type=RNA and --type=CHIP") }
if (is.null(opt$info)) { stop("[MISSING REQUIRED '--info=' ARGUMENT]: The info filename (of the matrix of data about each project) must be specified with --info=FILENAME.") }

#     opt$info      <- file.path(BASEDIR, "X.out.CHIP.3d_chipseq_first_entry_only_for_dupes.txt")
#     opt$info      <- file.path(BASEDIR, "X.out.RNA.4a_with_manual_fixes.txt")

GLOBAL_ERR <- c() # Problems are logged to this *GLOBAL* variable

#system0             <- function(...) { print0("[SYSTEM CALL]: ", ...); if (GLOBAL_DRY_RUN) { print0("(Not run--this is a DRY RUN") } else { system(paste0(...))} }
file.nonzero.exists <- function(f) { return (file.exists(f)&&file.info(f)$size>0) }
is.nothing          <- function(x) { return (is.null(x) || is.na(x) || toupper(x) %in% c("NULL","N/A","NA","NONE","")) }
print0              <- function(...) { cerr(paste0(...)) }
statusPrint         <- function(...) { print0("[Status Update] ", ...); }
verboseStatusPrint  <- function(...) { if (opt$verbose) { statusPrint(...); } }
errlog              <- function(..., fatal=FALSE) { msg<-paste0(...); print0(msg); warning(msg); GLOBAL_ERR <<- append(GLOBAL_ERR, msg); if (fatal) { stop(msg); } }
warnlog             <- function(...) { msg<-paste0(...); print0(msg); warning(msg); GLOBAL_ERR <<- append(GLOBAL_ERR, msg); }
die_if_file_missing <- function(f, ...) {
     if(!file.exists(f)) {
          errlog("[Fatal error: missing file}: ", f, " ", ..., fatal=T)
     }
     return(invisible(TRUE));
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
     # We have stored all the problems encountered into a list of strings named "GLOBAL_ERR" (global variable)
     # Now we will write all of these to a file.
     fileConn <- file(filename, open="w") # write to new file
     writeLines(GLOBAL_ERR, fileConn)
     close(fileConn)
     return(invisible(NULL))
}

remember_exid_failure <- function(exid_that_failed, comment="") {
     FAILED_IDS[[exid]] <<- comment # GLOBAL VARIABLE
     errlog("[", exid, " ERROR]: ", comment)
     return(invisible(NULL))
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
#      * Generates a file (opt$out) with the downstream a qsub commands for running bowtie/tophat alignments on all relevant input files
#              * Also adds to this file the HTSEQ commands for counting reads.
#              * For CHIP-seq, also adds the 'gem' commands

# ************ CURRENTLY: also has the option to run edgeR analysis on RNAseq data, if "SHOULD_RUN_EDGER" is specified.
#         * note that this is run in this script itself, and NOT via qsub! That's not ideal

# Sample input file:
# X.out.RNA.4a_with_manual_fixes.txt

software.list <- list("Tophat (splice-aware aligner)"                              =list(exe="tophat"     , version="2.1.1"        ,install="See details at: http://ccb.jhu.edu/software/tophat/index.shtml")
                    , "Bowtie (non-spliced aligner)"                               =list(exe="bowtie2"    , version="2.2.4"        ,install="See details at: http://bowtie-bio.sourceforge.net/index.shtml")
                    , "bwa (non-splice/transcript-aware alinger)"                  =list(exe="bwa"        , version="0.7.12-r1039" ,install="See details at: http://bio-bwa.sourceforge.net/")
                    , "java"                                                       =list(exe="java"       , version="1.8.0"        ,install="Install via package manager or on Oracle's web site: https://java.com/en/", install_comment="May also be possible to install using your package manager (e.g. 'yum install java'). Java is required to run 'gem.jar' and 'MarkDuplicates.jar'")
                    , "BCP (ChIP-seq peak caller for everything except TF binding)"=list(exe="BCP_HM"     , version="1.1"          ,install="Available at https://cb.utdallas.edu/BCP/", install_comment="Paper: http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1002613")  # There are two versions of BCP - BCP_TF and BCP_HM. The latter is for histone marks - what you are interested in. BCP takes bed files as input.
                    , "bam2bed (part of the 'bedops' suite)"                       =list(exe="bam2bed"    , version="2.4.20"       ,install="Available here: http://bedops.readthedocs.io/en/latest/content/installation.html", install_comment="This program is required only for running BCP in the CHIPseq pipeline--nothing else uses it. Paper: http://bioinformatics.oxfordjournals.org/content/28/14/1919.abstract")
                    , "GEM (motif-aware ChIP-seq peak caller for TF binding sites)"=list(exe="gem.jar"    , version="2.5"          ,install="See instructions at: http://groups.csail.mit.edu/cgs/gem/", install_comment="We were using version 2.5, but version 2.7+ is available now.")
                    , "htseq-count (reads -> gene-level counts)"                   =list(exe="htseq-count", version="0.6.0"        ,install="pip2.7 install HTseq", install_comment="Note: must be installed via pip (or other package manager). Do not just copy the binaries--it will not work.")
                    , "samtools"                                                   =list(exe="samtools"   , version="1.3"          ,install="yum install samtools (on Redhat/CentOS)", install_comment="Available through your package manager. Other examples: brew install samtools (Mac Homebrew), apt-get install samtools (Ubuntu)")
                    , "edgeR (R library for differential expression)"              =list(exe=NA           , version="3.14.0"       ,install="source('https://bioconductor.org/biocLite.R'); biocLite('edgeR')", install_comment="Available through Bioconductor.")
                    #, "test"=list(exe="testfdsafsda")
                      )

cerr("------------------------------")
cerr("Required software:")
for (desc in names(software.list)) {
     details <- software.list[[desc]]
     cerr(desc)
     SPACER <- "     "
     cerr(SPACER, "Executable name: ", details$exe, " (version ", details$version, ")")
     cerr(SPACER, SPACER, "To install: ", details$install)
     if (!is.null(details$install_comment)) { cerr(SPACER, SPACER, "            (", details$install_comment, ")") }
}
cerr("------------------------------")

for (desc in names(software.list)) {
     details <- software.list[[desc]]
     get_exe_path(details$exe, paste0("How to install: ", details$install))
}


BASEDIR <- getwd()

ALIGN_CMD_FILE <- "D1a_alignment_commands_autogenerated.txt"

SHOULD_USE_FAKE_SIMULATED_DATA <- FALSE   # Should we use fake simulated data? Set to TRUE to not use real files
SIM_STATUS                     <- ifelse(SHOULD_USE_FAKE_SIMULATED_DATA, "--FAKE_SIMULATED_DATA--", "")

ERR_LOG_FILE  <- file.path(BASEDIR, "Z.errors.txt") # <-- print the GLOBAL_ERR variable to this file afterward

FAILED_IDS    <- list()
OK_IDS        <- list()
GEM_RAM_GB    <- 25  # in gigabytes. Sometimes crsahes if it's < 10

GLOBAL_DRY_RUN                <- TRUE
SAMTOOLS_EXE                  <- get_exe_path("samtools")
TOPHAT_N_THREADS              <- 4
BOWTIE_N_THREADS              <- 4
SAMTOOLS_N_THREADS            <- 4
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

clear_file  <- function(filename) {
     if (!is.null(filename)) {
          connection <- file(filename, open="w")
          close(connection)
          stopifnot(file.exists(filename))
     } else {
          # do nothing
     }
     invisible(NULL) # prevent R from printing NULL to the console if there is no action taken
}

append_commands <- function(filename, lines, commented_out=FALSE) {
     # If filename is null, writes to STDOUT
     if (commented_out) { lines <- paste("### ", lines) } # should we actually print these as COMMENTED OUT lines?
     if (is.null(filename)) {
          writeLines(lines, stdout())
     } else {
          connection <- file(filename, open="a")
          writeLines(lines, connection);
          close(connection)
     }
     invisible(NULL) # prevent R from printing NULL to the console if there is no action taken
}

read_sample_info_from_file <- function(filename) {
     if (!file.exists(filename)) { stop(paste0("[ERROR: Missing input sample data file]: The input sample data file could NOT be found. Please make sure you specified the FULL PATH to it if you're running this on a cluster, and double check the '--info=FILENAME' option on the command line. The inaccessible file was specified at the following (seemingly nonexistent or inaccessible to this user) path: ", filename)) }
     REQUIRED_COLUMN_NAMES <- c("species", "file1", "file2", "proj_type", "expr_group", "chip_input", "pubtype") # These must all be the headers in 'filename'. Additional columns are ALLOWED.
     statusPrint("Reading the 'all' data matrix from input file ", filename)
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
     die_if_file_missing(fq1, "[ERR]: Input FASTQ file ", fq1, " (forward (pair 1) and/or unpaired) seems to be missing!")
     (is.nothing(fq2) || die_if_file_missing(fq2, "[ERR]: Input FASTQ file ", fq2, " (#2, reverse pair file) seems to be missing!"))

     isPaired <- !is.nothing(fq2)
     innerDistStr <- ifelse(isPaired, " --mate-inner-dist=150 ", " ")

     btie   <- agw_get_annotation(species, "bowtie_index")
     ###fa   <- agw_get_annotation(species, "fasta")

     all_reads_bam         <- file.path(outdir, "all.bam")
     filtered_final_bam    <- file.path(outdir, "accepted_hits.bam")

     cmd <- paste0("mkdir -p ", outdir)
     if (aligner %in% c("tophat", "tophat2")) {
          gtf      <- agw_get_annotation(species, "gtf")    # <-- only actually used by tophat
          pair2Str <- ifelse(isPaired, yes=paste0(" ", fq2), no="") # might be blank, if unpaired input
          cmd <- paste0(cmd, "\n"
                      , get_exe_path('tophat'), " -o ", outdir
                      , " --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search "
                      , " --GTF=", gtf, " --num-threads=", TOPHAT_N_THREADS, " ", innerDistStr, " ", btie, " ", fq1, " ", pair2Str)
          
     } else if (aligner %in% c("bowtie","bowtie2")) {
          inputFqString <- ifelse(isPaired, yes=paste0(" -1 ", fq1," -2 ", fq2), no=paste0(" -U ", fq1))
          cmd <- paste0(cmd, "\n"
                      , get_exe_path('bowtie2'), " --threads=", BOWTIE_N_THREADS, " -x ", btie, " ", inputFqString
                      , " | " , SAMTOOLS_EXE, " view -@ ", SAMTOOLS_N_THREADS, " -b -                            > ", all_reads_bam
                      , " && ", SAMTOOLS_EXE, " view -@ ", SAMTOOLS_N_THREADS, " -b -F 0x104 ", all_reads_bam, " > ", filtered_final_bam)
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

agw_get_gem_chipseq_cmd <- function(species, ctlBam=NA, expBam, outPrefix) {
     # inputbam means; the INPUT file for chipseq. There MIGHT possibly not be one.
     # expbam is the (required) experimental bam file. 
     threads <- 1 # gem uses a ton of threads for some reason, so 1 is even conservative
     qval    <- 2	      # sets the q-value threshold (-log10, so 2 = 0.01)
     kmin    <- 6	      # minimum kmer length. From Monkey.
     kmax    <- 13	      # maximum kmer length. From Monkey.
     genomeByChromFastaDir <- agw_get_annotation(species, "fa_by_chrom")
     chrSizesFile          <- agw_get_annotation(species, "chr_sizes")

     if (!is.nothing(ctlBam) && !file.nonzero.exists(ctlBam)) {
          errlog("The missing BAM file we were looking for was named <<", ctlBam, ">>. It was part of experiment with output prefix <", outPrefix, ">.")
          return(paste0("echo 'failure for", ctlBam, "'"))
          #stopifnot(is.na(ctlBam) || file.exists(ctlBam))
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
                    , ifelse(is.na(ctlBam), " ", paste0(" --ctrl ", ctlBam))  # can be omitted if there is NO input file
                    , " --f SAM --sl --outBED"
                    , " --out ", outPrefix) # outprefix creates a new folder with this name, respecting existing subdirectories
     return(agw_construct_cmd(gemCmd))
}

agw_get_bcp_chipseq_cmd <- function(species, ctlBam, expBam, outPrefix) {

     outPrefix <- gsub("[.](bam|sam|bed)$", "", outPrefix, ignore.case=T, perl=T) # remove any '.bam'/sam/bed from the end. Those are not correct suffixes. We will add the suffix ourselves
     
     # There are two versions of BCP - BCP_TF and BCP_HM. The latter is for histone marks - what you are interested in. BCP takes bed files as input.
     # Note: BCP_HM has a very confusing help method: you specify "BCP_HM -h 1" for help. The "1" is mandatory.
     #This is an example command I ran on a sample of H3K4me3 data:
     #BCP_HM -1 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1.bed -2 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed -f 200 -w 200 -p 0.001 -3 GM12878_H3K4me3_Rep1_EncodeAlign_results001_HM.bed
     #Note, there are changes to this command if you have Chip-exo data. You will need to replace  "Read_Distribution_default.txt" by "Read_Distribution_ChIP-exo.txt" and will need to add --smooth 3 option to the commas.
     # -Reuben
     #BCP_HM -1 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1.bed
     #       -2 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed
     #       -f 200 -w 200 -p 0.001 -3 GM12878_H3K4me3_Rep1_EncodeAlign_results001_HM.bed
     
     # These are intermediate files. The problem is that bcp cannot deal with BAM files directly, so we have to convert them.
     bed_exp <- "/04a.mapping/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1.bed" # "The ChIP-seq data set you want to input."
     bed_ctl <- "/wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed" # "The control/input data set you want to input."
     
     bam2bed_exp_cmd <- agw_construct_cmd(get_exe_path("bam2bed"), " --do-not-sort ", " < ", ctlBam, " > ", bed_exp)
     bam2bed_ctl_cmd <- agw_construct_cmd(get_exe_path("bam2bed"), " --do-not-sort ", " < ", expBam, " > ", bed_ctl)
     
     fragment_size  <- 200  # "The fragment size to which we extend the reads in pre-processing data step. Default:200bp."
     window_size    <- 200  # "The window size we apply the adjcaent window in pre-processing data step. Default:200bp."
     p_cutoff       <- 1e-3 # The p_value you want to use for remove false positive based on control data. Range:1e-2---1e-6, default is 1e-3.
     bed_out        <- paste0(outPrefix, ".bed") # The results data set with 5 columns you want to output.
     
     bcp_cmd <- agw_construct_cmd(get_exe_path("BCP_HM")
                               , " -1 ", bed_exp
                               , " -2 ", bed_ctl
                               , " -f ", fragment_size
                               , " -w ", window_size
                               , " -p ", p_cutoff
                               , " -3 ", bed_out)

     allCommands <- paste(c(bam2bed_exp_cmd, bam2bed_exp_cmd, bcp_cmd), collapse="\n") # it's three separate lines!
     return(allCommands)
}

agw_get_bam_path <- function(experiment_id, sample_id) {
     bam_path <- file.path(BAM_DIR, experiment_id, sample_id, paste0("accepted_hits.bam"))
     return(bam_path)
}

agw_get_count_path <- function(experiment_id, sample_id) {
     count_path <- file.path(COUNT_DIR, experiment_id, sample_id, paste0("htseq_count.",experiment_id,"-",sample_id, ".counts.txt"))
     return(list("path"=count_path, "dir"=file.path(COUNT_DIR, experiment_id, sample_id)))
}

agw_get_gem_output_path <- function(experiment_id, sample_id) {
     parentDir <- file.path(OUT_GEM_DIR, experiment_id, sample_id)
     prefixBase <- "gem"
     return(list("prefix"=file.path(parentDir, prefixBase)
               , "dir"=parentDir
               , "created_file"=file.path(parentDir, prefixBase, paste0(prefixBase,"_result.htm"))))
}

agw_get_bcp_output_path_for_sample <- function(experiment_id, sample_id) {
     parentDir <- file.path(OUT_BCP_DIR, experiment_id, sample_id)
     return(list("prefix"=paste0(parentDir, "/")
               , "dir"=parentDir))
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
          remember_exid_failure(exid, "No replicates, skip this comparison FOR NOW ONLY.")
     } else if (!all(sapply(count_filenames.vec, file.nonzero.exists))) {
          remember_exid_failure(exid, "EDGER IS MISSING SOME HTSEQ INPUT FILES for that experiment ID so we are skipping this differential expression.")
     } else if (file.exists(edgerOut)) {
          errlog(exid, " Looks like the edgeR out ALREADY EXISTS for experiment ID ", exid, ", in the location <", edgerOut, "> so we will not be re-running it.")
          OK_IDS[[exid]] <- "OK" # we (probably) did this properly
     } else if (length(count_filenames.vec) != length(df$expr_group)) {
          remember_exid_failure(exid, paste0("missing input files for experiment ", exid, ": The number of 'count filenames' is NOT equal to the length of the group vector. This means some experimental groups were not aligned and/or processed downstream."))
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
               remember_exid_failure(exid, "Warning: did not find any group comparisons for this experiment! This may be a programming error (?) or it may be due to unusual input data.")
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
     bam    <- agw_get_bam_path(experiment_id=exid, sample_id=sid)
     gtf    <- agw_get_annotation(species, "gtf"); stopifnot(file.nonzero.exists(gtf))
     ccc    <- agw_get_count_path(experiment_id=exid, sample_id=sid)$"path"
     cccDir <- agw_get_count_path(experiment_id=exid, sample_id=sid)$"dir"
     cerr("Validating that counts file <", ccc, "> exists for sample ", sid)
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

agw_handle_chip_peak_calls <- function(exid, sid, species, df_for_this_exid, cmdfile) {
     if (agw_sample_is_an_input(exid, sid, df_for_this_exid)) { # WE don't run GEM or BCP on an input sample by itself!
          verboseStatusPrint("No need to peak-call the following INPUT file for ChIP-seq (", exid, " / ", sid, ")")
          return(invisible(NULL)) # return early...
     }

     peak_type <- "TF" # or "HISTONE" or "NON-TF"

     
     inputBam  <- agw_get_chip_input_bam_for_sample(exid, sid, df_for_this_exid)
     expBam    <- agw_get_bam_path(exid, sid)
     
     if (peak_type == "TF") {
          ggg_gem          <- agw_get_gem_output_path(exid, sid)
          chip_dir_to_make <- ggg_gem$dir
          #if (is.na(inputBam)) {
          #     msg <- paste0("[", exid, " FAILURE]: Failure to find a matching INPUT sample for experiment ", exid, " / sample ", sid, " -- not running GEM")
          #     errlog(msg)
          #     gemCmd <- paste0("# ", msg)
          #} else {
          should_comment_out_since_outputs_exist <- file.nonzero.exists(ggg_gem$created_file) # if the output file ALREADY EXISTS, then comment out the command
          peak_call_cmd <- agw_get_gem_chipseq_cmd(species, ctlBam=inputBam, expBam=expBam, outPrefix=ggg_gem$prefix)

          
     } else if (peak_type == "HISTONE" || peak_type == "NON-TF") {
          ggg_bcp          <- agw_get_bcp_output_path_for_sample(exid, sid)
          chip_dir_to_make <- ggg_bcp$dir
          
          if ("OUTPUT ALREADY EXISTS" == TRUE) {
               # TODO: maybe don't re-run in this case 
               warning("this is a placeholder--can optimize this later")
          } else {

          }
          should_comment_out_since_outputs_exist <- FALSE          
          peak_call_cmd <- agw_get_bcp_chipseq_cmd(species=species, ctlBam=inputBam, expBam=expBam, outPrefix=ggg_bcp$prefix)
          
     } else {
          stop("ERROR: Unknown peak type. We can only hadnle 'TF' and 'HISTONE' right now.")
     }

     append_commands(cmdfile, agw_construct_cmd("mkdir -p ", chip_dir_to_make), commented_out=should_comment_out_since_outputs_exist)
     append_commands(cmdfile, agw_construct_cmd(peak_call_cmd)                , commented_out=should_comment_out_since_outputs_exist)
}

# =====================  DONE WITH FUNCTION DEFINITIONS ==========================

#MARK_DUP_JAR = "/data/applications/picard/picard-tools-1.114/MarkDuplicates.jar";
if (opt$type %in% c("CHIP")) {
     BAM_DIR     <- file.path(BASEDIR, "D2b_bowtie_output")
     ALIGNED_DIR <- file.path(getwd(), "D2b_bowtie_output")
     ALIGNER     <- "bowtie2"
} else if (opt$type %in% c("RNA")) {
     BAM_DIR     <- file.path(BASEDIR, "D2a_tophat_output")
     ALIGNED_DIR <- file.path(getwd(), "D2a_tophat_output")
     ALIGNER     <- "tophat"
} else {
     stop(paste0("invalid assay type. On the command line, the assay type was specified as '--type=", opt$type, "' -- however, the only valid types currently are 'RNA' and 'CHIP'. Please specify one of those valid options!"))
}

alldat      <- read_sample_info_from_file(opt$info)
#for (exid in "24R") { # unique(dat[,"exid"])) {
counter <- 0
max_counter <- length(get_unique_experiment_ids_from_b2frame(alldat))
for (exid in get_unique_experiment_ids_from_b2frame(alldat)) {
     counter <- counter+1

     file_safe_exid    <- make.names(gsub("[:;, ]", "_", exid)) # make this exid safe for use in a filename
     zero_padded_count <- sprintf(paste0("%0", floor(log10(max_counter)+1),"d"), counter) # <-- try to calculate the right number of zeros to pad with. e.g. if the max is "999", then pad to a total of at least 3 digits (001, 002, etc...)
     xout              <- paste0(opt$out, ".", zero_padded_count, ".", file_safe_exid, ".cmd.txt")     # <-- the output filename

     clear_file(xout) # Clear out the output file before we write to it, in case it already exists...
     append_commands(xout, paste0("echo 'Handling the experiment ID ", file_safe_exid, "'. This file contains the commands for ONLY that specific experiment ID"))
     
     subdat                      <- get_b2frame_for_exid_from_b2frame(alldat, exid)
     sample_ids_for_exid.vec     <- get_sids_for_exid_from_b2frame(df=subdat, experiment_id=exid)
     species                     <- get_species_from_b2frame(subdat)
     cerr("Handling experiment ", exid, " (", species, "). The full set of data for this experiment is:"); cerr(subdat)
     f1.vec <- subdat[["file1"]]
     f2.vec <- subdat[["file2"]]

     # Handle the alignments
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
          append_commands(xout, align_cmd) # align it!
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
          count_filenames.vec <- c() # Get real data
          for (sid in sample_ids_for_exid.vec) {
               bam <- agw_get_bam_path(experiment_id=exid, sample_id=sid)
               if (!file.nonzero.exists(bam)) {
                    errlog("[", exid, " FAILURE]: [MISSING BAM FILE in experiment ID ", exid, ", sample ID ", sid, ". Specifically, file <", bam, "> does not exist.")
                    exid_is_missing_a_sample <- TRUE
               } else {
                    if      (opt$type %in% c("RNA") ) { agw_handle_rna_counting_per_feature(exid=exid, sid=sid, species=species,                          cmdfile=xout) }
                    else if (opt$type %in% c("CHIP")) {          agw_handle_chip_peak_calls(exid=exid, sid=sid, species=species, df_for_this_exid=subdat, cmdfile=xout) }
                    else                              { stop(paste0("Unrecognized assay type, specifically: ", opt$type)) }
               }
          }
          if (exid_is_missing_a_sample) {
               remember_exid_failure(exid, paste0("[SKIPPING all remaining analysis steps for experiment ID ", exid, " due to a missing file issue earlier. See the logs for the exact sample ID that was missing."))
          } else {
               if (opt$type == "RNASEQ") {
                    agw_handle_rna_diff_expression(exid, sid, subdat, outdir=file.path(OUT_EDGER_DIR, exid), cmdfile=xout)
               }
          }
          #browser()
     }
}

write_global_errors_to_file(filename=ERR_LOG_FILE)

if (length(GLOBAL_ERR) > 0) {
     cerr("[ERROR] We encountered the following ", length(GLOBAL_ERR), " error messages:")
     cerr(paste(GLOBAL_ERR, collapse="\n"))
     cerr("---------------")
}
if (length(FAILED_IDS) > 0) {
     cerr("[ERROR] Failure for the following ", length(FAILED_IDS), " experiment IDs:")
     cerr(paste(names(FAILED_IDS), collapse=", "))
 } else {
     cerr("[OK] Commands were successfully written for every ID.")
}
cerr("[OK] with the following ", length(OK_IDS), " experiment IDs:")
cerr(paste(names(OK_IDS), collapse=", "))
cerr("---------------")


cerr("You should see a bunch of output files with the following prefix/path (one file per experiment):")
cerr("      Look for files matching this pattern ---> ls ", opt$out, "*")


'''