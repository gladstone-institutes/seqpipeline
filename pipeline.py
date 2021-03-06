#!/usr/bin/env python2

'''
Run this script with --help to see some examples use cases and command line option descriptions.

This is a python2.x script that generates a UNIX shell script, which then invokes external tools to handle:
   * Alignment of FASTQ -> BAM
   * RNA-SEQ: Counting of reads (using Subread FeatureCounts) & differential expression (via an edgeR script in R)
   * CHIP-SEQ: Peak-calling using GEM and/or BCP

It requires at least python2.7 due to the use of the 'argparse' library.

Note that a lot of external software is used; check the top of the source file for the (hard-coded) paths to
programs like java (the JAVA_PATH variable), GEM (GEM_JAR_PATH), etc...
'''

from __future__ import print_function # sets print("...") to have parens like in python 3
from __future__ import division # defaults to floating point division. No more "1/2 == 0"
import sys
import pdb #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques
import os.path
import argparse # requires python 2.7
#import re

global_file_paths_already_aligned = dict() # Remember which files we are aligning, so we don't try to double-align a file and overwrite it
global_annot = dict()
# ==================================================================================
# ==================================================================================
# ======= Paths to the programs we need and various other hard-coded things.
# ======= Programs are currently assumed to be on the user's $PATH ========
global_annot['hg19'   ] = "/data/info/genome/hg19_ensembl_igenome_with_chr/hg19.chr"
global_annot['mm9'    ] = "/data/info/genome/mm9_ensembl_igenome_with_chr/mm9.chr"
global_annot['danRer7'] = "/data/info/genome/danRer7_ensembl_igenome_ucsc/danRer7_ucsc"
global_annot['galGal4'] = "/data/info/genome/galGal4_ensembl_igenome_ucsc/galGal4_ucsc"
TOPHAT_PATH             = "tophat"
BOWTIE_PATH             = "bowtie2"
BAM2BED_PATH            = "bam2bed"
SAMTOOLS_PATH           = "samtools"
HTSEQ_COUNT_PATH        = "htseq-count"
RSCRIPT_PATH            = "Rscript" # <-- part of the default R install
EDGER_SCRIPT            = "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/run_edgeR_diff_expr.R" # this is a CUSTOM SCRIPT that we run
GEM_JAR_PATH            = "/data/applications/2015_06/bin/gem.jar" # peak caller for TF / narrow peaks with motifs
BCP_HM_PATH             = "BCP_HM" # Peak caller, good on broad peaks. Not the best for motifs, but still pretty good (see Reuben et al's paper)
JAVA_PATH               = "java" # needed for GEM.jar, the GEM ChIP peak calling program. GEM is best at narrow peaks with particular sequence motifs.

DEFAULT_ALIGN_BOWTIE_DIR    = "Pipeline_01b_Align_Bowtie_Dir"
DEFAULT_ALIGN_TOPHAT_DIR    = "Pipeline_01t_Align_Tophat_Dir"
DEFAULT_HTSEQ_COUNT_DIR     = "Pipeline_02_HTSeq_Count_Dir"
DEFAULT_EDGER_DIFF_EXPR_DIR = "Pipeline_03_EdgeR_Diff_Expr_Dir"
DEFAULT_BCP_PEAK_DIR        = "Pipeline_05b_BCP_Peak_Dir"
DEFAULT_GEM_PEAK_DIR        = "Pipeline_05g_GEM_Peak_Dir"

CHIP_BCP_STR = "BCP"
CHIP_GEM_STR = "GEM"

GEM_READ_DIST_FILE = "/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/A-2016-08-August/github_seqpipeline/resources/GEM_Read_Distribution_default.txt"
print("Warning: note that the gem read distribution file is a HARD CODED file path right now!")

GEM_RAM_GB        = 25  # in gigabytes. Sometimes crashes if it's < 10

TOPHAT_N_THREADS   = 4 # Somewhat arbitrary
BOWTIE_N_THREADS   = 4
SAMTOOLS_N_THREADS = BOWTIE_N_THREADS
# You probably won't need to change anything below... hopefully
# ==================================================================================
# ==================================================================================

EXIT_CODE_IF_MISSING_FILE        = 49   # just some arbitrary non-zero exit code number
LITERAL_BACKSLASH                = "\\" # <-- reduces to just one backslash
LITERAL_DQUOTE                   = "\""
ASSAY_RNA  = 1 # <-- Arbitrary. could theoretically be an Enum
ASSAY_CHIP = 2
DEFAULT_OUTPUT_BASEDIR         = "./" # write files to the current directory if no one says otherwise
DELIM_FOR_ARGS_TO_EDGER_SCRIPT = ',' # should be a comma normally. Do not change this unless the edgeR R script (WHICH IS A SEPARATE FILE!!!!) also changes. Do not change this without changing the R file too!




opt = None # <-- "opt" is a global variable that stores the cmd line options. This should really be the ONLY non-constant global!

def enquote(s):
    # should we check to see if s has quotes in it already? Maybe.
    return("\"" + str(s) + "\"")

def verbosePrint(m):
    if (opt.verbose): printStderr(m)
    return

def progressPrint(m):
    printStderr(m) # print even if we are NOT in verbose mode
    return

def argAssert(thing, message): # for things that involve user input somehow
    if not thing: dieBadArgs(message)

def xAssert(thing, message="No specific message"): # for detecting programming bugs only
    if not isinstance(message, str):
        raise "wrong message argument--must be a string"
    if not thing: raise(Exception(message))

def isSameLenOrNone(x, y):
    return (x is None) or (y is None) or (len(x) == len(y))

def printStderr(*args, **keyargs):
    print(*args, file=sys.stderr, **keyargs)

def pathWithBase(basedir, path):
    return path if basedir is None else os.path.join(basedir, path)

def splitFileList(string_of_files, delim, basedir=None):
    # Splits up a string
    if string_of_files is None:
        return None
    assertType(string_of_files, str)
    assertType(delim, str)
    list_of_files = [pathWithBase(basedir, fff) for fff in string_of_files.split(delim)] # python list comprehension---add the basedir to the beginning of each file path
    return list_of_files

def getListItemUnlessNone(theList, index):
    return None if (theList is None) else theList[index]

#def tar_up_that_annotation():
#    '''Just a function for generating the list of files to be tarred up, if we want to send this annotation to someone else. Not really production-ready code here!'''
#    for assembly in global_annot:
#        a = getAnnot(assembly)
#        for subitem in a:
#            path = a[subitem]
#            if subitem == 'bowtie_index':
#                path += "*.bt2"
#                pass
#            else:
#                if not (os.path.isfile(path) or os.path.isdir(path) or os.path.isfile(path + ".1.bt2")):
#                    print("**WARNING**: file or path above does not appear to exist! (This one: " + path + ")")
#                    pass
#                pass
#            print(path)
#            pass
#        pass
#    return

def dieBadArgs(errMsg):
    sys.exit("ERROR: Problem with the command line arguments! Specifically: " + errMsg + ". Try using '--help' to see all the possible arguments.")
    return

def handleCmdLineArgs():
    DESC = '''
A program for running a standard ChIP- or RNA-seq pipeline.
Written for Python 2.7+ by Alex Williams, 2016. (Python 2.7 is needed for the 'argparse' library)

What this program does:

1) a) First, you provide the details about an experiment and the associated samples names and input FASTQ files.
   b) Then, you run this 'pipeline.py' program. An example command is listed below in section (2)
   c) 'pipeline.py' does not actually perform any bioinformatics analysis directly---instead,
      it generates a shell script with a list of commands.
   d) You then would want to run that shell script (let's call it 'script.sh') either
      directly ('bash script.sh') or by submitting it to a compute
      cluster (e.g. 'qsub --options-here script.sh')
   e) When the script runs, it will take quite some time. The script will (slowly) generate a bunch of output
      directories, each with one step of the analysis process. The default names are:
         * Pipeline_01t_Alignment_Tophat
         * Pipeline_01b_Alignment_Bowtie
         * Pipeline_02_HTSeq_Count_Dir
         * Pipeline_03_EdgeR_Diff_Expr_Dir
         * ... etc.

Important caveats to be aware of:
   * You will need to specify FULL PATHS to everything if you are planning on running the output script on a **cluster**!

2) You run a command line invocation like:

     RNA-SEQ example:
     python2  pipeline.py  --script="script_test.sh" --outdir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/" --basedir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/" --experiment-id="Test_Experiment" --sample-ids="A1,A2,B1,B2" --rna-samples=a1.mm9.chr19.fq.gz,a2.mm9.chr19.fq.gz,b1.mm9.chr19.fq.gz,b2.mm9.chr19.fq.gz --groups=1,1,2,2 --species=mm9
        * Note that ***FULL PATHS** for 'outdir' and 'basedir' are specified since this example is for a cluster.

     ChIP-SEQ example:
python2  pipeline.py \\
--script="chip_test.sh" \\
 --outdir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/" \\
--basedir="/data/projects/kp-600-b2b-osono-data-pipeline-run-feb-16/B-2016-11-November/test_data/" \\
--experiment-id="Test_Experiment_ChIP" \\
--sample-ids="ALPHA,BETA" \\
--groups=1,2 \\
--chip-samples=a1.mm9.chr19.fq.gz,b1.mm9.chr19.fq.gz \\
 --chip-inputs=a2.mm9.chr19.fq.gz,b2.mm9.chr19.fq.gz \\
--species=mm9 \\


3) Details: the '--outdir' is the location where you want the script to generate OUTPUT
   Details: the '--script' is the name of the actual script to generate. One output script is generated per invocation of this program.

Note that this is actual test data (a small set of reads from mouse chr19) that you can download:
      a1.mm9.chr19.fq.gz   (these reads were pre-selected to align to mm9, so 100% of them should align.)
      a2.mm9.chr19.fq.gz   (you can of course align them to any mouse genome build, however)
      a3.mm9.chr19.fq.gz
      b1.mm9.chr19.fq.gz
      b2.mm9.chr19.fq.gz
      b3.mm9.chr19.fq.gz

3) The output is a shell script that you can then submit on your cluster.
'''

    EPILOG = '''
Version history: (none yet)
'''

    # Check for certain files that are maybe-required first...
    for fff in [EDGER_SCRIPT, GEM_JAR_PATH, global_annot['hg19'], global_annot['mm9'], global_annot['danRer7'], global_annot['galGal4'] ]:
        if not (os.path.isfile(fff) or os.path.isfile(fff + ".1.bt2")):
            print("**WARNING**: Note that we can't seem to find the following (probably hard-coded) path <" + fff + ">. Probably this is because it's HARD CODED on the filesystem right now!")
            #sys.sleep(1)
            pass
        pass
    
    parser = argparse.ArgumentParser(description=DESC, epilog=EPILOG, add_help=True, formatter_class=argparse.RawTextHelpFormatter)
    #, formatter_class=argparse.RawTextHelpFormatter)

    #parser.add_argument('--sup', '-s', dest='sup', default=None, metavar="\"UMI_LEFT,LINKER_LEFT,UMI_RIGHT,LINKER_RIGHT\"", type=str, help="Default: if unspecified (or set as --sup=GUESS), we will just guess the UMI lengths by looking at the consensus sequence of the first 1000 reads (and the linker wil be assumed to be ZERO bases--it will be included as part of the 'payload' sequence). 'Sup' stands for 'Supplementary' lengths (note: numbers, NOT actual sequence) for the UMI_LEFT, LINKER_LEFT, UMI_RIGHT, and LINKER_RIGHT.\nExample: --sup=\"8,2,2,5\"\n   * for a sample with an 8+5 UMI and linkers of length 2.")

    parser.add_argument('-v',  '--verbose', dest='verbose', action='store_true', help="Verbose debugging")
    parser.add_argument('--groups', '-g', type=str, default=None, dest='groups', help="Specify experiment groups (numeric). Must be the same length as the number of files, and must be sequentially numbered with the lowest group being 1. Can be 'NA' or blank for samples that are not associated with a group (like ChIP inputs). OK example: --groups='1,1,2,3,,,4,4' (Note the two blank samples!).")

    parser.add_argument('--basedir', '-b', type=str, default=None, dest='basedir', help="Optional. Specify a base directory to prepend to all paths.!")

    parser.add_argument('--chip-samples', '--c1', type=str, default=None, dest='chip_samples', help="ChIP-seq only: Specify ChIP-seq sample files, comma-delimited. Must match the order of the corresponding --chip_inputs files!")
    parser.add_argument('--chip-inputs' , '--p1', type=str, default=None, dest='chip_inputs', help="ChIP-seq only: Specify ChIP-seq input files, comma-delimited. Must match the order of the corresponding --chip_samples files!")
    parser.add_argument('--peak-caller' , '--pc', type=str, default=None, dest='peak_caller', help="ChIP-seq only: Specify which peak caller to use. Must be either 'gem' or 'bcp' (case-insensitive). No default.!")
    parser.add_argument('--rna-samples', '--r1', type=str, default=None, dest='rna_samples', help="RNA-seq only: Specify RNA-seq input files, comma-delimited.")

    parser.add_argument('--sample-mates', '--c2', '--r2', type=str, default=None, dest='sample_mates', help="The mate pair associatd with either the CHIP or RNA samples. Comma-delimited just like the normal samples. Example: --r1=DrugX.fq,DrugY.fq --r2=DrugXpair2.fq,DrugYpair2.fq")
    parser.add_argument('--input-mates', '--p2', type=str, default=None, dest='input_mates', help="The mate pair associated with the second end of a pair in a CHIP sample. Comma-delimited just like the normal samples.")

    parser.add_argument('--species', '-s', type=str, default=None, dest='species', help="Required: The short name for the species.")

    parser.add_argument('--sample-ids', '--sids', type=str, default=None, dest='sample_ids', help="Required: The short sample ID names, comma-delimited.")

    parser.add_argument('--experiment-id', '--eids', type=str, default=None, dest='exper_id', help="Required: The single experiment ID (name).")
    
    parser.add_argument('--delim', type=str, default=',', dest='delim', help="Default file/group delimiter. Normally should be a comma, unless you're doing something very unusual.")

    parser.add_argument('--script', type=str, default=None, dest='scriptname', help="The output script file to write, e.g. --script='/path/to/script.sh' ")
    parser.add_argument('--outdir', '-o'   , type=str, default=None, dest='outdir', help="Optional. Specify a base output directory. Must already exists. Default is the current directory.")

    parser.add_argument('--debug-tar-annotations', dest="debug_tar_annotations", action='store_true', help="A debug command only for setting up this software. Generates a tar command that can be used to save the annotation data in one place.")

    global opt # <-- critical, since we modify (global variable) "got" below. Do not remove this!
    opt = parser.parse_args()

    #if (opt.debug_tar_annotations):
    #    tar_up_that_annotation()
    #    print("Exiting early since this is just a command for tarring files.")
    #    return

    assay = None

    argAssert(opt.species is not None, "A species must be specified! It must be one species for ALL samples currently. Example species: 'mm9' or 'danRer7' or 'hg19'.")
    argAssert(opt.rna_samples is None or (opt.chip_samples is None and opt.chip_inputs is None), "You cannot specify BOTH RNA-seq and ALSO ChIP-seq at the same time! Pick either RNA-seq or ChIP-seq")
    argAssert(opt.groups is not None, "You must specify the groups that each sample belongs to, or either 'NA' or a blank value for files that are not direclty associated with a group (such as ChIP inputs).")

    opt.outdir = DEFAULT_OUTPUT_BASEDIR if opt.outdir is None else opt.outdir

    if (opt.rna_samples is not None):
        assay = ASSAY_RNA
        samp1 = splitFileList(opt.rna_samples , delim=opt.delim, basedir=opt.basedir)
        samp2 = splitFileList(opt.sample_mates, delim=opt.delim, basedir=opt.basedir)
        inp1 = None # Rna-seq doesn't have 'input' files
        inp2 = None # Rna-seq doesn't have 'input' files
    elif (opt.chip_samples is not None):
        argAssert(opt.chip_inputs is not None, "If you specify --chip-samples, you have to also specify --chip-inputs, but it appears that was NOT specified.")
        assay = ASSAY_CHIP
        samp1 = splitFileList(opt.chip_samples, delim=opt.delim, basedir=opt.basedir)
        samp2 = splitFileList(opt.sample_mates, delim=opt.delim, basedir=opt.basedir)

        inp1  = splitFileList(opt.chip_inputs, delim=opt.delim, basedir=opt.basedir)
        inp2  = splitFileList(opt.input_mates, delim=opt.delim, basedir=opt.basedir)
        # Apparently we're running a ChIP-seq project
        argAssert(len(inp1) == len(samp1), "Somehow your --chip-inputs (which is length " + str(len(inp1)) + ") was a different length from the --chip-samples (which is length " + str(len(samp1)) + "! Fix this.")
        argAssert(isSameLenOrNone(inp1, inp2), "Forward and reverse MATE PAIRS must have the same number of files!")
    else:
        dieBadArgs("You have to AT LEAST have an RNA-seq sample to process or a ChIP-seq sample to process. Specify at least one!")
        pass

    argAssert(opt.sample_ids is not None, "Sample IDs must be specified as a comma-delimited list. For example: --sample-ids='S123,S191,S990'. One per fastq file pair.")
    sample_id_list = opt.sample_ids.split(opt.delim)
    argAssert(len(samp1) == len(sample_id_list), "The number of sample IDs and the number of input fastq files is DIFFERENT, which indicate a problem")

    argAssert(opt.exper_id is not None, "An experiment ID (--experiment-id='Something') is required, but one was not specified")
    argAssert(isSameLenOrNone(samp1, samp2), "Forward and reverse MATE PAIRS must have the same number of files!")

    groups = opt.groups.split(opt.delim)
    
    are_groups_ok = [(g.isdigit() and int(g) > 0) or g == '' or g is None or g.upper() == 'NA' for g in groups]
    argAssert(all(are_groups_ok), "Some of your input groups were not in the valid set of being either a number >= 1, the literal text 'NA', or totally blank. Fix the input groups!")

    argAssert(len(groups) == len(samp1), "Somehow your '--groups' (length " + str(len(groups)) + ") was not equal in length to your specified number of sample files (which had a length of " + str(len(samp1)) + ")!")
    verbosePrint("Groups: " + str(groups))

    # Finally, actually generate the commands...

    argAssert(opt.species in global_annot, "Your specified --species (a genome assembly) was not in our recognized list, unfortunately. Check your capitalization and spaces.")

    xcmd = OurScript()

    base_bam_subdirectory_only = DEFAULT_ALIGN_TOPHAT_DIR if (assay == ASSAY_RNA) else DEFAULT_ALIGN_BOWTIE_DIR # this changes based on whether we're using Tophat or Bowtie

    eee = Experiment(expName=opt.exper_id, species=opt.species, sid_list=sample_id_list
                     , fq1_list=samp1, fq2_list=samp2 # Normal fastq files for the samples
                     , inp1_list=inp1, inp2_list=inp2 # Inputs for ChIPseq
                     ,      base_bam_dir=pathWithBase(opt.outdir, base_bam_subdirectory_only)
                     ,    base_count_dir=pathWithBase(opt.outdir, DEFAULT_HTSEQ_COUNT_DIR)
                     ,    base_edger_dir=pathWithBase(opt.outdir, DEFAULT_EDGER_DIFF_EXPR_DIR)
                     , base_bcp_peak_dir=pathWithBase(opt.outdir, DEFAULT_BCP_PEAK_DIR)
                     , base_gem_peak_dir=pathWithBase(opt.outdir, DEFAULT_GEM_PEAK_DIR)
                 )

    xcmd.append("echo 'First, we are checking to make sure all the FASTQ input files exist...'")
    xcmd.appendCheckForRequiredFiles(file_tup=eee.getTupleOfAllFqFiles())
    xcmd.append("")

    if assay == ASSAY_RNA:
        xcmd.append("###### Handling RNA-seq here ")
        handleRNA(groups=groups, script_obj=xcmd, exper=eee)
    elif assay == ASSAY_CHIP:
        xcmd.append("###### Handling ChIP-seq here ")
        argAssert(opt.peak_caller is not None and opt.peak_caller.upper() in ("BCP", "GEM"), "If you are running chIP-seq peak calls, then be sure to specify the peak caller as either BCP or GEM, using --peak-caller=BCP or --peak-caller=GEM .")
        handleCHIP(groups=groups, script_obj=xcmd, exper=eee, caller=opt.peak_caller)
    else:
        sys.exit("Something went wrong (programming bug)---this assay type is not recognized!")
        pass

    xcmd.writeToDisk(scriptname=opt.scriptname)
    
    return  # End of command-line-reading function


def assertType(x, oktype, none_ok=False):
    if none_ok and x is None:
        return # this is ok!
    if (not isinstance(oktype, type)):
        raise Exception("The 'oktype' argument must be a 'type'.")
    # Require that everything in the list_of_things is one of the "oktype"
    xAssert(isinstance(x, oktype), "We expected the variable to be of type "+str(oktype)+", but instead it was a "+str(type(x))+"")

def withoutNone(aaa):
    if aaa is None: # turn a 'None' input into an empty list
        aaa = ()
    xAssert(isinstance(aaa, (tuple, list)), "problem! incorrect argument to 'withoutNone', which must be a tuple or list, or None.")
    return tuple([x for x in aaa if x is not None]) # remove any 'None' elements from the list

def maybeNoneDict(keys, values):
    if values is None:
        return dict.fromkeys(keys) # keys are normal, but values are all 'None' if a single "none" was passed in as 'values'
    else:
        return dict(zip(keys, values)) # just normally initialize a python dictionary with keys & matching values



class Experiment(object):
    '''A class that keeps track of all our samples for an experiment.
    This is the class that can answer questions like "what are the FASTQ files associated with a specific sample?"
    Also this is where we remember all the sample names.
    This is basically where ALL THE DATA about an experiment is stored.'''
    def __init__(self, expName, species, sid_list, fq1_list, fq2_list, inp1_list, inp2_list, base_bam_dir, base_count_dir, base_edger_dir, base_bcp_peak_dir, base_gem_peak_dir):
        xAssert(isinstance(sid_list, (list, tuple)))
        xAssert(isinstance(expName, (str)))
        self.expName = expName
        self.sids    = sid_list # list of sample IDs
        self.species = {key:species for key in self.sids} # dict
        self.f1 = maybeNoneDict(self.sids, fq1_list)
        self.f2 = maybeNoneDict(self.sids, fq2_list) # f2 = mate pair corresponding to f1
        self.i1 = maybeNoneDict(self.sids, inp1_list) # <-- INPUT (control) for ChIP. Not relevant for RNA-seq
        self.i2 = maybeNoneDict(self.sids, inp2_list) # mate pair of i1, if any ("None" if there isn't one)
        self.base_bam_dir   = base_bam_dir
        self.base_count_dir = base_count_dir
        self.base_edger_dir = base_edger_dir
        self.base_bcp_peak_dir = base_bcp_peak_dir
        self.base_gem_peak_dir = base_gem_peak_dir
        pass

    def sampleHasChIPControl(self, sampName):
        return self.i1[sampName] is not None
    
    def getBamDirForSample(self, sampName):
        # Example:  my_bam_files/EXPERIMENT_27/SAMPLE_18/
        return os.path.join(self.base_bam_dir, self.expName, sampName) # Note that this makes TWO subdirectories!

    def getBamDirForChIPControlOfSample(self, sampName):
        # Example:  my_bam_files/EXPERIMENT_27/SAMPLE_18/
        CHIP_CONTROL_SUFFIX="_CHIP_CONTROL"
        return os.path.join(self.base_bam_dir, self.expName, sampName+CHIP_CONTROL_SUFFIX) # Note that this makes TWO subdirectories!
                        
    def getSampleNames(self): # returns: list
        return self.sids

    def getFqPairForSample(self, sampName): # returns (2-element tuple)
        return (self.f1[sampName], self.f2[sampName])

    def getFqControlsForSample(self, sampName): # returns (2-element tuple)
        return (self.i1[sampName], self.i2[sampName]) # The inputs/controls for ChIPSeq

    def getTupleOfAllFqFiles(self): # Just a list of all *unique* fastq files, not in any particularly significant order. No "none"s
        return withoutNone(tuple(set(self.f1.values()).union(self.f2.values()).union(self.i1.values()).union(self.i2.values())))

    def getSpeciesForSample(self, sampName): # returns: string (single assembly name)
        return self.species[sampName]

    def getAlignedBamForSample(self, sampName):  # returns: string (single file path)
        return os.path.join(self.getBamDirForSample(sampName), "accepted_hits.bam")

    def getAlignedChIPControlBamForSample(self, sampName):  # returns: string (single file path)
        # Returns a file path, OR rturns None if this sample doesn't have a corresponding chip input.
        return os.path.join(self.getBamDirForChIPControlOfSample(sampName), "accepted_hits.bam") if self.sampleHasChIPControl(sampName) else None

    def getCountDirForSample(self, sampName): # returns: string (a single directory full path)
        return os.path.join(self.base_count_dir, self.expName, sampName)

    def getCountFileForSample(self, sampName): # returns: string (single file path)
        return os.path.join(self.getCountDirForSample(sampName), "htseq_count."+self.expName+"."+sampName+".matrix.txt")

    def getAllCountFiles(self): # returns: **list** of file paths
        return [self.getCountFileForSample(x) for x in self.getSampleNames()] # <-- list comprehension: return a LIST of all count file paths

    def getEdgeRDiffExprDir(self): # returns: string (single full file path)
        return os.path.join(self.base_edger_dir)

    def getEdgeRDiffExprFile(self): # returns: string (single full file path)
        return os.path.join(self.getEdgeRDiffExprDir(), "edgeR." + self.expName + ".out.txt")

    def getBcpPeakDir(self, sampName): # returns: string (single full file path)
        return os.path.join(self.base_bcp_peak_dir, self.expName, sampName)

    def getBcpPeakCallBedFile(self, sampName): # returns: string (single full file path)
        return os.path.join(self.getBcpPeakDir(sampName), (sampName+"_BCP_peaks.bed"))

    def getGemPeakDir(self, sampName): # returns: string (single full file path)
        return os.path.join(self.base_gem_peak_dir, self.expName, sampName)

    def getGemOutHTMLFile(self, sampName): # returns: string (single full file path)
        return os.path.join(self.getGemPeakDir(sampName), sampName + "_result.htm")


class OurScript(object):
    """Script text for writing to a file"""
    def __init__(self):
        self.lines = ["#!/bin/bash -u"
                      , "set -e"
                      , "set -o pipefail"]
        self.lines.append("")
        self.lines.append("### You can run this script with a command like one of these:")
        self.lines.append("###                                           * bash ./this_script_name.sh")
        self.lines.append("###             Or submit it to a cluster:    * qsub --some-options ./this_script_name.sh")
        self.lines.append("")
        pass

    def append(self, text):
        self.lines.append(text) # add exactly one line
        return

    def appendMkDir(self, dirpath):
        self.append("mkdir -p " + enquote(dirpath) + "")
        return

    def appendCheckForRequiredFiles(self, file_tup=()):
        # File list can be a list of files or NONE (none just means we don't do anything)
        assertType(file_tup, tuple, none_ok=True)
        file_tup = withoutNone(file_tup) # make sure it's a [ ] list, and not just like ONE file someone passed in
        for f in file_tup:
            self.append("")
            self.append("# ===== Checking for a required file")
            self.append("if [[ ! -e " + enquote(f) + " ]]; then")
            self.append("    echo '[ERROR] Cannot find the following required file, which we expect to exist: <" + f + ">" + "'")
            self.append("    exit " + str(EXIT_CODE_IF_MISSING_FILE)) # exit with a somewhat arbitrary code number if there is a missing file
            self.append("fi")
            self.append("")
            pass
        return

    def appendCommandUnlessFilesExist(self, check_files=(), cmd=None, required_prereqs=(), required_output=() ):
        # cmd will ONLY be executed if there is one or more filenames in "check_files" that do not exist.
        # Note: check_files absolutely CANNOT be empty or None!
        check_files      = tuple(check_files)       if isinstance(check_files, str)      else check_files      # wrap the string in a tuple just in case...
        required_prereqs = tuple(required_prereqs)  if isinstance(required_prereqs, str) else required_prereqs # wrap the string in a tuple just in case...
        required_output  = tuple(required_output)   if isinstance(required_output, str)  else required_output  # wrap the string in a tuple just in case...
        assertType(cmd, str)
        assertType(check_files, tuple)
        assertType(required_prereqs, tuple)
        assertType(required_output, tuple)

        self.appendCheckForRequiredFiles(required_prereqs) # make sure any prereqs (if any) actually exist

        check_files = withoutNone(check_files) # make sure it's a [ ] list, and not just like ONE file someone passed in
        xAssert(len(check_files) >= 1)
        multiIfString = "if " + " && ".join(  [ "[[ -e " + enquote(g) + " ]]" for g in check_files ] ) # <-- python list comprehension that generates results like "if [[ -e myFile ]] && [[ -e file2 ]] && [[ -e file3 ]]"


        checked_files_str = ", ".join(check_files) # <-- python list comprehension that generates results like "if [[ -e myFile ]] && [[ -e file2 ]] && [[ -e file3 ]]"

        self.append("\n#-----------------------")
        self.append(multiIfString)
        self.append("then")
        self.append("    echo '[SKIPPING] All the specified output files (" + checked_files_str + ") already exist--so we are skipping this step' ")
        self.append("else")
        self.append("    echo '[RUNNING]...'")
        self.append(cmd)
        self.append("fi")
        self.append("#-----------------------\n")

        self.appendCheckForRequiredFiles(required_output) # make sure we generated the required files (if any)

        return
        
    def writeToDisk(self, scriptname=None):
        if (scriptname is None):
            print("\n".join(self.lines)) # just print to stdout
        else:
            scriptname = scriptname
            with open(scriptname, 'w') as fff:
                fff.write("\n".join(self.lines) + "\n")
                pass
            pass
            progressPrint("Wrote script data to the following filename: " + scriptname)
        return


def handleRNA(groups, script_obj, exper):
    xAssert(isinstance(exper, (Experiment)))
    # Step 1/3: Alignment....
    for sampName in exper.getSampleNames():
        (s1, s2) = exper.getFqPairForSample(sampName)
        generateAlignmentCmd(fastq1=s1, fastq2=s2, species=exper.getSpeciesForSample(sampName), outdir=exper.getBamDirForSample(sampName), outbam=exper.getAlignedBamForSample(sampName), aligner=TOPHAT_PATH, script_obj=script_obj)
        pass

    # STEP 2/3: Counts....
    for sampName in exper.getSampleNames():
        theGtf = getAnnot(exper.getSpeciesForSample(sampName), "gtf")
        theBam = exper.getAlignedBamForSample(sampName)
        generateHtseqCountCmd(bam=theBam, gtf=theGtf, outdir=exper.getCountDirForSample(sampName), outfile=exper.getCountFileForSample(sampName), script_obj=script_obj)
        pass
    
    # STEP 3/3: Differential expression
    generateEdgeRCmd(count_files=exper.getAllCountFiles(), groups=groups, outdir=exper.getEdgeRDiffExprDir(), outfile=exper.getEdgeRDiffExprFile(), script_obj=script_obj) # Note: this will call a custom edgeR script, which is also in this same repository.
    return

def handleCHIP(groups, script_obj, exper, caller=None):
    # Alignment...
    assertType(caller, str)
    for sampName in exper.getSampleNames():
        species = exper.getSpeciesForSample(sampName)

        (s1, s2) = exper.getFqPairForSample(sampName)
        expr_bam_dir = exper.getBamDirForSample(sampName)       
        generateAlignmentCmd(fastq1=s1, fastq2=s2, species=species, outdir=expr_bam_dir, outbam=exper.getAlignedBamForSample(sampName), aligner=BOWTIE_PATH, script_obj=script_obj)

        (p1, p2) = exper.getFqControlsForSample(sampName)
        chip_control_bam_dir = exper.getBamDirForChIPControlOfSample(sampName)
        generateAlignmentCmd(fastq1=p1, fastq2=p2, species=species, outdir=chip_control_bam_dir, outbam=exper.getAlignedChIPControlBamForSample(sampName), aligner=BOWTIE_PATH, script_obj=script_obj)
        # Note that although we generate redundant alignment commands for controls appearing multiple times, only one of those actually gets run, so things should be A-OK
        pass

    # Peak calling...
    for sampName in exper.getSampleNames():
        if caller.upper() == CHIP_GEM_STR.upper():  generateGemPeakCmd(exper=exper, sampName=sampName, script_obj=script_obj)
        elif caller.upper() == CHIP_BCP_STR.upper(): generateBcpPeakCmd(exper=exper, sampName=sampName, script_obj=script_obj)
        else:                          raise Exception("Programming bug -- unsupported peak caller")
        pass

    return


def getAnnot(assembly, key=None):
    # Sort of a weird wrapper to an array. Probably not the most elegant way to do this, but should be OK for our low-performance needs. Easy to change in the future!
    # 'None' is how you get the whole data structure
    xAssert(key in [None, 'gtf', 'bowtie_index', 'fasta', 'fa_by_chrom', 'chr_sizes'], "Invalid annotation requested") # sanity check
    vals = dict()
    vals['bowtie_index'] = global_annot[assembly]
    vals['gtf'         ] = global_annot[assembly] + ".gtf"
    vals['fasta'       ] = global_annot[assembly] + ".fa"
    vals['fa_by_chrom' ] = global_annot[assembly] + "_chromosome_fastas"
    vals['chr_sizes'   ] = global_annot[assembly] + ".chrom.sizes.txt"
    return vals if key is None else vals[key]

def generateHtseqCountCmd(bam, gtf, outdir, outfile, script_obj):
    sortedbam                 = bam + ".sorted_by_name.bam"
    samtools_sort_by_name_cmd = " ".join([SAMTOOLS_PATH, "sort -n", bam, ">", sortedbam])
    ht_params                 = " --format=bam --order=name --mode=intersection-nonempty --stranded=no --minaqual=10 --type=exon --idattr=gene_id "
    ht_cmd                    = " ".join([HTSEQ_COUNT_PATH, ht_params, sortedbam, gtf, ">", outfile])
    script_obj.appendMkDir(outdir) # this directory must exist before we try moving any files to it...
    script_obj.appendCommandUnlessFilesExist(check_files=(sortedbam,), cmd=samtools_sort_by_name_cmd, required_prereqs=(bam,)      , required_output=(sortedbam,)) # sort by NAME and not position! That's just a thing that HTSeq needs/prefers if it's to handle paired-end files properly
    script_obj.appendCommandUnlessFilesExist(check_files=(outfile,)  , cmd=ht_cmd                   , required_prereqs=(sortedbam,), required_output=(outfile,))
    return

def generateBcpPeakCmd(exper, sampName, script_obj):
    '''BCP is a good all-purpose peak caller. It works especially well for histone-related (broad) peaks. We do not use it for the narrow TF-related peaks, where GEM is preferred.'''
    species            = exper.getSpeciesForSample(sampName)
    expBam             = exper.getAlignedBamForSample(sampName)
    ctrlBam            = exper.getAlignedChIPControlBamForSample(sampName)
    outdir             = exper.getBcpPeakDir(sampName)
    bed_final_out      = exper.getBcpPeakCallBedFile(sampName)
    bcp_fragment_size  = 200     # "The fragment size to which we extend the reads in pre-processing data step. Default:200bp."
    bcp_window_size    = 200     # "The window size we apply the adjcaent window in pre-processing data step. Default:200bp."
    bcp_p_cutoff       =   0.001 # The p_value you want to use for remove false positive based on control data. Range:1e-2---1e-6, default is 1e-3.

    #     bed_out        <- paste0(outPrefix, ".bed") # The results data set with 5 columns you want to output.
    
    # BCP_HM docs are accessible with this unusual command:   BCP_HM -h 1
    #     -1      The ChIP-seq data set you want to input.
    #     -2      The control/input data set you want to input.
    #     -f      The fragment size to which we extend the reads in pre-processing data step. Default:200bp.
    #     -w      The window size we apply the adjacent window in pre-processing data step. Default:200bp.
    #     -p      The p_value you want to use for remove false positive based on control data. Range:1e-2---1e-6, default is 1e-3.
    #     -3      The results data set with 5 columns you want to output.
    #     -h      Use "-h 1" to display the help (not just -h)
    exp_as_bedfile  = expBam  + ".converted_to.bed"
    ctrl_as_bedfile = ctrlBam + ".converted_to.bed"

    # BCP requires BED inputs, not BAM inputs. So we need to convert them using the common tool 'bam2bed'
    exp_to_bed_cmd  = ''.join([BAM2BED_PATH, " --do-not-sort ", " < " + expBam , " | cut -f 1-3  > ", exp_as_bedfile])
    ctrl_to_bed_cmd = ''.join([BAM2BED_PATH, " --do-not-sort ", " < " + ctrlBam, " | cut -f 1-3  > ", ctrl_as_bedfile])
    # Note that we cut the files down to just columns 1 through 3! We do not need any of the extra columns

    bcp_cmd  = BCP_HM_PATH
    bcp_cmd += " -1 " + str(exp_as_bedfile)   # BED file
    bcp_cmd += " -2 " + str(ctrl_as_bedfile)  # BED file
    bcp_cmd += " -f " + str(bcp_fragment_size)
    bcp_cmd += " -w " + str(bcp_window_size)
    bcp_cmd += " -p " + str(bcp_p_cutoff)
    bcp_cmd += " -3 " + str(bed_final_out)

     # There are two versions of BCP - BCP_TF and BCP_HM. The latter is for histone marks - what you are interested in. BCP takes bed files as input.
     # Note: BCP_HM has a very confusing help method: you specify "BCP_HM -h 1" for help. The "1" is mandatory.
     #This is an example command I ran on a sample of H3K4me3 data:
     #BCP_HM -1 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1.bed -2 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed -f 200 -w 200 -p 0.001 -3 GM12878_H3K4me3_Rep1_EncodeAlign_results001_HM.bed
     #Note, there are changes to this command if you have Chip-exo data. You will need to replace  "Read_Distribution_default.txt" by "Read_Distribution_ChIP-exo.txt" and will need to add --smooth 3 option to the commas.
     # -Reuben
     #BCP_HM -1 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1.bed
     #       -2 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed
     #       -f 200 -w 200 -p 0.001 -3 GM12878_H3K4me3_Rep1_EncodeAlign_results001_HM.bed


    script_obj.appendMkDir(outdir) # this directory must exist before we try moving any files to it...
    script_obj.appendCheckForRequiredFiles(file_tup=(expBam, ctrlBam))

    # convert the EXPERIMENT bam file to a bed file
    script_obj.appendCommandUnlessFilesExist(check_files=(exp_as_bedfile,) , cmd=exp_to_bed_cmd , required_output=(exp_as_bedfile,))

    # convert the CONTROL bam file to a bed file
    script_obj.appendCommandUnlessFilesExist(check_files=(ctrl_as_bedfile,), cmd=ctrl_to_bed_cmd, required_output=(ctrl_as_bedfile,))

    # Finally, run the experiment on the input files
    script_obj.appendCommandUnlessFilesExist(check_files=(bed_final_out,)  , cmd=bcp_cmd        , required_output=(bed_final_out,))

    #     outPrefix <- gsub("[.](bam|sam|bed)$", "", outPrefix, ignore.case=T, perl=T) # remove any '.bam'/sam/bed from the end. Those are not correct suffixes. We will add the suffix ourselves
    #     
    #     # There are two versions of BCP - BCP_TF and BCP_HM. The latter is for histone marks - what you are interested in. BCP takes bed files as input.
    #     # Note: BCP_HM has a very confusing help method: you specify "BCP_HM -h 1" for help. The "1" is mandatory.
    #     #This is an example command I ran on a sample of H3K4me3 data:
    #     #BCP_HM -1 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1.bed -2 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed -f 200 -w 200 -p 0.001 -3 GM12878_H3K4me3_Rep1_EncodeAlign_results001_HM.bed
    #     #Note, there are changes to this command if you have Chip-exo data. You will need to replace  "Read_Distribution_default.txt" by "Read_Distribution_ChIP-exo.txt" and will need to add --smooth 3 option to the commas.
    #     # -Reuben
    #     #BCP_HM -1 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1.bed
    #     #       -2 /home/rthomas/ChipSeqChallenge/RealData/HistoneData/Monkey/Results/04a.mapping/wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed
    #     #       -f 200 -w 200 -p 0.001 -3 GM12878_H3K4me3_Rep1_EncodeAlign_results001_HM.bed
    #     
    #     # These are intermediate files. The problem is that bcp cannot deal with BAM files directly, so we have to convert them.
    #     bed_exp <- "/04a.mapping/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1.bed" # "The ChIP-seq data set you want to input."
    #     bed_ctl <- "/wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bed" # "The control/input data set you want to input."
    #     
    #     bam2bed_exp_cmd <- agw_construct_cmd(get_exe_path("bam2bed"), " --do-not-sort ", " < ", ctlBam, " > ", bed_exp)
    #     bam2bed_ctl_cmd <- agw_construct_cmd(get_exe_path("bam2bed"), " --do-not-sort ", " < ", expBam, " > ", bed_ctl)
    return  # end of def generateBcpPeakCmd(...)

def generateGemPeakCmd(exper, sampName, script_obj):
    '''Gem is best for TF-related peaks (ones with sequence motifs)'''
    species         = exper.getSpeciesForSample(sampName)
    outdir          = exper.getGemPeakDir(sampName)
    expBam          = exper.getAlignedBamForSample(sampName)
    ctrlBam         = exper.getAlignedChIPControlBamForSample(sampName)
    gemReadDistFile = GEM_READ_DIST_FILE
    
    GEM_THREADS      = 1 # gem uses a ton of threads for some reason, so 1 is even conservative
    GEM_QVAL         = 2	      # sets the q-value threshold (-log10, so 2 = 0.01)
    GEM_KMIN         = 6	      # minimum kmer length. From our default UCSF settings used in the 'Monkey' pipeline
    GEM_KMAX         = 13	      # maximum kmer length. From our default UCSF settings used in the 'Monkey' pipeline
    genomeByChromFastaDir = getAnnot(species, "fa_by_chrom")
    chrSizesFile          = getAnnot(species, "chr_sizes")

    cmd  = JAVA_PATH + " -Xmx" + str(GEM_RAM_GB) + "G" + " -jar " + GEM_JAR_PATH
    cmd +=      " --t " + str(GEM_THREADS)
    cmd +=      " --q " + str(GEM_QVAL)
    cmd +=  " --k_min " + str(GEM_KMIN) +  " --k_max " + str(GEM_KMAX)
    cmd +=      " --d " + gemReadDistFile  # readdist file / read distribution file
    cmd +=      " --g " + chrSizesFile
    cmd += " --genome " + genomeByChromFastaDir
    cmd +=   " --expt " + expBam
    cmd +=  (" --ctrl " + ctrlBam) if ctrlBam is not None else "" # can be omitted if there is NO input file
    cmd += " --f SAM " # <-- "SAM" is required in order to accept SAM or BAM files, otherwise BED files are the default
    cmd += " --sl --outBED "
    cmd += " --out " + outdir # outprefix creates a new folder with this name, respecting existing subdirectories.

    expectedOutFile = exper.getGemOutHTMLFile(sampName)
    script_obj.appendMkDir(outdir) # this directory must exist before we try moving any files to it...
    script_obj.appendCheckForRequiredFiles(file_tup=(expBam, ctrlBam, chrSizesFile, genomeByChromFastaDir, gemReadDistFile))
    script_obj.appendCommandUnlessFilesExist(check_files=(expectedOutFile,), cmd=cmd, required_output=(expectedOutFile,))
    return # end of def generateGemPeakCmd(...)

def generateEdgeRCmd(count_files, groups, outdir, outfile, script_obj):
    cmd  = RSCRIPT_PATH + " " + EDGER_SCRIPT
    cmd += " --verbose "
    cmd += " --groups=" + enquote(DELIM_FOR_ARGS_TO_EDGER_SCRIPT.join(groups))
    cmd += " --countfiles=" + enquote(DELIM_FOR_ARGS_TO_EDGER_SCRIPT.join(count_files))
    cmd += " --outfile=" + enquote(outfile)
    script_obj.appendMkDir(outdir) # this directory must exist before we try moving any files to it...
    script_obj.appendCheckForRequiredFiles(file_tup=(count_files,))
    script_obj.appendCommandUnlessFilesExist(check_files=(outfile,), cmd=cmd, required_prereqs=(), required_output=(outfile,))
    return


def generateAlignmentCmd(fastq1, fastq2, species, outdir, outbam, aligner, script_obj):
    '''Note that this is will only generate ONE alignment command per file, even if it's called multiple times.'''
    xAssert(outbam is not None, "Somehow the output bam is NONE!")
    fqkey = " | ".join([str(fastq1), str(fastq2), str(species), str(aligner)]) # just a unique key for this pair of reads, species, and aligner
    if fqkey in global_file_paths_already_aligned:
        verbosePrint("[Skipping...]: Already generated an alignment command for the file/pair " + str(fastq1) + " and " + str(fastq2) + " with species " + str(species) + " and the aligner " + str(aligner))
        return # exit early
    
    global_file_paths_already_aligned[fqkey] = True # remember that we saw this file...
    isPaired            = fastq2 is not None
    innerDistStr        = " --mate-inner-dist=150 " if isPaired else " "
    bowtieIndex         = getAnnot(species, "bowtie_index")
    all_reads_bam       = os.path.join(outdir, "all.bam") # an intermediate file
    sortbam             = outbam # <-- final output is sorted and aligned reads

    if aligner == TOPHAT_PATH:
        gtf  = getAnnot(species, "gtf")    # <-- only actually used by tophat
        cmd  = TOPHAT_PATH + " -o " + outdir 
        cmd += " " + "--min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search "
        cmd += " " + "--GTF=" + gtf
        cmd += " " + "--num-threads=" + str(TOPHAT_N_THREADS)
        cmd += " " + innerDistStr
        cmd += " " + bowtieIndex
        cmd += " " + fastq1 + " "
        cmd += " " + fastq2  if isPaired  else "" # might or might not have a second file in the pair!
        TOPHAT_HARDCODED_BAM_OUT_NAME = "accepted_hits.bam" # tophat always generates files with this output name, no matter what.
        tophat_outbam = os.path.join(outdir, TOPHAT_HARDCODED_BAM_OUT_NAME)
        cmd += "\n" + "test -e " + sortbam + " || mv -f " + tophat_outbam + " " + sortbam # move it to the FINAL bam location. Which might be the same! "test" is so that we ONLY MOVE IT if the target does not already exist. Otherwise do not move it. So if the filenames are the same, then don't move it, since that generates an error code.
    elif aligner == BOWTIE_PATH:
        inputFqString = " ".join(["-1", enquote(fastq1), "-2", enquote(fastq2)])   if isPaired else   " ".join(["-U", enquote(fastq1)])
        cmd  = BOWTIE_PATH
        cmd += " " + "--threads=" + str(BOWTIE_N_THREADS)
        cmd += " " + "-x " + bowtieIndex
        cmd += " " + inputFqString + " "
        cmd += " | "  + SAMTOOLS_PATH + " view -@ " + str(SAMTOOLS_N_THREADS) + " -b - " +                        " > " + all_reads_bam
        cmd += " && " + SAMTOOLS_PATH + " view -@ " + str(SAMTOOLS_N_THREADS) + " -b -F 0x104 " + all_reads_bam + " "           # samtools view -b -F 0x104 in.bam > mapped_primary.bam This is how you get ONLY primary mapped reads
        cmd += " | samtools sort - " + " > " + sortbam
        cmd += " && " + " /bin/rm " + all_reads_bam # delete the intermediate (unsorted) "all_reads" bam file
    else:
        raise Exception("Unrecognized aligner! We currently only know about tophat and bowtie2. Note that this is case-sensitive. This is a CODING bug and should not be related to a user-specified file. Your specified aligner was: " + aligner)
        pass

    script_obj.appendMkDir(outdir)
    script_obj.appendCommandUnlessFilesExist(check_files=(sortbam,), cmd=cmd, required_prereqs=(fastq1, fastq2), required_output=(sortbam,))
    return

# Must come at the VERY END!
if __name__ == "__main__":
    handleCmdLineArgs()
    pass

