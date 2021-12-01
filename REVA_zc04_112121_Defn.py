#!/usr/bin/env python -i

## Version REVA_zc04_112121
PreREVA1 = False

## Modules to import
import array
from time import time, strftime, localtime, asctime
try:
    import cPickle
except:
    import pickle as cPickle
import numpy as np
from itertools import chain
try:
    from itertools import imap, izip, izip_longest 
except:
    imap=map
    izip=zip
    from itertools import zip_longest as izip_longest
try:
    temprange1 = xrange(10)
except:
    xrange = range
import gzip
import zipfile
from operator import ne
from sys import version,argv,platform,getsizeof
import os
import subprocess
import re
from random import choice
from glob import glob ## fixed a missing import statement [fix 07-26-21, thanks to Lamia]
def ham1(a,b):
    return(sum(imap(ne,a,b))+abs(len(a)-len(b)))  ## quick hamming distance (with penalties for any mismatch in length

## Constants relevant to the execution of the program
klen1=25  ## the Kmer length used in searches.  Recommended that this be an odd number to avoid any individual k-mer being its own antisense (which can lead to ambiguities in counting) 
BitsIndexed1=25 ## REVA uses a combination of indexing and sorting to keep track of k-mers.  BitsIndexed determines how many bits are accounted for by each.  Default is 25.
RefSeqFile1='No_File_Specified_For_Reference_Please_Specify_Your_FastA_Reference_Sequence_File_In_Command_Line' ##'miminal_hg38_Nsreplaced.fa' #### /Users/firelab08/Desktop/BigData/CircleFinder/ws220.fa' ##'OP50Mock.fa' ## Reference sequence in FastA format
IndexFile1='default' ###'/Users/firelab08/Desktop/BigData/CircleFinder/ws220_kIndex_Kmer25_IndexedBits25_FirstOccurencePlusMultiplicities.hqa2' ## A premade index file that will speed things up
UCSCAssembly1 = 'default'  ## this is, in particular the species identifier for the UCSC browser.  Current values are ce11 for worm and hg# for human
UCSCLinks1 = True ## Setting this to False skips UCSC link output
CircleMax1=100000 ##Maximum Circle Length that we'll look for (everything larger is assumed to be a structural rearrangement)
DeletionMax1=1000000 ##Maximum deletion Length that we'll look for (everything larger is assumed to be a structural rearrangement)
InsertionMin1=4 ## Minimum Circle Length: operationally this is a sequencing error filter.  Although it gets set in principal to the smallest circle.  The smallest simple deletion that could come from sequencing error is fine(4 is pretty much okay) 
Tn5DupMax1=16 ## Maximum tagmentation overlap for detecting singly-tagmented circles (16 should be fine) 
Tn5DupMin1=0  ## Minimum tagmentation overlap for detecting singly-tagmented circles (0 should be fine) 
SeparationGranularity1=10000 ## Granularity of storage for ciclular events (histogram or array)
ReadSeparationMax1=3000 ## Largest fragment we expect to be able to capture and sequence on the flow cell
SeparationMin1=0 ## Ignore all read pairs where and end fails to map, map to different chromosomes, or where fragment length is < this value
SeparationMax1=0 ## Ignore all read pairs where and end fails to map, map to different chromosomes, or where fragment length is > this value
AvoidR1Ambiguous1=False ## True will ignore any read pair without a unique-mapping k-mer in R1
AvoidR2Ambiguous1=False ## True will ignore any read pair without a unique-mapping k-mer in R2
AvoidBothAmbiguous1=False ## True will ignore any read pair where neither R1 or R2 has a unique matching k-mer
AvoidEitherAmbiguous1=False ## True will ignore any read pair where either R1 or R2 failes to have a unique matching k-mer
ReportInterval1=20000 ## Lines between reporting intervals (text messages that indicate progress... not relevant to output of the progrma)
LinesToProcess1= -1 ## 100000 ## number of lines to process in the input files (-1 for everything)
VerboseOutput1=True  #Several extra parameters reported for potential rearranged reads and read pairs
DomainLength1=3000000 ##Length limit for all instances of an individual k-mer.  All instances of the k-mer must be within this limit and on the same chromosome for the repeat to be considered "focal".  Otherwise, the repeat is considered "dispersed".
Short1,Long1 = 100,200 ##Cutoffs in total span for a read pair to be considered short, medium, or long.  Default values are 100 for Short1, 100 for Long1.
MinRepeatEvalue1=float('1.0E-6') ## Minimum e-value for repeat elements (from a DFAM database file) to be considered as identified
CircularChromosomeLimit1=13500000 ## Above this size all chromosomes are considered linear
OtherSVCandidates1 = False ## Report translocation and inversion candidates (junctions between chromosomes or between sense and antisense strands respectively)
MPLS1 = False  ## if MPLS1 is False, use high stringency matching (two nonoverlapping unique k-mers miust match to call junction
test1='TestSVsUnc13R1.fastq.gz'
test2='TestSVsUnc13R2.fastq.gz'
Feature1D={}
FeatureTags1=[]
Feature1 = False
FiL0=[] ##'/Users/firelab08/Desktop/BigData/ReadEvaluation/VC2010SingleWorms/VC2010-1-PCR4_S4_L001_R1_001.fastq.gz']## aardvark ##[]  ## A list of fastq Read files to match
FiD0=[] ##'R1' or 'R2' Indicates that this is read1 or read2
FiN0=[] ## Potential 'Narratives' for each read file.  Default is the file name but can include additional information after a '#' indicator
GFF1 = [] ## List of GTF/GFF data files for coverage calculation
t0=time()
MyCode='default'
Mnemonic1=[]
chrbase1='default'
DFAM1=''
DFAMCount1=False
ai1=0
FullCoverage1=True
BriefCoverage1=True
CoverageAccumulationInterval1=200000 ## Number of read pairs between updates of cumulative coverage zero for no coverage updates
FullCoverageInterval1=0 ## Zero instructs REVA to only calculate coverage once at the end of the cycle, otherwise the program will output aggretgate coverage after each CoverageInterval1 reads
FullCoverageByFile1=False ## True outputs the cumulative coverage after each examined file
KmerBinCoverageUnique1 = True
KmerBinCoverageAll1 = True

KeepDoubleExtension1 = False  ## when looking for 5' and 3' extensions, setting this to false throws away any sequences extended at both ends.  This is to avoid contaminating sequences that might match in one or a few central k-mers
SnpFilter1 = 2 ## this will filter out sequences that have an extended match followed by mismatch followed by match as likely not having extensions.  Setting this to 2 will filter out any apparent extension (5' or 3') which has 2 or more matches following the first mismatch
R1Buffer5 = 0 ## automatically trim this number of bases from the 5' end of each R1 read (after barcode elimination)
R1Buffer3 = 0 ## automatically trim this number of bases from the 3' end of each R1 read (after linker elimination)
R2Buffer5 = 0 ## automatically trim this number of bases from the 5' end of each R2 read (after barcode elimination)
R2Buffer3 = 0 ## automatically trim this number of bases from the 3' end of each R2 read (after linker elimination)
ReportBuffers1 = False ## Are the buffer sequences reported to the Junctions file (useful for Unique-Molecular-Identifier applications)

## Reporting of relative read coverage/starts/ends for regions and features 
## What are we counting
ReportReads1 = True ## Report the number of reads with start/end/k-mer in a given region
ReportPositions1 = True ## Report the number of positions in the sequence with start/end/kmer mapping to that position

ReportStarts1 = True ## Report "Starts" (literally extrapolated first match positions identified by first mapping kmer in each read) 
ReportEnds1 = True ## Report "Ends" (literally extrapolated last match positions identified by first mapping kmer in each read)
ReportCoverage1 = True ## Report coverage by every observed k-mer in each read  (not relevant for PreReva)  

##Seq Index Tools (including plurality-- which finds the most common index combination and looks at only reads with that index combination)
SeqIndexMode1 = 0  # 0 for no filtering on sequence index, 1 for plurality [most common], 2 for minority [least common], eventually holds the sequences of one or both required indices as a tuple
SeqIndex1 = ''  ## Narrative list of seq indices to either include (SeqIndexMode=1) or exclude (SeqIndexMode=2)
SeqIndexD1 = {} ## Dictionary of sequence indices to include or exclude (depending on SeqIndexMode1)
def GetSeqIndexFromReadID1(id1):  ## Input is description line in FastQ file, output is the index sequence (upper case)
    return id1.strip().split(':')[-1].upper()  ## assumes index is the bases following the terminal colon in the 

## How to handle two strands
ReportSense1 = True ## Report sense reads separately
ReportAntisense1 = True ## Report antisense reads separately
ReportBothStrands1 = False ## Report sum of sense and antisense reads

## How to handle repeated sequences
ReportUnique1 = True ## Report unique regions
ReportRepeats1 = True ## Report Local, Chromosomal, and Dispersed repeated regions (not relevant for PreReva)

FindStructuralAnomalies1=True ## Output a file of possible structural anomalies either from paired end or split read inference
CombinationList1=0 ## Set CombinationList1 to 1 to make a raw list of read pair start combinations and counts for all positions,
                   ## set CombinationList1 to 2 for just sequences in Base-ByBase output
                   ## set CombinationList1 to 3 for just paired end reads that meet criteria for a Tn5 tagmentation (8-10 base overlap with circular topology) 
RequireKMerOnly1=True ## Setting this to true outputs potential structural variants based only on k-mer matches not imposing a secondary requirement that all other bases on either side of the potential variant match.  Setting this to False is mostly of use with highly accurate reference genomes, accurate sequencing, and relatively short read lengths
MatePairSubstitutionMax1=1 ## When aligning sequences and calling variants using mate pairs how much variation (substitution tolerance) is allowed between read and genome
SplitReadSubstitutionMax1=1 ## when aligning sequences at the ends of possible variants, how many substitutions are allowed
MatePairIndelCountMax1=1 ## When aligning sequences and calling variants using mate pairs how manuy indels are allowed between read and genome
MatePairIndelLengthMax1=3 ## When aligning sequences and calling variants using mate pairs what is the maximum length of indels is allowed between read and genome
SplitReadIndelCountMax1=1 ## when aligning sequences at the ends of possible variants, how many indels are allowed (maximum value, inclusive of the number given)
SplitReadIndelLengthMax1=3 ## when aligning sequences at the ends of possible variants, how many bases of indels are allowed
CoverageDType1=np.uint16
CoverageSEMax1=255   ## maximum coverage that will be recorded in looking for starts and ends in basebybase output.  Setting to 255 conserves memory, setting to higher number will eat up a bit of memory but allow accurate numbers in high coverage regions
MaxCoverage1=2**16-1
MaxCoverHist1=100
MaxMultiplicityReported1=63 ## If multiplicity is reported, should we max out at some value (usually 63).  Setting at 63 results in best memmory usage.  For larger numbers set at 2**(2*n)-1 or 16383 to use two bytes per position, 1073741823 for four bytes and larger for 8 bytes
MultiplicityDataType1=np.uint8
IndexConstructionBins1 = 0
Debug1 = False ## For debug purposes only read and operate with the first chromosome/FastA entry from the reference file
TempFileUID1 = choice('QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm')  ## This will be a seven character unique ID for temp files

BarcodeSeqR1 = ''
BarcodeRequireR1 = False
BarcodeSeqR2 = ''
BarcodeRequireR2 = False
BarcodeLenR1 = 0
BarcodeLenR2 = 0
AnyBarcodeR1 = False
AnyBarcodeR2 = False
LinkerR1 = ''
LinkerRequireR1 = False
LinkerR2 = ''
LinkerRequireR2 = False
LinkerLenR1 = 0
LinkerLenR2 = 0

BinByBin1 = True
ChromosomeByChromosome1 = True

BaseByBase1 = False
BaseByBaseOutFile1 = 'default'
BaseByBaseRegionsFiles1 = [] ## will take from infile if available... if a GTF/GFF is specified in the BaseByBaseCall, will use that; if there is no GTF/GFF file in the BaseByBase call but one used in the command line, will use that, otherwise will use all positions.
BaseByBaseCategories1 = {}
BaseByBaseTags1 = []
BaseByBasePositions1 = []
BaseByBaseColumns1 = 'CPFOBRMKSET'
for vL in BaseByBaseColumns1:
    vars()['bbb'+vL+'1'] = True

MaxExtension1 = 4  ## the maximun length of recorded extension for base-by-base analysis
UCSCBuffer1 = 1000 ## Buffer to add to UCSC display around any gene that is shown

FivePrimeExtensionDisallowed1 = False  ## Filtering based on homology match for R1.  This eliminates any sequence where the first base doesn't match
MinHomology1 = 0  ## Allows option to filter based on minimum homology length for R1 default is zero (no filtering)
MaxHomology1 = 999999999 ## Allows option to filter based on maximum homology length for R1 default is 999999999 (essentially no filtering)
StartHomology1 = ''  ## Allows option to filter based on start of homology sequence (e.g. 'G' insists first base of homology is a 'G') default is no filtering

SRRList1 = [] ## A List of SRR IDs for analysis from the short read archive
SRRTempList1 = [] ## A local file name list for SRR (list of files to eventually be deleted)
MetaData1 = True ## True instructs REVA to display all available metadata prominently in each output file.

DeleteNs1 = False ## Seting this to true does a prefiltering of experimental sequence (not reference genome in which all Ns are deleted.  Default is to convert all Ns to Gs.
FirstKTable1 = False ## Report a table with the first position of a unique k-mer in each sequence
FirstKByStrand1 = False ## Report separate firstK values for sense and antisense matches
InterimKTable1 = False ## Do interim K table reports -- note that this will turn off the global firstKTable.
klen1d = 0 ## This will hold the max value for unmatched 5' or 3' end length (default will be 2*klen1)

R1Only = False ## Report only data from the R1 file?
R2Only = False ## Report only data from the R2 file?

FastQDumpProgram1 = 'fastq-dump' ## Location of a program that will download files from SRA if needed, can be reset by setting FastQDump=<path to fastq-dump executable>.
for i in range(8):
    TempFileUID1 += choice('QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm')
ScratchPath1 = os.path.join(os.getcwd(), 'REVATempFiles')
if not(os.path.isdir(ScratchPath1)):
    os.mkdir(ScratchPath1)

def R2FileTest1(n):  ## Look for any occurence of 'r2' in a name with no digit afterwards
    n0 = n.lower()
    if n0.endswith('r2'):
        return True
    if 'r2' in n0:
        n1 = n0.split('r2')
        for n2 in n1[1:]:
            if not(n2[0].isdigit()):
                return True
    return False
## parse command line arguments

CircleMaxSet1 = False
DeletionMaxSet1 = False ## Keep track of whether these have been set and if so don't reset them automatically if OtherSV is set to True

while ai1<len(argv):
    a1=argv[ai1]
    ai1+=1
    a1=a1.replace('"','').replace("'","")
    if a1=='choose' or a1.lower()=='-c' or a1.lower()=='c':
        from Tkinter import Tk
        root=Tk()
        root.attributes("-topmost", True)
        from tkFileDialog import askopenfilenames
        R1File1=askopenfilenames(title='Files with R1 reads',initialdir=os.getcwd(), filetypes=[("gz files","*.gz"),("Fastq files","*.fastq")])
        root.destroy()
        FiL0.extend(R1File1.split('#')[0])
        FiN0.extend(R1File1)
        FiD0.extend('R1'*len(R1File1))
        continue
    if a1.lower()=='exon' or a1.lower()=='-e' or a1.lower()=='e':
        Feature1D['exon'] = len(Feature1D)+1
        continue
    if a1.lower()=='gene' or a1.lower()=='-g' or a1.lower()=='g':
        Feature1D['gene'] = len(Feature1D)+1
        continue
    if not('=' in a1) and a1.lower().endswith('.py'):
        MyCode1=a1.strip()
        continue
    if not('=' in a1) and a1.lower().startswith('deleten'):
        DeleteNs1 = True
        continue
    if not('=' in a1) and a1.lower().startswith('firstkt'):
        FirstKTable1 = True
        continue
    if not('=' in a1) and a1.lower().startswith('interimkt'):
        InterimKTable1 = True
        continue
    if not('=' in a1) and (a1.lower().startswith('firstkstr') or a1.lower().startswith('firstkbystr')):
        FirstKByStrand1 = True
        continue
    if a1.lower().startswith('avoidr1') and not('false'  in a1.lower()):
        AvoidR1Ambiguous1 = True
        continue
    if a1.lower().startswith('avoidr2') and not('false'  in a1.lower()):
        AvoidR2Ambiguous1 = True
        continue
    if a1.lower().startswith('avoideither') and not('false'  in a1.lower()):
        AvoidR1Ambiguous1 = True
        AvoidR2Ambiguous1 = True
        continue
    if a1.lower().startswith('avoidboth') and not('false'  in a1.lower()):
        AvoidBothAmbiguous1 = True
        continue
    if a1.lower().startswith('r1only') and not('false'  in a1.lower()):
        R1Only = True
        continue
    if a1.lower().startswith('r2only') and not('false'  in a1.lower()):
        R2Only = True
        continue
    if not('=' in a1) and (a1.lower().endswith('.fa') or a1.lower().endswith('.fasta') or a1.lower().endswith('.fa.gz') or a1.lower().endswith('.fasta.gz')  or a1.lower().endswith('.fa.zip') or a1.lower().endswith('.fasta.zip')):
        RefSeqFile1=a1.strip()
        continue
    if not('=' in a1) and (a1.lower().endswith('.gtf') or a1.lower().endswith('.gtf.gz') or a1.lower().endswith('.gff') or a1.lower().endswith('.gff.gz')  or a1.lower().endswith('.gff3') or a1.lower().endswith('.gff3.gz')):
        GFF1.append(a1.strip())
        continue
    if not('=' in a1) and (a1.lower().endswith('files') or a1.lower().endswith('files.txt')):
        R1File1=a1
        Mnemonic1.append(a1)
        for fx1 in R1File1.read().replace(',',' ').split():
            if os.path.isfile(fx1):
                FiL0.append(fx1.split('#')[0])
                FiN0.append(fx1)
                iD0 = 'R1'
                if R2FileTest1(fx1):
                    if os.path.isfile(fx1.replace('R2','R1')) or os.path.isfile(fx1.replace('r2','r1')):
                        iD0 = 'R2'
                    if os.path.isfile(fx1.replace('R2','R1',1)) or os.path.isfile(fx1.replace('r2','r1',1)):
                        iD0 = 'R2'
                    if os.path.isfile(fx1[::-1].replace('2R','1R',1)[::-1]) or os.path.isfile(fx1[::-1].replace('2r','1r',1)[::-1]):
                        iD0 = 'R2'
                FiD0.append(iD0)
        continue
    if not('=' in a1) and (a1.lower().startswith('srr')):
        for a222 in a1.split(','):
            SRRList1.append(a222)
        continue
    if  not('=' in a1) and a1.lower().startswith('reportread'):
        ReportReads1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportposition'):
        ReportPositions1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportstart'):
        ReportStarts1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportend'):
        ReportEnds1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportcover'):
        ReportCoverage1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportsens'):
        ReportSense1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportantisens'):
        ReportAntisense1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportboth'):
        ReportBothStrands1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportuniqu'):
        ReportUnique1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportrepeat'):
        ReportRepeats1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('basebybase'):
        BaseByBase1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('binbybin'):
        BinByBin1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('chromosomebychromosome'):
        ChromosomeByChromosome1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('ucsc'):
        UCSCLinks1 = True
        continue
    if not('=' in a1) and (a1.lower().endswith('fastq.gz') or a1.lower().endswith('.fastq')):
        if 'R1' in a1 and not('R2' in a1):
            R1File1=a1.strip().replace(',',' ').split()
            for R1File1a1 in R1File1:
                FiL0.append(R1File1a1.split('#')[0])
                FiN0.append(R1File1a1)
                FiD0.append('R1')
            continue
        if 'R2' in a1 and not('R1' in a1) and not(('SRR2' in a1) and (a1.count('R2')==1)):
            R2File1=a1.strip().replace(',',' ').split()
            for R1File1a1 in R1File1:
                FiL0.append(R2File1a1.split('#')[0])
                FiN0.append(R2File1a1)
                FiD0.append('R2')
            continue
        else:
            R1File1=a1.strip().replace(',',' ').split()
            for R1File1a1 in R1File1:
                FiL0.append(R1File1a1.split('#')[0])
                FiN0.append(R1File1a1)
                FiD0.append('')
            continue
    if not('=' in a1) and (a1.lower().endswith('dfam.hits.gz') or a1.lower().endswith('dfam.hits') or a1.lower()=='dfam' or a1.lower()=='-dfam'):
        if not('.hits' in a1):
            if 'ws220' in RefSeqFile1.lower():
                DFAM1='ce10_dfam.hits.gz'
            elif 'hg38' in RefSeqFile1.lower():
                DFAM1='hg38_dfam.hits.gz'
        else:
            DFAM1=a1
        DFAMCount1=True
        continue
    if not('=' in a1) and a1.lower().endswith('wv2.9re'):
        IndexFile1=a1.strip()
        continue
    if a1.lower()=='-1' or a1.isdigit():
        LinesToProcess1=int(a1.strip())
        continue
    if a1.lower().startswith('nocover'):
        FullCoverage1=False
        FullCoverageInterval1=0 
        FullCoverageByFile1=False 
        continue
    if a1[0]=='-':
        a11=a1.strip()[1:].lower()
        a22=argv[ai1].strip().replace('"','').replace("'","")
        ai1+=1
    else:
        a11=a1.split('=')[0].strip().lower()
        a22=a1.split('=')[-1].strip()
    if a11.startswith('dfam'):
        DFAM1=a22
        DFAMCount1=True
    elif a11.startswith('gfffile') or a11.startswith('featurefile') or a11.startswith('gtffile') or a11.startswith('gffdatafile') or a11.startswith('gtfdatafile'):
        for a222 in a22.split(','):
            GFF1.append(a222)
    elif a11.startswith('srr') or a11.startswith('sra'):
        for a222 in a22.split(','):
            SRRList1.append(a222)
        continue
    elif a11.startswith('featuretag'):
        for a222 in a22.split(','):
            FeatureTags1.extend(a222.lower().split(','))
    elif a11.startswith('featurecategory') or a11.startswith('gfffeature') or a11.startswith('gtffeature') or a11.startswith('feature'):
        for a222 in a22.split(','):
            Feature1D[a222] = len(Feature1D)+1
    elif a11.startswith('r1buffer5'):
        R1Buffer5=int(a22)
    elif a11.startswith('r1buffer3'):
        R1Buffer3=int(a22)
    elif a11.startswith('r2buffer5'):
        R2Buffer5=int(a22)
    elif a11.startswith('r2buffer3'):
        R2Buffer3=int(a22)
    elif a11.startswith('snpfilter'):
        SnpFilter1=int(a22)
    elif a11.startswith('fastqdump') or a11.startswith('fastq-dump'):
        if os.path.isfile(a22) and not(os.path.isdir(a22)):
            FastQDumpProgram1 = a22
        elif os.path.isdir(a22):
            FastQDumpProgram1 = os.path.join(a22,'fastq-dump')
            if not(os.path.isfile(FastQDumpProgram1)):
                FastQDumpProgram1 = os.path.join(a22,'bin','fastq-dump')
                if not(os.path.isfile(FastQDumpProgram1)):
                    os.environ["PATH"]=os.getenv("PATH")+':'+a22
    elif a11.startswith('r1'):        
        for R1File1 in a22.strip().replace(',',' ').split():
            if '*' in R1File1:
                Mnemonic1.append(a22.replace('*','_all_'))
                for fx1 in glob(R1File1):
                    FiL0.append(fx1)
                    FiN0.append(fx1)
                    FiD0.append('R1')
                    if R2FileTest1(fx1):
                        if os.path.isfile(fx1.replace('R1','R2')):
                            Mnemonic1[-1] = Mnemonic1[-1].replace(('R1',''))
                        if os.path.isfile(fx1.replace('R1','R2'),1):
                            Mnemonic1[-1] = Mnemonic1[-1].replace(('R1',''),1)
                        if os.path.isfile(fx1[::-1].replace('1R','2R')[::-1],1):
                            Mnemonic1[-1] = Mnemonic1[-1][::-1].replace('1R','')[::-1]
            elif R1File1.lower().endswith('files') or R1File1.lower().endswith('files.txt'):
                Mnemonic1.append(a22)
                for fx1 in R1File1.read().replace(',',' ').split():
                    if os.path.isfile(fx1):
                        FiL0.append(fx1.split('#')[0])
                        FiN0.append(fx1)
                        FiD0.append('R1')
            else:           
                FiL0.append(R1File1.split('#')[0])
                FiN0.append(R1File1)
                FiD0.append('R1')
    elif a11.startswith('r2'):
        for R2File1 in a22.strip().replace(',',' ').split():
            if '*' in R2File1:
                if not(a22.replace('*','_all_').replace('R1','R2') in Mnemonic1):
                    Mnemonic1.append(a22.replace('*','_all_'))
                for fx1 in glob(R2File1):
                    FiL0.append(fx1)
                    FiN0.append(fx1)
                    FiD0.append('R2')
            elif R2File1.lower().endswith('files') or R2File1.lower().endswith('files.txt'):
                Mnemonic1.append(a22)
                for fx1 in R2File1.read().replace(',',' ').split():
                    if os.path.isfile(fx1):
                        FiL0.append(fx1.split('#')[0])
                        FiN0.append(fx1)
                        FiD0.append('R1')
            else:           
                FiL0.append(R2File1.split('#')[0])
                FiN0.append(R2File1)
                FiD0.append('R2')
    elif a11.startswith('ref'):
        RefSeqFile1=a22
    elif a11.startswith('indexfile'):
        IndexFile1=a22
    elif a11.startswith('seqindex'):
        if a22.lower().startswith('p') or a22.lower().startswith('1') :
            SeqIndexMode1 = 1   ## plurality -- keeps only read pairs with the most prevalent combination of 5' and 3' indices
        elif a22.lower().startswith('m') or  a22.lower().startswith('2'):
            SeqIndexMode1 = 2   ## minority -- keeps only read pairs with something other than the most prevalent combination of 5' and 3' indices
        else:
            for SeqIndex1s in a22.split('),('):
                SeqIndex01 = SeqIndex1s.split(',')[0].strip().replace('"','').replace("'",'').replace(")",'').replace("(",'')
                if len(a22.split(','))>1:
                    SeqIndex02 = SeqIndex1s.split(',')[1].strip().replace('"','').replace("'",'').replace(")",'').replace("(",'')
                else:
                    SeqIndex02 = ''
                SeqIndexD1[(SeqIndex01.upper(),SeqIndex02.upper())] = 0        
    elif a11.startswith('bit'):
        BitsIndexed1=int(a22)
    elif a11.startswith('minre'):
        MinRepeatEvalue1=int(a22)
    elif a11.startswith('minhomol'):
        MinHomology1=int(a22)
    elif a11.startswith('maxhomol'):
        MaxHomology1=int(a22)
    elif a11.startswith('coveragesemax'):
        if int(a22)==0 or int(a22)==-1:
            CoverageSEMax1=2**32-1
        elif int(a22)<=64:
            CoverageSEMax1=2**int(a22)-1        
        else:
            CoverageSEMax1=int(a22)
    elif a11.startswith('circlemax'):
        CircleMax1=int(a22)
        CircleMaxSet1 = True
    elif a11.startswith('deletionmax'):
        DeletionMax1=int(a22)
        DeletionMaxSet1 = True
    elif a11.startswith('circlemin'):
        InsertionMin1=int(a22)
    elif a11.startswith('insertionmin'):
        InsertionMin1=int(a22)
    elif a11.startswith('tn5dupmax'):
        Tn5DupMax1=int(a22)
    elif a11.startswith('tn5dupmin'):
        Tn5DupMin1=int(a22)
    elif a11.startswith('ucscbuffer'):
        UCSCBuffer1=int(a22)
    elif a11.startswith('reportinterval'):
        ReportInterval1=int(a22)
    elif a11.startswith('separationmin'):
        SeparationMin1=int(a22)
    elif a11.startswith('separationmax'):
        SeparationMax1=int(a22)
    elif a11.startswith('separation') or a11.startswith('gran'):
        SeparationGranularity1=int(a22)
    elif a11.startswith('readsep'):
        ReadSeparationMax1=int(a22)
    elif a11.startswith('domain'):
        DomainLength1=int(a22)
    elif a11.startswith('short'):
        Short1=int(a22)
    elif a11.startswith('long'):
        Long1=int(a22)
    elif a11.startswith('chrbase'):
        chrbase1=a22
    elif a11.startswith('starthomology'):
        StartHomology1 = a22.strip('"').strip("'").upper()
    elif a11.startswith('basebybasecolumns'):
        BaseByBase1 = True
        bbOptionList1 = ['C','P','F','O','B','R','M','K','S','E','T']
        for bbOs in bbOptionStrings1:
            if not('.' in bbOs):
                for vL in bbOptionList1:
                    vars()['bbb'+vL+'1'] = False
                for vS in bbOptionString1:
                    vars()['bbb'+vS+'1'] = True
    elif a11.startswith('basebybaseregions') or a11.startswith('basebybasefile') or a11.startswith('basebybasefea'):
        BaseByBase1 = True
        BaseByBaseRegionsFiles1.extend(a22.split(','))
    elif a11.startswith('basebybasecategor'):
        BaseByBase1 = True
        for a222 in a22.split(','):
            BaseByBaseCategories1[a222] = len(BaseByBaseCategories1)+1
    elif a11.startswith('basebybasetag'):
        BaseByBase1 = True
        for a222 in a22.split(','):
            BaseByBaseTags1.append(a222.strip('"').strip("'"))
    elif a11.startswith('basebybaseposition'):
        BaseByBase1 = True
        BaseByBasePositions1 =  a22.split(',')                     
    elif a11.startswith('basebybase') :
        if a22.lower().startswith('f'):
            BaseByBase1=False
        else:
            BaseByBase1=True
    elif a11.startswith('binbybin') :
        if a22.lower().startswith('f'):
            BinByBin1=False
        else:
            BinByBin1=True
    elif a11.startswith('othersv') :
        if a22.lower().startswith('f'):
            OtherSVCandidates1=False
        else:
            OtherSVCandidates1=True
            if not(CircleMaxSet1):
                CircleMax1 = 2**63
            if not(DeletionMaxSet1):
                DeletionMax1 = 2**63
    elif a11.startswith('mpls') or a11.startswith('matepairlowstring'):
        if a22.lower().startswith('f'):
            MPLS1=False
        else:
            MPLS1=True
    elif a11.startswith('chromosomebychromosome') :
        if a22.lower().startswith('f'):
            ChromosomeByChromosome1=False
        else:
            ChromosomeByChromosome1=True
    elif a11.startswith('reportbuffers') :
        if a22.lower().startswith('f'):
            ReportBuffers1=False
        else:
            ReportBuffers1=True
    elif a11.startswith('reportread') :
        if a22.lower().startswith('f'):
            ReportReads1=False
        else:
            ReportReads1=True
    elif a11.startswith('reportposition') :
        if a22.lower().startswith('f'):
            ReportPositions1=False
        else:
            ReportPositions1=True
    elif a11.startswith('reportstart') :
        if a22.lower().startswith('f'):
            ReportStarts1=False
        else:
            ReportStarts1=True
    elif a11.startswith('reportend') :
        if a22.lower().startswith('f'):
            ReportEnds1=False
        else:
            ReportEnds1=True
    elif a11.startswith('reportcover') :
        if a22.lower().startswith('f'):
            ReportCoverage1=False
        else:
            ReportCoverage1=True
    elif a11.startswith('reportsens') :
        if a22.lower().startswith('f'):
            ReportSense1=False
        else:
            ReportSense1=True
    elif a11.startswith('reportanti') :
        if a22.lower().startswith('f'):
            ReportAntisense1=False
        else:
            ReportAntisense1=True
    elif a11.startswith('reportboth') :
        if a22.lower().startswith('f'):
            ReportBothStrands1=False
        else:
            ReportBothStrands1=True
    elif a11.startswith('reportuniqu') :
        if a22.lower().startswith('f'):
            ReportUnique1=False
        else:
            ReportUnique1=True
    elif a11.startswith('reportrepeat') :
        if a22.lower().startswith('f'):
            ReportRepeats1=False
        else:
            ReportRepeats1=True
    elif a11.startswith('fullcover') :
        if a22.lower().startswith('f'):
            FullCoverage1=False
        else:
            FullCoverage1=True
            FullCoverageByFile1=False 
    elif a11.startswith('gffstart') or a11.startswith('featurestart') or a11.startswith('gtfstart') :
        ReportStarts1=True
        if a22.lower().startswith('f'):
            ReportStarts1=False
    elif a11.startswith('featureend') or a11.startswith('gffend') or a11.startswith('gtfend') :
        ReportEnds1=True
        if a22.lower().startswith('f'):
            ReportEnds1=False
    elif a11.startswith('featurerepeat') or a11.startswith('gffrepeat') or a11.startswith('gtfrepeat') :
        ReportRepeats1=True
        if a22.lower().startswith('f'):
            ReportRepeats1=False
    elif a11.startswith('briefcover'):
        BriefCoverage1=True
        if a22.lower().startswith('f'):
            BriefCoverage1=False
    elif a11.startswith('metadata'):
        MetaData1=True
        if a22.lower().startswith('f'):
            MetaData1=False
    elif a11.startswith('deleten'):
        DeleteNs1=True
        if a22.lower().startswith('f'):
            DeleteNs1=False
    elif a11.startswith('firstkt'):
        FirstKTable1=True
        if a22.lower().startswith('f'):
            FirstKTable1=False
    elif a11.startswith('interimkt'):
        InterimKTable1=True
        if a22.lower().startswith('f'):
            InterimKTable1=False
    elif a11.startswith('firstkstr') or a11.startswith('firstkbystr') :
        FirstKByStrand1=True
        if a22.lower().startswith('f'):
            FirstKByStrand1=False
    elif a11.startswith('fiveprimeex'):
        FivePrimeExtensionDisallowed1=True
        if a22.lower().startswith('f'):
            FivePrimeExtensionDisallowed1=False
    elif a11.startswith('keepdouble'):
        KeepDoubleExtension1=True
        if a22.lower().startswith('f'):
            KeepDoubleExtension1=False
    elif a11.startswith('findstructure') :
        FindStructuralAnomalies1=True
        if a22.lower().startswith('f'):
            FindStructuralAnomalies1=False
    elif a11.startswith('ucsclinks') :
        UCSCLinks1=True
        if a22.lower().startswith('f'):
            UCSCLinks1=False
    elif a11.startswith('combination'):
        if a22.lower().startswith('a') or a22.lower().startswith('c') or a22.lower().startswith('1'):  ## 'all' or no entry -- record all paired read locations
            CombinationList1=1
        elif a22.lower().startswith('b') or a22.lower().startswith('i') or a22.lower().startswith('2'):  ## use base-by-base positions to record
            CombinationList1=2
            BaseByBase1 = True
        elif a22.lower().startswith('t') or a22.lower().startswith('3'):  ## only record for potential Tn5 insertions with 8-10b duplication
            CombinationList1=3
    elif a11.startswith('requirek') or a11.startswith('kmerrequire') :
        RequireKMerOnly1=True
        if a22.lower().startswith('f'):
            RequireKMerOnly1=False
    elif a11.startswith("verbose"):
        VerboseOutput1=True
        if a22.lower().startswith('f'):
            VerboseOutput1=False
    elif a1.lower().startswith('coveragebyfile'):        
        FullCoverageByFile1=True 
        if a22.lower().startswith('f'):
            FullCoverageByFile1=False
    elif a1.lower().startswith('kmerbincoverageunique'):
        KmerBinCoverageUnique1= True
        if a22.lower().startswith('f'):
            KmerBinCoverageUnique1=False
    elif a1.lower().startswith('kmerbincoverageall'):
        KmerBinCoverageAll1 = True
        if a22.lower().startswith('f'):
            KmerBinCoverageAll1=False
    elif a11.startswith('barcoderequire') or a11.startswith('requirebarcode') :
        if '2' in a11:
            if 'f' in a22.lower():
                BarcodeRequireR2=False
            else:
                BarcodeRequireR2=True
        else:
            if 'f' in a22.lower():
                BarcodeRequireR1=False
            else:
                BarcodeRequireR1=True
    elif a11.startswith('barcode'):
        if '2' in a11:
            BarcodeR2 = a22.upper()
            if BarcodeR2.count('N')==len(BarcodeR2):
                AnyBarcodeR2 = True
            BarcodeLenR2 = len(BarcodeR2)
        else:
            BarcodeR1 = a22.upper()
            if BarcodeR1.count('N')==len(BarcodeR1):
                AnyBarcodeR1 = True
            BarcodeLenR1 = len(BarcodeR1)
    elif a11.startswith('linkerrequire') or a11.startswith('requirelinker')  :
        if '2' in a11:
            if 'true' in a22.lower():
                LinkerRequireR2=True
            elif 'f' in a22.lower():
                LinkerRequireR2=False
        else:
            if 'true' in a22.lower():
                LinkerRequireR1=True
            if 'f' in a22.lower():
                LinkerRequireR1=False
    elif a11.startswith('linker')  :
        if '2' in a11:
            LinkerR2 = a22.upper()
            LinkerLenR2 = len(LinkerR2)
        else:
            LinkerR1 = a22.upper()
            LinkerLenR1 = len(LinkerR1)
    elif a11.startswith('matepairsub') :
        MatePairSubstitutionMax1=int(a22)
    elif a11.startswith('matepairindelcount') :
        MatePairIndelCountMax1=int(a22)
    elif a11.startswith('matepairindellength') :
        MatePairIndelLengthMax1=int(a22)
    elif a11.startswith('splitreadsub') :
        SplitReadSubstitutionMax1=int(a22)
    elif a11.startswith('splitreadindelcount') :
        SplitReadIndelCountMax1=int(a22)
    elif a11.startswith('splitreadindellength') :
        SplitReadIndelLengthMax1=int(a22)
    elif a11.startswith('maxexten') :
        MaxExtension1=int(a22)
    elif a11.startswith('lines'):
        LinesToProcess1=int(a22)
    elif a11.startswith('firstkmax') or a11.startswith('maxfirstk'):
        klen1d=int(a22)
    elif a11.startswith('ucscassembl'):
        UCSCAssembly1=a22
    elif a11.startswith('k'):
        klen1=int(a22)
    elif a11.startswith('multiplicity') :
        MaxMultiplicityReported1=int(a22)
    elif a11.startswith('coverageaccumulationinterval'):
        CoverageAccumulationInterval1=int(a22)
    elif a11.startswith('fullcoverageinterval'):
        FullCoverageInterval1=int(a22)
        FullCoverageByFile1=True
        FullCoverage1=True
    elif a11.startswith('maxcoverhist'):
        MaxCoverHist1=int(a22)
    elif a11.startswith('fullcoveragebyfile'):
        CoverageAccumulationInterval1=int(a22)
        FullCoverageByFile1=True
        FullCoverage1=True
    elif a11.startswith('circularchromosomelimit'):
        CircularChromosomeLimit1=int(a22)
    elif a11.startswith('coveragedtype'):
        dt11=np.dtype(a22.split('.')[-1])
        CoverageDType1=dt11
    elif a11.startswith('indexcon'):
        IndexConstructionBins1=int(a22)
    elif a11.startswith('debug'):
        Debug1=True

if PreREVA1:
    ReportRepeats1 = False  ## No repeat data to report for PreREVA, so repeat columns are not useful
    ReportCoverage1 = False
    
if len(Feature1D)>0:
    Feature1All = False
else:
    Feature1All = True
RefAbbrev1=os.path.basename(RefSeqFile1).split('.')[0]    
if IndexFile1=='default':
    if MaxMultiplicityReported1 == 63:
        IndexFile0='REVAIndex_'+RefAbbrev1+'_k'+str(klen1)+'_b'+str(BitsIndexed1)+'_wv2.9re'
    else:
        IndexFile0='REVAIndex_'+RefAbbrev1+'_k'+str(klen1)+'_b'+str(BitsIndexed1)+'_m'+str(MaxMultiplicityReported1)+'_wv2.9re'
    if Debug1:
        IndexFile0 = 'Debug'+IndexFile0
    if os.path.isfile(IndexFile0):
        IndexFile1=IndexFile0
    else:
        IndexFile1=os.path.join( os.path.dirname(RefSeqFile1) , IndexFile0 ) ## try to find/make the index in the folder with the Reference Sequence File
        if not ( os.path.isfile(IndexFile1) ):
            try:
                TryIndexOpen1=open(IndexFile1,mode='wb')  ##make sure we have the ability to open a file at that location, otherwise default to the current folder
                TryIndexOpen1.close()
            except:
                IndexFile1=IndexFile0
                    
if UCSCAssembly1=='default':
    if ('ws220' in RefSeqFile1.lower()):
        UCSCAssembly1='ce10'
    elif ('ws235' in RefSeqFile1.lower()):
        UCSCAssembly1='ce11'
    elif ('hg38' in RefSeqFile1.lower()):
        UCSCAssembly1='hg38'
if chrbase1=='default':
    chrbase1='chr'

pypy1=False
if 'pypy' in version.lower():
    pypy1=True

if SRRList1:
    for srr1 in SRRList1:
        FiL0.append(srr1.split('#')[0])
        FiN0.append(srr1)
        FiD0.append('R1')
        
if len(FiL0)>0 and not(Mnemonic1):
    for fil0,fin0 in zip(FiL0,FiN0):
        mn0 = os.path.basename(fil0).split('.')[0]
        if '#' in fin0:
            mn0 += fin0.split('#',1)[1]
        mn1 = mn0.replace('R1','R2')
        mn2 = mn0.replace('R2','R1')
        if mn1 in Mnemonic1:
            mi1 = Mnemonic1.index(mn1)
            Mnemonic1[mi1] = mn0.replace('R1','R1R2')
        elif mn2 in Mnemonic1:
            mi2 = Mnemonic1.index(mn2)
            Mnemonic1[mi2] = mn0.replace('R1','R1R2')
        elif not(mn0 in Mnemonic1):
            Mnemonic1.append(mn0)
if Mnemonic1:
    Mnemonic1='_'.join(Mnemonic1)
else:
    Mnemonic1='Diagnostic'
OutFileNameBase='REVA_'+Mnemonic1+'_'+RefAbbrev1+'_'+strftime("D_%m_%d_%y_T_%H_%M_%S",localtime())

LogFile1="LogSummary_"+OutFileNameBase+'.tdv'
def LogNote1(note,LogFileName):
    LogFile=open(LogFileName,mode='a')
    LogFile.write(note+'\t'+'; Time='+"{0:.2f}".format(time()-t0)+' sec'+'\t'+strftime("D_%m_%d_%y_T_%H_%M_%S",localtime())+' \n')
    LogFile.close()
    print(note.split('#')[0].replace('\t',' ').strip(',') + '; Time='+"{0:.2f}".format(time()-t0)+' sec')

def HeaderTranspose(hT1):
    hT2 = hT1.split('\t')
    hT0 = '<!--\tOutput_Key\t\t-->\n'
    hT0 += '<!--\tColumnNumber\tColumnHeader\t-->\n'
    for iT2,nT2 in enumerate(hT2):
        nT2 = nT2.strip()
        if not(nT2.startswith('<!')):
            hT0 += '<!--\t'+str(iT2)+'\t'+nT2+'\t-->\n'
    hT3 = argv
    if MetaData1:
        hT4 = '<!--\tRunning REVA With Parameters\t-->\n'
        for hT5 in hT3:
            if hT5:
                hT4 += '<!--\t'+hT5+'\t\t-->\n'
        hT4 += '<!--\t\t\t-->\n'
    else:
        hT4 =''
    return hT4+hT0+'<!--\t\t\t-->\n'
    
if BaseByBase1:
    BaseByBaseD1 = {} ## key is a chromosome,position [1-based] pair, value is a list of readouts to be joined [string, integer, and real]
    BaseByBaseHeader1 = []
    if (BaseByBaseRegionsFiles1 == []) and not(BaseByBasePositions1) and GFF1:
        BaseByBaseRegionsFiles1 = GFF1
    if BaseByBaseRegionsFiles1 == []:
        bbbF1 = False
        bbbO1 = False
    if BaseByBaseOutFile1 == 'default':
        BaseByBaseOutFile1 = 'BaseByBase_'+OutFileNameBase+'.tdv'
    FullCoverage1 = True
    if bbbC1:
        BaseByBaseHeader1.append('Chromosome__'+RefAbbrev1)
    if bbbP1:
        BaseByBaseHeader1.append('PositionInChromosome__'+RefAbbrev1)
    if bbbF1:
        BaseByBaseHeader1.append('Feature_Name')
    if bbbO1:
        BaseByBaseHeader1.append('PositionInFeature')
    if bbbB1:
        BaseByBaseHeader1.append('Base__'+RefAbbrev1)
    if bbbR1:
        BaseByBaseHeader1.append('RepeatCharacter_centered-kMer_ULCD__'+RefAbbrev1)
    if bbbM1:
        BaseByBaseHeader1.append('Multiplicity_centered-kMer__'+RefAbbrev1)
    if bbbK1:
        BaseByBaseHeader1.append('Sense_Coverage_centered-kMer__'+Mnemonic1)
        BaseByBaseHeader1.append('Asense_Coverage_centered-kMer__'+Mnemonic1)
    if bbbS1:
        BaseByBaseHeader1.append('Sense_Starts__'+Mnemonic1)
        BaseByBaseHeader1.append('Asense_Starts__'+Mnemonic1)
        ReportStarts1 = True
    if bbbE1:
        BaseByBaseHeader1.append('Sense_Ends__'+Mnemonic1)
        BaseByBaseHeader1.append('Asense_Ends__'+Mnemonic1)
        ReportEnds1 = True
    if bbbT1:
        BaseByBaseHeader1.append('Extensions__'+Mnemonic1)
    BaseByBaseHeader1 = '\t'.join(BaseByBaseHeader1)
    
def uidfinder(L,excludeprefix=('WBG',)):
    '''derives a provisional UID from a line in a GTF file, very 'rough' but shoud be useful in getting abbreviated names for features'''
    L1 = L.split('\t')
    ut = ''
    if len(L1)>2:
        ut = L1[2]
    ul = ''
    exon1 = ''
    if len(L1)>8:
        idnamelist = [L1[2]]
        ul = L1[8].split(';')
        for u0 in ul:
            us = u0.strip().split()
            if len(us)>1:
                u1 = us[0].strip()
                u2 = us[1].strip('" ').strip("'")
                if u1.endswith('id') or u1.endswith('name') or u1.endswith('biotype'):
                    NewName = True
                    for idi1,idn1 in enumerate(idnamelist):
                        if u2 in idn1:
                            NewName = False
                            break
                        if idn1 in u2:
                            NewName = False
                            idnamelist[idi1] = u2
                        if idn1.startswith('ENS') and u2.startswith('ENS'):
                            if re.search('G\d',idn1):
                                idnamelist[idi1] = u2
                            if re.search('T\d',idn1) and re.search('E\d',u2):
                                idnamelist[idi1] = u2
                            NewName = False                           
                        if '.' in idn1 and '.' and u2:
                            p1 = idn1.split('.')[0]
                            p2 = u2.split('.')[0]
                            if p1==p2:
                                idnamelist[idi1] = sorted((idn1,u2))[1]
                                NewName = False                           
                    if NewName:
                        for pre1 in excludeprefix:
                            if u2.startswith(pre1):
                                NewName=False
                    if NewName:
                        if not(u2.startswith('protein')):
                            idnamelist.append(u2)
                if u1=='exon_number':
                    exon1 = '.e'+str(u2)
            if exon1 and not exon1 in ''.join(idnamelist):
                numpointID = len([d for d in idnamelist if '.' in d])
                if numpointID==1:
                    for i,d in enumerate(idnamelist):
                        if '.' in d:
                            idnamelist[i] += exon1
        return '/'.join(idnamelist)
    else:
        return ''
    
def FileInfo1(FileID):
    if type(FileID)==str:
        Fn1=FileID
    else:
        Fn1=FileID.name
    s11=os.stat(Fn1)
    return ','.join([Fn1,
                     '#',
                      'Path='+os.path.abspath(Fn1),
                        'Size='+str(s11[6]),
                        'Accessed='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[7])),
                        'Modified='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[8])),
                        'Created='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[9])),
                        'FullInfo='+str(s11)])
LogNote1("Running REVA with parameters:"+' '.join(argv)+' #Python Version'+version,LogFile1)
LogNote1("Python Flavor/Version: "+version,LogFile1)
def LogOpeningFile1(FileID):
    LogNote1("Opening "+FileInfo1(FileID),LogFile1)
def LogClosingFile1(FileID):
    LogNote1("Closed "+FileInfo1(FileID),LogFile1)
def LogRunningFile1(FileID):
    LogNote1("Running "+FileInfo1(FileID),LogFile1)
LogRunningFile1(MyCode1)

VirtualFPath1 = {}
VirtualDtype1 = {}
def npush(Variable0):
    if type(Variable0) ==str:
        VList = [Variable0,]
    else:
        VList = Variable0
    for Variable1 in VList:
        if not(Variable1 in VirtualFPath1):
            fp1 = os.path.join(ScratchPath1,Variable1+'_'+TempFileUID1+'_reva.tmp')
            VirtualFPath1[Variable1] = fp1
            if type(globals()[Variable1]) == np.ndarray:
                VirtualDtype1[Variable1] = globals()[Variable1].dtype
                open(fp1,mode='wb').write(globals()[Variable1])                
            else:
                pf0 = open(fp1, mode='wb')
                pp0 = cPickle.Pickler(pf0)
                pp0.dump(globals()[Variable1])
                pf0.close()
                VirtualDtype1[Variable1] ='pickle'
        globals()[Variable1] = None

def npull(Variable0):
    if type(Variable0) ==str:
        VList = [Variable0,]
    else:
        VList = Variable0
    for Variable1 in VList:
        fp1 = VirtualFPath1[Variable0]
        dt1 = VirtualDtype1[Variable0]
        if  type(dt1)==str:
            pf0 = open(fp1, mode='rb')
            pp0 = cPickle.Unpickler(pf0)
            globals()[Variable1] = pp0.load()
            pf0.close()
        else:
            globals()[Variable1] =  np.frombuffer(open(fp1,mode='rb').read(),dtype=VirtualDtype1[Variable1])

def GCContent1(s):
    GCount = s.count('G')
    ACount = s.count('A')
    TCount = s.count('T')
    CCount = s.count('C')
    GATCCount = GCount+ACount+TCount+CCount
    if GATCCount==0:
        return 0.5
    else:
        return float(GCount+CCount)/GATCCount
    
if 'darwin' in platform:
    LogNote1('Trying to Run \'caffeinate\' on MacOSX to prevent the system from sleeping',LogFile1)
    try:
        Coffee_process1=subprocess.Popen('caffeinate')
    except:
        LogNote1("Couldn't start 'caffeinate', you may need to manually set your mac to avoid dozing (System Preferences, Energy)",LogFile1)
else:
    LogNote1('Trying to Run \'caffeine\' to prevent the system from sleeping',LogFile1)
    try:
        Coffee_process1=subprocess.Popen('caffeine')
    except:
        LogNote1("Couldn't start 'caffeine', you may need to manually set your mac to avoid dozing (System Preferences, Energy)",LogFile1)
        
   
BitsSorted1=klen1*2-BitsIndexed1
Bins1=2**BitsIndexed1
BitMask1=2**BitsSorted1-1

AllBase1=['G','A','T','C']
Numbase1={AllBase1[0]:0,AllBase1[1]:1,AllBase1[2]:2,AllBase1[3]:3,'N':0}
## e4 is a set of 4**x exponents
e4=[]
ne4=[]
eList32=(1,2,4,8,16,32)
for i in range(32):
    e4.append(4**i)
    ne4.append(4**i-1)
ex2=np.array([2**i for i in range(64)],dtype=np.uint64)
def means1(ListN0,ListN1):
    M1=[]
    for g0,g1 in zip(ListN0,ListN1):
        if g0==0:
            M1.append(0.0)
        else:
            M1.append((1.0*g1)/g0)
    return M1
def sdeviations1(ListN0,ListN1,ListN2):
    StD1=[]
    for g0,g1,g2 in zip(ListN0,ListN1,ListN2):
        if g0==0:
            StD1.append(0.0)
        else:
            StD1.append( ( (1.0*g2)/g0-(1.0*g1*g1)/(g0*g0) )**0.5)
    return StD1
def str02(ListF1):
    FT1=[]
    for r1 in ListF1:
        FT1.append('{0:.2f}'.format(r1))
    return FT1
def per03(ListF1):
    FT1=[]
    for r1 in ListF1:
        FT1.append('{0:.3f}'.format(100.0*r1))
    return FT1
def prependtotal(ListF1):
    return list(map(int,[sum(ListF1),]+ListF1))
def prependtotal2D(ListF1):
    Totals1=[]
    for i11 in xrange(len(ListF1[0])):
        Totals1.append(int(sum([ Fl[i11] for Fl in ListF1 ])))
    return [Totals1,]+ListF1

## Purge Files from a directory
def PurgeREVATempFiles(dir1):
    LogNote1('Starting Purge of Unused Files from Temporary Directory '+ScratchPath1,LogFile1)
    try:
        for f1 in os.listdir(dir1):
            if f1.endswith('_'+TempFileUID1+'_reva.tmp') or f1.endswith('_'+TempFileUID1+'.fasta.gz'):
                os.remove(os.path.join(ScratchPath1,f1))
        LogNote1('Finished Purge of Unused Files from Temporary Directory '+ScratchPath1,LogFile1)
    except Exception as e:
        LogNote1('Purge of temporary files may have failed.  Reason='+str(e)+' --- you may want to do this manually when program is finished.  Directory='+ScratchPath1,LogFile1)

## Open and unpack the index of unique sequences (three sorted lists of integers with 2*klen1 bit list of values-converted-to binary, sorted, then two lists that are the source sequence (chromosome for C. elegna) and position
## fastfind is the routine that returns the index for any given binary sequence k-mer representation
AntisenseB1={'A':'T','T':'A','G':'C','C':'G'}
Tr2=''  ## Tr2 allows a quick translation of sequence to a numerical array.
for i in range(256):
    if chr(i) in 'AGCTagct':
        Tr2+=chr(Numbase1[chr(i).upper()])
    else:  ##in case there are unusual characters in the sequence, they will be converted to Gs
        Tr2+=chr(Numbase1['N'])
TrN1=''  ## Tr2 allows a quick translation of sequence to a numerical array.
for i in range(256):
    if chr(i) in 'AGCTagct':
        TrN1+=chr(1)
    else:  ##in case there are unusual characters in the sequence, they will be converted to Gs
        TrN1+=chr(0)
def seqnum(s):
    ''' converts an input sequence s (upper case string) into a numpy array of values from 0 to 3.  Translation between sequence and numbers from AllBase1'''
    return np.frombuffer(s.translate(Tr2).encode('ascii'),dtype=np.uint8)
def seqNotN(s):
    ''' converts an input sequence s (upper case string) into a numpy array of values from 0 to 3.  Translation between sequence and numbers from AllBase1'''
    return np.frombuffer(s.translate(TrN1).encode('ascii'),dtype=np.bool)

## antisense-- returns the reverse complement of argument string 
filterminus=''
ASB11="AaCcNn*nNgGtT"
for i in range(256):
    if chr(i) in ASB11:
        filterminus+=ASB11[12-ASB11.find(chr(i))]
    else:
        filterminus+='N'  ## switched from filterminus+='' 5/08/21 to avoid crashes from antisense and sense sequences having different lengths
def antisense(s):
    '''return an antisense and filtered version of any sequence)'''
    return s.translate(filterminus)[::-1]

BitDepthD1={np.dtype('bool'):1,
    np.dtype('int8'):1,
    np.dtype('int16'):2,
    np.dtype('int32'):4,
    np.dtype('int64'):8,
    np.dtype('uint8'):1,
    np.dtype('uint16'):2,
    np.dtype('uint32'):4,
    np.dtype('uint64'):8,
    np.dtype('float16'):2,
    np.dtype('float32'):4,
    np.dtype('float64'):8,
    np.dtype('complex64'):8,
    np.dtype('complex128'):16}

def AltFromFile(F_file,dtype,count):  ## PyPy Numpy has no .fromfile method for arrays, so I made this very simple equivalent for 1-D arrays.  Could fail if the index was written on another platform 
    ## F_file needs to be open with 'wb' mode, 'dtype is a numpy data type, and count is the number of entries to read)
    h=F_file.read( BitDepthD1[dtype] * count )
    return np.frombuffer(h,dtype=dtype)
def AltToFile (Array1,F_file):  ## PyPy Numpy has no .fromfile method for arrays, so I made this very simple equivalent for 1-D arrays.  Could fail if the index was written on another platform 
    ## F_file needs to be open with 'wb' mode, Array1 is a numpy array.  May fail if file is written on one platform and read on another)
    F_file.write(Array1)

dtype1=np.uint16
if klen1>8:
    dtype1=np.uint32
if klen1>16:
    dtype1=np.uint64

b_a1 = 12  ## bit length for individual counts of events.  2**b1 is the maximum sequence length for analysis
         ## for anything like 'current' technology (2018) the implicit maximum of 4096 bases per read seems quite adequate
b_a2 = b_a1*2
b_a3 = b_a1*3
b_a4 = b_a1*4
b_a5 = b_a1*5
def unpackB(v0):
    'output parsed list (score,match,mismatch,indelS,indelE) from input of raw alignment value,  used by nw1 to deliver integral values of score, matches, mismatches, etc'
    ## This is designed to work with sequences up to 1kb with up to 99 distinct insertions/deletions
    score = v0 >> b_a4
    v1 = v0-(score<<b_a4)
    match = v1>>b_a3
    v2 = v1-(match<<b_a3)
    mismatch = v2>>b_a2
    v3 = v2-(mismatch<<b_a2)
    indelE = v3>>b_a1
    indelS = v3-(indelE<<b_a1)
    return [score,match,mismatch,indelS,indelE]

def SimpleAlign1(s01,s02):
    '''S01 is a short string to decorate with upper case at matches, second is a reference, p02 is one-based; + for sense, - for antisense; output is an alignment string with caps at matched positions followed by number of matches and mismatches'''
    salign1 = ''
    match1 = 0
    mismatch1 = 0
    for b01,b02 in izip_longest(s01,s02,fillvalue='X'):
        if b01.upper()==b02.upper():
            salign1 += b01.upper()
            match1+=1
        else:
            salign1 += b01.lower()
            mismatch1+=1
    return salign1, match1, mismatch1
    
## Match and mismatch score values.  These can be changes as needed but won't affect that much in this program
matchV=1  ## reward for single base match (default is 1)
mismatchV=-1  ## "reward" for single base mismatch (default is -1.. so a penalty)
indelSV=-3  ## "reward" for starting a gap (default is -3)
indelEV=-1  ## "reward" for extending a gap (default is -3)

maxLa1 = 2**b_a1  ## Maximum length for aligned segment-- default is 2**12==4096, will truncate beyond that to avoid troubles
scoreI = 2**b_a4
matchI = matchV*scoreI+2**b_a3
mismatchI = mismatchV*scoreI+2**b_a2
indelEI = indelEV*scoreI+2**b_a1
indelSI = indelSV*scoreI+1
defaultStart = -2**(b_a5-1)
def nw1(s1,s2,span=12):  ## s1 and s2 are sequence strings or byte arrays.  No case checking done so all upper or lower to shart
    '''Alignment Quality Check-- input is two sequences plus a 'span' distance, output is a tuple of score, match number, mismatch number, gap number, total gap length'''
    ## The alignment differs in a few ways from standard Smith=Waterman/Needleman-Wunsch.
    ## First, only base positions within a range of +/- span from each other are considered for alignment
    ## so for span=1, a base is only considered for alignment to the previous or next base.  Thus any gaps can't be longer than length span
    ## and net gap length is similarly limited
    ## Second, no penalty is given for end-gaps within the range span.  So if span=3 the total net allowed gap distance is 3 in either direction at any point-- but gaps at the beginning or end count toward this total but not in the eventual score or number of reported gaps.  This avoids
    ## an issue where a gap can be double-penalized if a SW alignment algorithm insists on both termini matching
    ## Third, the open-gap penalty is larger than the etend gap penalty, with the distinction being dependent on
    ## the highest current value for alignment for the previous base-pair combination.  For a small
    ## number of cases, this results in not properly closing two adjacent deletions which might by SW alignment be otherwise
    ## preferred.  Allowing an additional gap should avoid most problems with this.
    ## the subroutine is coded in pure python (not Numpy) as this was faster in all tests.
    ## the seuqential nature of NW/SW algorithm prevents bulk operations facilitated by Numpy and numpy otherwise imposes a substantial overhea
    ## PyPy runs this code really really fast, as would recoding in C using Cython
    l1 = min(maxLa1,len(s1))
    l2 = min(maxLa1,len(s2))
    X = l1+1
    Y = 2*span+1
    d1 = [0]*(X*Y) ## d1[(1+n1)*Y+d1] is the highest alignment value for bases zero to n1 of S1 with bases zero to (1+n1+(d1-span-1) of S2.  It is also the value such an alignment contributes to a junction score from the S1 side (thorugh position n1) if the junctino point is set just after position n1 in S1 (zero based)
    e1 = [0]*(X*Y)  ## programmer's note.. you are never going to believe me, but this is really faster with a native Python data structure and direct comparisons (rather than max statements) than with Numpy.  It really is.  You can recode it with Numpy and it will run, but just 20x slower.   
    vMax = defaultStart
    j1 = -span-1
    i1 = 0
    for i2 in xrange(Y,X*Y):  ##i1+1 is the position in the sequence, i1 is the position in array
        j1 += 1
        if j1 > span:
            j1 = -span
            i1 += 1
        if j1+i1<0 or j1+i1>=l2:
            continue
        if s1[i1]==s2[j1+i1]:
            f1 = d1[i2-Y]+matchI
        else:
            f1 = d1[i2-Y]+mismatchI
        if j1 ==-span:
            f2 = defaultStart
        elif e1[i2-1]==2:
            f2 = d1[i2-1] + indelEI
        else:
            f2 = d1[i2-1] + indelSI
        if j1 == span:
            f3 = defaultStart
        elif e1[i2-Y+1]==3:
            f3 = d1[i2-Y+1] + indelEI
        else:
            f3 = d1[i2-Y+1] + indelSI
        if (f1>f2 and f1>f3):
            d1[i2] = f1
            e1[i2] = 1                                
        elif f2>f3:
            d1[i2] = f2
            e1[i2] = 2                
        else:
            d1[i2] = f3
            e1[i2] = 3
        if (d1[i2]>vMax) and (i1==l1-1 or i1+j1==l2-1):
            vMax = d1[i2]
    return unpackB(vMax)

MultiplicityDataType1=np.uint8
Local1 = 2**6  ## coding for local repeats
Chromosomal1 = 2**7 ## coding for global repeats 
if MaxMultiplicityReported1>63:
    MultiplicityDataType1=np.uint16
    Local1 = 2**14  ## coding for local repeats
    Chromosomal1 = 2**15 ## coding for global repeats 
if MaxMultiplicityReported1>(2**14-1):
    MultiplicityDataType1=np.uint32
    Local1 = 2**30  ## coding for local repeats
    Chromosomal1 = 2**31 ## coding for global repeats 
if MaxMultiplicityReported1>(2**30-1) or (MaxMultiplicityReported1<=0):
    MultiplicityDataType1=np.uint64
    Local1 = 2**62  ## coding for local repeats
    Chromosomal1 = 2**63 ## coding for global repeats
Mult1 = Local1-1  ## Mask for bits relevant to multiplicity (loses two highest bits, of Multiplicity array, which are used for Local and Chromosomal respectively)
    
LogNote1('starting fastA file read',LogFile1)
LogOpeningFile1(RefSeqFile1)
if RefSeqFile1.endswith('.gz'):
    if version.startswith('2.'):
        F0=gzip.open(RefSeqFile1,mode='r')
    else:
        F0=gzip.open(RefSeqFile1,mode='rt')
elif RefSeqFile1.endswith('.zip'):
    ZipFile1 = zipfile.ZipFile(RefSeqFile1)  ##.readline()
    F0 = []
    for zfn1 in ZipFile1.namelist():
        if zfn1[0]!='_':
            F0.extend(ZipFile1.open(zfn1,mode='r').read().splitlines())
else:
    if version.startswith('2.'):
        F0=open(RefSeqFile1,mode='rU')
    else:
        F0=open(RefSeqFile1,mode='r')
SD1=[] ## Sense Sequences
AD1=[] ## AntiSense Sequences
LD1=[] ## Length of each sequence
NameA1=[] ## Names of each reference DNA entity by number
NameD1={} ## Numbers for each reference DNA entity by name
CircD1=[] ## Cicular status of each input DNA entity by number if anything less than 13.5MB is cir
for L0 in F0:
    L1=L0.strip()
    if L1.startswith('>'):
        if len(SD1)>0 and Debug1:  ## DEbug mode-- only look at first chromosomal unit in FastA reference file
            break
        NameA1.append(L1[1:].strip().split()[0])
        if "circular" in L1.lower():
            CircD1.append(True)
        elif "linear" in L1.lower():
            CircD1.append(False)
        else:
            CircD1.append(-1)            
        sn1=NameA1[-1]
        NameD1[sn1]=len(NameA1)-1
        SD1.append([])
    else:
        SD1[-1].append(L1.upper().replace('U','T'))
for sn1 in range(len(SD1)):
    SD1[sn1]=''.join(SD1[sn1])   ##SD1[chr][1] is the sense sequence of chr1, SD1.  N's will be converted to "G" to preserve absolute position.
    AD1.append(antisense(SD1[sn1]))
    LD1.append(len(SD1[sn1]))
    if LD1[sn1]<CircularChromosomeLimit1:
        if CircD1[sn1] == -1: 
            CircD1[sn1] = True
    else:
        if CircD1[sn1] == -1: 
            CircD1[sn1] = False
    if CircD1[sn1]:
        SD1[sn1]=SD1[sn1]+SD1[sn1][:klen1-1]
        AD1[sn1]=AD1[sn1]+AD1[sn1][:klen1-1]

TotalRefBases1=sum(LD1)
LongestContig1=max(LD1)
dtype2=np.int32  ## data type (signed) for any reference to positions (sense + antisense-) in the individual chromosomes 
if LongestContig1>=2**31-2*klen1-2*Tn5DupMax1:  ## 2*klen1+2*Tn5DupMax1 is a "safety factor to make sure there is no overflow
    dtype2=np.int64
dtype4=np.uint16
if BitsSorted1>16:
    dtype4=np.uint32
if BitsSorted1>32:
    dtype4=np.uint64
dtype5=np.uint16
if BitsIndexed1>16:
    dtype5=np.uint32
if BitsIndexed1>32:
    dtype5=np.uint64

NumSeq1=len(NameA1)
NameARange0=range(NumSeq1)
NameARange1=range(NumSeq1+1)
if not(RefSeqFile1.endswith('zip')):
    F0.close()
    LogClosingFile1(F0)

if IndexConstructionBins1 == 0:
    if TotalRefBases1<120000000:
        IndexConstructionBins1 = 1
    elif TotalRefBases1<400000000:
        IndexConstructionBins1 = 4
    elif TotalRefBases1<1600000000:
        IndexConstructionBins1 = 16
    else:
        IndexConstructionBins1 = 64        
qssD = 1+ (4**klen1//IndexConstructionBins1)

LogNote1('finished fastA file read',LogFile1)
LogNote1('Getting ready to try unpickle/unpack of Index File: '+IndexFile1,LogFile1)
try:
    pf1=open(IndexFile1,mode='rb')
    LogOpeningFile1(pf1)
    pp1=cPickle.Unpickler(pf1)
    RefSeqFileR1=pp1.load()
    klenR1=pp1.load()
    ulenR1=pp1.load()
    V1dtype=pp1.load()
    P1dtype=pp1.load()
    I1dtype=pp1.load()
    NameA1=pp1.load()
    MultiplicityDataType1=pp1.load()
    MaxMultiplicityReported1=pp1.load()
    Local1=pp1.load()
    Chromosomal1=pp1.load()
    FastIndexDataType1=pp1.load()
    TotalUnique1=pp1.load()
    TotalLocal1=pp1.load()
    TotalChromosomal1=pp1.load()
    TotalDispersed1=pp1.load()
    if pypy1:
        Sort_V1dedup=AltFromFile(pf1,dtype=V1dtype,count=ulenR1)
        Sort_P1dedup=AltFromFile(pf1,dtype=P1dtype,count=ulenR1+1)
        Sort_I1dedup=AltFromFile(pf1,dtype=I1dtype,count=ulenR1+1)
        Multiplicities1=AltFromFile(pf1,dtype=MultiplicityDataType1,count=ulenR1+1)
        FastIndex1=AltFromFile(pf1,dtype=FastIndexDataType1,count=2**BitsIndexed1+1)
    else:
        Sort_V1dedup=np.fromfile(pf1,dtype=V1dtype,count=ulenR1)
        Sort_P1dedup=np.fromfile(pf1,dtype=P1dtype,count=ulenR1+1)
        Sort_I1dedup=np.fromfile(pf1,dtype=I1dtype,count=ulenR1+1)
        Multiplicities1=np.fromfile(pf1,dtype=MultiplicityDataType1,count=ulenR1+1)
        FastIndex1=np.fromfile(pf1,dtype=FastIndexDataType1,count=2**BitsIndexed1+1)
    pf1.close()
    LogClosingFile1(pf1)
    LogNote1('Finished unpicle/unpack',LogFile1)
except Exception as e:
    LogNote1('Exception encountered in finding/reading the index file from disk: '+str(e),LogFile1)
    LogNote1('Creating the index from scratch (takes extra time but not a problem otherwise)',LogFile1)
    ## IndexConstructionBins1 = How many bins to divide data in for sorting
    ## The code below is a very rough attempt to avoid suboptimal setting.  It could be tweaked in various ways,
    ## or IndexConstructionBins1 can be set from the command line.
    qssT = 0
    if IndexConstructionBins1>1:
        IdF1 = [open(os.path.join(ScratchPath1,'Id_'+str(i)+'_'+TempFileUID1+'_reva.tmp'),mode='wb') for i in xrange(IndexConstructionBins1)]
        ValF1 = [open(os.path.join(ScratchPath1,'Val_'+str(i)+'_'+TempFileUID1+'_reva.tmp'),mode='wb') for i in xrange(IndexConstructionBins1)]
        PosF1 = [open(os.path.join(ScratchPath1,'Pos_'+str(i)+'_'+TempFileUID1+'_reva.tmp'),mode='wb') for i in xrange(IndexConstructionBins1)]
    else:
        IdTx1 = []
        ValTx1 = []
        PosTx1 = []
    if NumSeq1<256:
        IdAdtype=np.uint8
    elif NumSeq1<65536:
        IdAdtype=np.uint16
    elif NumSeq1<4294967296:
        IdAdtype=np.uint32
    else:
        IdAdtype=np.uint64 
    for Ni1 in NameARange0:
        sL1=len(SD1[Ni1]) ## number of bases for analysis: featureLen+klen1-1 for a circular feature, featureLen for a linear feature
        sL2=sL1-klen1+1 ##number of k-mers that can be extracted
        sL3=sL1-(CircD1[Ni1]*(klen1-1))  ## the actual number of bases in the feature
        s=np.zeros(sL1,dtype=dtype1)+seqnum(SD1[Ni1])
        a=np.zeros(sL1,dtype=dtype1)+seqnum(AD1[Ni1])
        n=seqNotN(SD1[Ni1])
        if sL2>=0:
            sA1=np.zeros(sL2,dtype=dtype1)
            sN1=np.ones(sL2,dtype=np.bool)
            aA1=np.zeros(sL2,dtype=dtype1)
            sP1=np.arange(1,sL2+1,dtype=dtype2)        
            aP1=np.arange(-sL3,sL2-sL3,dtype=dtype2)
            i1=0
            for kin1 in eList32:
                if kin1>klen1:
                    break
                if klen1 & kin1:                   
                    if i1==0:
                        sA1 += s[klen1-kin1:]
                        aA1 += a[klen1-kin1:]
                        sN1 &= n[klen1-kin1:]
                    else:
                        sA1 += s[klen1-kin1-i1:-i1]*e4[i1]
                        aA1 += a[klen1-kin1-i1:-i1]*e4[i1]
                        sN1 &= n[klen1-kin1-i1:-i1]
                    i1+=kin1
                if kin1<32:                
                    s = s[:-kin1]*e4[kin1]+s[kin1:]
                    a = a[:-kin1]*e4[kin1]+a[kin1:]
                    n = n[:-kin1] & n[kin1:]
            NFilteredLength1 = np.sum(sN1)
            ValA1x = np.empty(NFilteredLength1*2, dtype=dtype1)
            PosA1x = np.empty(NFilteredLength1*2, dtype=dtype2)                
            ValA1x[0::2] = sA1[sN1]
            PosA1x[0::2] = sP1[sN1]
            if CircD1[Ni1]:  ## added to properly permute the junctional k-mers in the antisense orientation for proper merger with sense elements
                aA2 = np.concatenate( [aA1[-(klen1-1):] , aA1[:-(klen1-1 )]] )
                aP2 = np.concatenate( [aP1[-(klen1-1):] , aP1[:-(klen1-1 )]] )
                ValA1x[1::2] = aA2[::-1][sN1]
                PosA1x[1::2] = aP2[::-1][sN1]
            else:
                ValA1x[1::2] = aA1[::-1][sN1]
                PosA1x[1::2] = aP1[::-1][sN1]
            IdA1x=np.full(NFilteredLength1*2,Ni1,dtype=IdAdtype)
            if IndexConstructionBins1>1:
                for qssi in xrange(IndexConstructionBins1):
                    if qssi==0:
                        qssMask = np.nonzero( ValA1x<(qssi+1)*qssD )[0]
                    elif qssi<IndexConstructionBins1-1:
                        qssMask = np.nonzero( (ValA1x>=qssi*qssD) & (ValA1x<(qssi+1)*qssD) )[0]
                    else:
                        qssMask = np.nonzero( (ValA1x>=qssi*qssD) )[0]
                    IdF1[qssi].write(IdA1x[qssMask])
                    ValF1[qssi].write(ValA1x[qssMask])
                    PosF1[qssi].write(PosA1x[qssMask])
            else:
                IdTx1.append(IdA1x)
                ValTx1.append(ValA1x)
                PosTx1.append(PosA1x)         
            qssT += 2*len(sP1)
    s=None;  sA1=None; a=None; aA1=None; sP1=None; aP1=None ##Manual Garbage Collection here and below may help avoid memory overuse
    IdA1x = None; ValA1x = None; PosA1x =  None
    LogNote1('Initial Conversion completed; '+str(klen1)+'-mer positions='+str(qssT),LogFile1)
    if IndexConstructionBins1>1:
        ValA1x=None; PosA1x=None; IdA1x=None
        for qssf1 in ValF1+PosF1+IdF1:
            qssf1.close()
        IdF2 = open(os.path.join(ScratchPath1,'IdAll'+'_'+TempFileUID1+'_reva.tmp'),mode='wb')
        ValUF2 = open(os.path.join(ScratchPath1,'ValU'+'_'+TempFileUID1+'_reva.tmp'),mode='wb')
        ValLF2 = open(os.path.join(ScratchPath1,'ValL'+'_'+TempFileUID1+'_reva.tmp'),mode='wb')
        PosF2 = open(os.path.join(ScratchPath1,'PosAll'+'_'+TempFileUID1+'_reva.tmp'),mode='wb')
        MultF2 = open(os.path.join(ScratchPath1,'Mult'+'_'+TempFileUID1+'_reva.tmp'),mode='wb')
    else:
        IdTx1 = np.concatenate(IdTx1)
        ValTx1 = np.concatenate(ValTx1)
        PosTx1 =  np.concatenate(PosTx1)                
    TotalKmers1 = 0
    for qssi in xrange(IndexConstructionBins1):
        if IndexConstructionBins1>1:
            IdTx1 = np.frombuffer(open(os.path.join(ScratchPath1,'Id_'+str(qssi)+'_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=IdAdtype)
            ValTx1 = np.frombuffer(open(os.path.join(ScratchPath1,'Val_'+str(qssi)+'_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=dtype1)
            PosTx1 = np.frombuffer(open(os.path.join(ScratchPath1,'Pos_'+str(qssi)+'_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=dtype2)
        Sort_i1 = np.argsort(ValTx1,kind='mergesort')
        LogNote1('Sort completed: '+str(qssi),LogFile1)
        ValTx1 = ValTx1[Sort_i1]
        LogNote1('Array ValTx1 rearranged: '+str(qssi),LogFile1)
        PosTx1 = PosTx1[Sort_i1]
        LogNote1('Array PosTx1 rearranged: '+str(qssi),LogFile1)
        IdTx1 = IdTx1[Sort_i1]
        LogNote1('Array IdTx1 rearranged: '+str(qssi),LogFile1)
        u1=np.concatenate([[True],np.not_equal(ValTx1[1:],ValTx1[:-1])])  ##true at first occurence of each unique k-mer    
        Sort_V1dedup=ValTx1[u1]
        ValTx1=None
        Sort_P1dedup=PosTx1[u1]
        Sort_I1dedup=IdTx1[u1]
        u1=np.concatenate([u1,[True]])  ##true at last occurence of each k-mer
        LogNote1('Deduplication completed: '+str(qssi),LogFile1)
        Multiplicities1=np.nonzero(u1)[0]
        Multiplicities1=Multiplicities1[1:]-Multiplicities1[:-1]
        Multiplicities1=np.minimum(MaxMultiplicityReported1,Multiplicities1)
        Multiplicities1=Multiplicities1.astype(MultiplicityDataType1)
        LogNote1('Multiplicities calculated: '+str(qssi),LogFile1)
        u1=u1[1:]
        Chromosomal1a = (Multiplicities1>1) & np.equal(Sort_I1dedup,IdTx1[u1])## Boolean Array that is True if all instances are on the same chromosome ( or to be precise, the same entry in the original fasta reference file )
        Sort_I1=None
        Local1a = Chromosomal1a & ( ( abs(PosTx1[u1])-abs(Sort_P1dedup) )<=DomainLength1)  ## Boolean Array that is True if all instances of a k-mer are within MaxLocusSpan1 distance on the same chromosome
        Chromosomal1a = (Chromosomal1a & ~Local1a)
        Multiplicities1 = Multiplicities1 + (Chromosomal1a * Chromosomal1).astype(MultiplicityDataType1) + (Local1a * Local1).astype(MultiplicityDataType1)
        u1=None; Sort_P1=None
        LogNote1('Span test completed: '+str(qssi),LogFile1)
        TotalKmers1 += len(Sort_V1dedup)
        if IndexConstructionBins1>1:
            IdF2.write(Sort_I1dedup)
            PosF2.write(Sort_P1dedup)
            MultF2.write(Multiplicities1)
            ValUF2.write((Sort_V1dedup>>BitsSorted1).astype(dtype5))
            ValLF2.write((Sort_V1dedup & BitMask1).astype(dtype4))
            Sort_I1dedup = None
            Sort_P1dedup = None
            Sort_V1dedup = None
            Multiplicities1 = None
            LogNote1('TempFile writes completed: '+str(qssi),LogFile1)
    LogNote1('All deduplication completed; '+str(klen1)+'-mer positions='+str(TotalKmers1),LogFile1)
    LogNote1('TerminalCatenateStart',LogFile1)
    if IndexConstructionBins1>1:
        MultF2.write(np.zeros(1,dtype=MultiplicityDataType1))
        IdF2.write(np.array([NumSeq1],dtype=IdAdtype))
        PosF2.write(np.zeros(1,dtype=dtype2))
        IdF2.close(); ValUF2.close(); ValLF2.close(); PosF2.close(); MultF2.close()
    else:
        Multiplicities1 = np.concatenate((Multiplicities1,np.zeros(1,dtype=MultiplicityDataType1)))
        Sort_I1dedup = np.concatenate((Sort_I1dedup, np.array([NumSeq1],dtype=IdAdtype)))
        Sort_P1dedup = np.concatenate((Sort_P1dedup, np.zeros(1,dtype=dtype2)))
    LogNote1('TerminalCatenateEnd',LogFile1)
    
    LogNote1('FastIndex assembly starting',LogFile1)
    FastIndexDataType1=np.uint32  ## data type (unsigned) for any reference to positions in the sorted arrays (FastIndex) 
    if TotalKmers1>=2**32-2:
        FastIndexDataType1=np.uint64
    if IndexConstructionBins1>1:
        u3 = np.frombuffer(open(os.path.join(ScratchPath1,'ValU_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=dtype5)
    else:
        u3 = (Sort_V1dedup>>BitsSorted1).astype(dtype5)
    FastIndex1=np.zeros(2**BitsIndexed1,dtype=FastIndexDataType1) ## FastIndex1[ix1] is 1+the first position in the Sort_V1dedup array where Sort_V1dedup[iy1] >> BitsSorted1>=ix1  (zero if never ==ix1)
    u4=np.concatenate(([True],np.not_equal(u3[:-1],u3[1:]),[True]))
    u5=np.nonzero(u4)[0]
    u6=(u5[1:]-u5[:-1]).astype(FastIndexDataType1)
    u5=None
    FastIndex1[u3[u4[1:]]]=u6
    u3=None;u4=None;u6=None
    FastIndex1=np.cumsum(np.concatenate([np.zeros(1,dtype=FastIndex1.dtype),FastIndex1]))
    LogNote1('FastIndex completed',LogFile1)
    if IndexConstructionBins1>1:
        Sort_V1dedup=np.frombuffer(open(os.path.join(ScratchPath1,'ValL_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=dtype4)
        Sort_P1dedup=np.frombuffer(open(os.path.join(ScratchPath1,'PosAll_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=dtype2)
        Multiplicities1=np.frombuffer(open(os.path.join(ScratchPath1,'Mult_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=MultiplicityDataType1)
        Sort_I1dedup=np.frombuffer(open(os.path.join(ScratchPath1,'IdAll_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=IdAdtype)
        PurgeREVATempFiles(ScratchPath1)
    else:
        Sort_V1dedup = (Sort_V1dedup & BitMask1).astype(dtype4)
    LogNote1("Starting to calculate Coverage Overhead",LogFile1)
    TotalUnique1=[]
    TotalLocal1=[]
    TotalChromosomal1=[]
    TotalDispersed1=[]
    for z in range(NumSeq1):
        TotalUnique1.append( np.count_nonzero( (Sort_I1dedup==z) & (Multiplicities1==1))//2 )
        TotalLocal1.append( np.count_nonzero( (Sort_I1dedup==z) & ((Multiplicities1 & Local1)>0 ))//2 )
        TotalChromosomal1.append( np.count_nonzero( (Sort_I1dedup==z) & ((Multiplicities1 & Chromosomal1)>0 ))//2 )
        TotalDispersed1.append( np.count_nonzero( (Sort_I1dedup==z) & (Multiplicities1>1) & (Multiplicities1<Local1) )//2 )
        ## Note that for even k-mer lengths, the repeated total counts (Local, Chromosomal, and Dispersed) may be undercounted by approximately half the number of repeats in each category that are themselves perfect palindromes.  Fixing this would require considerable CPU time and not really provide substantial value...  This is not an issue for odd k-mer lengths
    TotalUnique1=prependtotal(TotalUnique1)
    TotalLocal1=prependtotal(TotalLocal1)
    TotalChromosomal1=prependtotal(TotalChromosomal1)
    TotalDispersed1=prependtotal(TotalDispersed1)
    LogNote1("Finished Calculating Coverage Overhead",LogFile1)
    V1dtype=Sort_V1dedup.dtype
    P1dtype=Sort_P1dedup.dtype
    I1dtype=Sort_I1dedup.dtype
    
    try:
        ## the following lines store the index as a disk file (making startup somewhat faster)
        LogNote1('Will try to write the index to disk',LogFile1)
        pf1=open(IndexFile1,mode='wb')
        LogOpeningFile1(pf1)
        pp1=cPickle.Pickler(pf1)
        pp1.dump(RefSeqFile1)
        pp1.dump(klen1)
        pp1.dump(len(Sort_V1dedup))
        pp1.dump(V1dtype)
        pp1.dump(P1dtype)
        pp1.dump(I1dtype)
        pp1.dump(NameA1)
        pp1.dump(Multiplicities1.dtype)
        pp1.dump(MaxMultiplicityReported1)
        pp1.dump(Local1)
        pp1.dump(Chromosomal1)
        pp1.dump(FastIndex1.dtype)
        pp1.dump(TotalUnique1)
        pp1.dump(TotalLocal1)
        pp1.dump(TotalChromosomal1)
        pp1.dump(TotalDispersed1)

        if pypy1:
            AltToFile(Sort_V1dedup,pf1)
            AltToFile(Sort_P1dedup,pf1)
            AltToFile(Sort_I1dedup,pf1)
            AltToFile(Multiplicities1,pf1) 
            AltToFile(FastIndex1,pf1) 
        else:
            Sort_V1dedup.tofile(pf1)
            Sort_P1dedup.tofile(pf1)
            Sort_I1dedup.tofile(pf1)
            Multiplicities1.tofile(pf1) 
            FastIndex1.tofile(pf1)
        pf1.close()
        LogClosingFile1(pf1)
        LogNote1('Success in writing index to disk',LogFile1)
    except Exception as e:
        LogNote1('Exception Encountered writing index to disk: '+str(e),LogFile1)
        LogNote1('****IMPORTANT: DUE TO INCOMPLETE WRITE, YOU WILL NEED TO DELETE FILE '+IndexFile1+' BEFORE RERUNNING PROGRAM',LogFile1)
        LogNote1('Unable to write indexes to disk (program may run but will need to assemble indexes next time',LogFile1)
ulen1=len(Sort_V1dedup)
NameA1.append('NotFound')
LD1.append(0)
SD1.append('')
AD1.append('')
CircD1.append(False)

##Some Summary Measures
TotalReadPairs1=0
LocatableReadPairs1=0
DoubleLocatableReadPairs1=0
DefinitiveDoubleLocatableReadPairs1=0
PA1=[]
for i in range(27):
    PA1.append( [ [0]*(1 + len(SD1[x])//SeparationGranularity1 ) for x in NameARange1] )

if KmerBinCoverageUnique1 or KmerBinCoverageAll1:
    UniqueKNumber1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    UniqueKObservedS1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    UniqueKObservedA1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    UniqueKCountS1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    UniqueKCountA1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    UniqueKCountS2 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    UniqueKCountA2 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    PA1.extend([UniqueKNumber1, UniqueKObservedS1, UniqueKObservedA1, UniqueKCountS1, UniqueKCountA1, UniqueKCountS2, UniqueKCountA2])
if KmerBinCoverageAll1:
    LocalKNumber1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    LocalKMultiplicity1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    LocalKObservedS1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    LocalKObservedA1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    LocalKCountS1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    LocalKCountA1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    LocalKCountS2 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    LocalKCountA2 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    PA1.extend([LocalKNumber1, LocalKMultiplicity1, LocalKObservedS1, LocalKObservedA1, LocalKCountS1, LocalKCountA1, LocalKCountS2, LocalKCountA2])
    ChromosomalKNumber1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    ChromosomalKMultiplicity1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    ChromosomalKObservedS1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    ChromosomalKObservedA1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    ChromosomalKCountS1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    ChromosomalKCountA1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    ChromosomalKCountS2 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    ChromosomalKCountA2 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    PA1.extend([ChromosomalKNumber1, ChromosomalKMultiplicity1, ChromosomalKObservedS1, ChromosomalKObservedA1, ChromosomalKCountS1, ChromosomalKCountA1, ChromosomalKCountS2, ChromosomalKCountA2])
    DispersedKNumber1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    DispersedKMultiplicity1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    DispersedKObservedS1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    DispersedKObservedA1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    DispersedKCountS1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    DispersedKCountA1 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    DispersedKCountS2 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    DispersedKCountA2 = [ np.zeros(1+len(SD1[x])//SeparationGranularity1, dtype=np.uint64) for x in NameARange1]
    PA1.extend([DispersedKNumber1, DispersedKMultiplicity1, DispersedKObservedS1, DispersedKObservedA1, DispersedKCountS1, DispersedKCountA1, DispersedKCountS2, DispersedKCountA2])

SeparationArray1=[[0]*(Tn5DupMax1+ReadSeparationMax1+1) for x in NameARange1]
MultiplicityD1=[{} for x in range(NumSeq1+1)]
D0={}  ## Archives unique combinations of position and chromosome

## Id bits for D1 dictionary.  Setting each bit indicates that the relevant read type has already been encountered
##bWellPositionedBoth=2
##bWellPositionedR1=4
##bWellPositionedR2=8
##bUnique=16
##bFocal=32
##bChromosomal=64
##bDispersed=128

D1={}

RecurrenceD0=[{} for x in range(NumSeq1+1)]
RecurrenceD1=[{} for x in range(NumSeq1+1)]
CoverageSEDType1 = np.uint8
if CoverageSEMax1>255:
    CoverageSEDType1 = np.uint16
if CoverageSEMax1>65535:
    CoverageSEDType1 = np.uint32
if CoverageSEMax1>2**32-1:
    CoverageSEDType1 = np.uint64
if FullCoverage1:
    CoverageK1=np.zeros(ulen1+1,dtype=CoverageDType1)
    if ReportStarts1:
        CoverageS1=np.zeros(ulen1+1,dtype=CoverageSEDType1)
    if ReportEnds1:
        CoverageE1=np.zeros(ulen1+1,dtype=CoverageSEDType1)

def ExplicitIndex1(Featurenumbers,Starts,Ends):
    fn1=max([max(u) for u in Featurenumbers if len(u)>0]) ## largest feature number
    fDtype1=np.uint64
    if fn1<255:
        fDtype=np.uint8
    elif fn1<65535:
        fDtype=np.uint16
    elif fn1<4**16-1:
        fDtype=np.uint32
    eI=[np.zeros(LD1[i]+klen1,dtype=fDtype) for i in xrange(NumSeq1)]
    for c in range(len(Featurenumbers)):
        for f,s,e in izip(Featurenumbers[c],Starts[c],Ends[c]):
            s=max(1,s)
            e=min(len(eI[c]),e)
            eI[c][s-1:e]=f+1  ## the entry in the array is actually featurenumber+1 so that zero can be "no feature"
    return eI
def FastFindB1(aV,N1):
    lastL=FastIndex1[N1>>BitsSorted1]    ## i=lastL-1 is the first position in aV that might hold a matching k-mer
    lastH=FastIndex1[(N1>>BitsSorted1)+1]    ## i=lastH is the first position in aV that definitely does not hold a matching k-mer
    mk1 = ( lastL!=lastH )
    lastLa=lastL[mk1]
    lastHa=lastH[mk1]
    N2=N1[mk1] & BitMask1
    while np.any(lastLa+1<lastHa):
        lastMa=(lastLa+lastHa) >> 1
        vlastMa=aV[lastMa]
        lastLa[N2>=vlastMa]=lastMa[N2>=vlastMa]
        lastHa[N2<vlastMa]=lastMa[N2<vlastMa]
    lastLa[aV[lastLa]!=N2]=len(aV)
    lastL=np.full(len(N1),ulen1,dtype=dtype1)
    lastL[mk1]=lastLa    
    return lastL
    
def FastFindS1(aV,N1):
    lastL=FastIndex1[N1>>BitsSorted1]
    lastH=FastIndex1[(N1>>BitsSorted1)+1]
    N2=N1 & BitMask1
    if lastL==lastH:
        return -1
    if lastL+1==lastH:
        if aV[lastL]==N2:
            return lastL
        else:
            return len(aV)
    while lastL<lastH:
        i=(lastL+lastH) >> 1
        if aV[i]<N2:
            if lastL==i:
                return len(aV)
            lastL=i
        elif aV[i]>N2:
            lastH=i
        else:
            return i
    return len(aV)        

DFAM1Count1=False
if DFAM1:
    LogNote1("DfamFile Specified... Making an explicit list of known repeats for categorization",LogFile1)

    DFAM1Count1=True
    LogOpeningFile1(DFAM1)
    if DFAM1.endswith('.gz'):
        if version.startswith('2.'):
            DFAMFile1=gzip.open(DFAM1,mode='r')
        else:
            DFAMFile1=gzip.open(DFAM1,mode='rt')           
    else:
        if version.startswith('2.'):
            DFAMFile1=open(DFAM1,mode='rU')
        else:
            DFAMFile1=open(DFAM1,mode='r')
    DfamS1=[[] for c in NameARange1]  ## feature start list for each chromosome
    DfamE1=[[] for c in NameARange1]  ## feature end list
    DfamI1=[[] for c in NameARange1]  ## feature index

    DfamN1={}   ## keys are feature names, values are feature index value
    DfamX1=[]   ## Values are feature names for each index
    DfamD1=[]   ## Values are narrative feature description for each index
    DfamW1=[]   ## Values are DFAM-WebID for each index
    DfamW2=[]   ## Values are DFAM-Web  Link for each index
    DfamL1=[]   ## Values are length of HMM model in DFAm database for each repeat
    DfamT1=[]   ## Values are total number of identified instances for each repeat
    DfamB1=[]   ## Values are total total number of alignable bases for each repeat (summed over all instances)
    
    DfamC1=[]  ## count of first-kmer hits to each indexed repeat type
    SortList1=[]
    for L0 in DFAMFile1:
        L1=L0.strip().split('\t')
        if len(L1)<5 or L1[0].startswith('#'):
            continue
        ch1=L1[0]
        ev1=float(L1[4]) ## evalue
        if ch1 in NameD1 and ev1<=MinRepeatEvalue1:
            cn1=NameD1[ch1]
            nam1=L1[2]
            if not(nam1 in DfamN1):
                DfamN1[nam1]=len(DfamX1)
                DfamL1.append(int(L1[8]))
                DfamD1.append(L1[16])
                DfamW1.append(L1[1])
                DfamW2.append('=HYPERLINK("http://dfam.org/entry/'+L1[1]+'")')
                DfamX1.append(nam1)
                DfamT1.append(0)
                DfamC1.append(0)
                DfamB1.append(0)
            rnum1=DfamN1[nam1]
            DfamT1[rnum1]+=1
            DfamB1[rnum1]+=abs(int(L1[7])-int(L1[6]))
            sta1=min(int(L1[10]),int(L1[11]))
            end1=max(int(L1[10]),int(L1[11]))
            SortList1.append((ev1,cn1,sta1,end1,rnum1))
    DFAMFile1.close()
    SortList1=sorted(SortList1,reverse=True)
    fnum1=0
    for ev1,cn1,sta1,end1,rnum1 in SortList1:
        DfamS1[cn1].append(sta1)
        DfamE1[cn1].append(end1)
        DfamI1[cn1].append(rnum1)
    for c in range(len(DfamS1)):
        DfamS1[c]=np.asarray(DfamS1[c],dtype=np.uint32)
        DfamE1[c]=np.asarray(DfamE1[c],dtype=np.uint32)
        DfamI1[c]=np.asarray(DfamI1[c],dtype=np.uint32)
    DfamQ1=ExplicitIndex1(DfamI1,DfamS1,DfamE1)
    LogNote1("Done making DFam index",LogFile1)
    
def FastFindP1(c,p):
    lastL=0
    lastH=len(Dfs1[c])
    while lastL<lastH-1:
        i=(lastL+lastH)//2
        if Dfe1[c][i]<p:
            lastL=i
        elif Dfs1[c][i]>p:
            lastH=i
        else:
            return i
    return -1
def FastFindR1(c,p):
    lastL=0
    lastH=len(DfamS1[c])
    while lastL<lastH-1:
        i=(lastL+lastH)//2
        if DfamE1[c][i]<p:
            lastL=i
        elif DfamS1[c][i]>p:
            lastH=i
        else:
            return DfamI1[c][i]
    return -1
alphalow1='abcdefghijklmnopqrstuvwxyz\n\t\n '

def Wormbaselink1(gn0):  ## link to a gene
    if '.' in gn0:
        gn1='.'.join(gn0.split('.')[:2])
    else:
        gn1=gn0
    gn2=gn1.rstrip(alphalow1)
    if gn2 in GeneNameD1:
        gn3=GeneNameD1[gn2]
    else:
        gn3=gn0
    return '=HYPERLINK("http://wormbase.org/db/get?name='+gn3+';class=Gene")'
    
if FindStructuralAnomalies1:
    F3=open("JunctionEvents_"+OutFileNameBase+'.tdv',mode='w')
    LogOpeningFile1(F3)
F5=open("IncidenceSummary_"+OutFileNameBase+'.tdv',mode='w')
LogOpeningFile1(F5)
F6=open("BinByBinEventSummary_"+OutFileNameBase+'.tdv',mode='w')
LogOpeningFile1(F6)
    
refV1=e4[:klen1][::-1]
maskV1=e4[klen1-1]-1
Numbase2={}
filelinenum=0
for c in 'ACGTN':
    Numbase2[c]=Numbase1[c]*(4**(klen1-1))
Dashes1000='-'*1000

## Corrected mySeq routine 06_20_17
def mySeq(Seq,ASeq,Cir,sS,sO,sD):
    ## sequence array, circular (true/false), start (one-based), offset, distance (how long of a sequence to return)
    ## start position (one based, negative for antisense)
    ## offset (constant to add to start position)
    ## distance (number of bases) to obtain
       
    if sS>0:
        lS1=len(Seq) 
        s11=sS+sO-1
        s12=sS+sO+sD-1
        S=Seq
    else:
        lS1=len(ASeq) 
        s11=lS1+sS+sO-Cir*(klen1-1)
        s12=lS1+sS+sO+sD-Cir*(klen1-1)
        S=ASeq
    if lS1==0:
        return ''
    if s11>=lS1:
        if Cir:
            s11 = s11 % lS1
            s12 = s11+sD-1
        else:
            return ''
    elif s12<=0:
        if Cir:
            s12 = s12 % lS1
            s11 = s12-sD
        else:
            return ''
    s11M=max(0,s11)
    s12M=min(lS1,s12)
    s11D=s11M-s11
    s12D=s12-s12M
    rL1=S[s11M:s12M]
    if s11D>0:
        if Cir:
            while s11D>=lS1:
                rL1=S+rL1
                s11D-=lS1
            if s11D:
                rL1=S[-s11D:]+rL1
        else:
            rL1='-'*int(s11D)+rL1
    if s12D>0:
        if Cir:
            while s12D>=lS1:
                rL1+=S
                s12D-=lS1
            if s12D:
                rL1+=S[:s12D]
        else:
            rL1+='-'*s12D
    return rL1

def UCSClink1(chrU,StartU,EndU):  ## link to a region
    if type(chrU)==int:
        chrU=NameA1[chrU]
    if chrU.startswith(chrbase1):
        return '=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?db='+UCSCAssembly1+'&position='+chrU+'%3A'+str(StartU)+'%2D'+str(EndU)+'&simpleRepeat=pack")'
    else:
        return 'No_Link'
## Candidate structural variation events are listed by
## -- a three character code,
## -- four characters to indicate the relevant read directions,
## -- and an indicator integer
## Example: 'LDJ_1s1s_411': A junction in R1, in which both ends of the junction are in the sense strand and located 411 bases apart 
## First character= Scale of event
##   'L': Local (e.g., Scale <100kb),
##   'C': Chromosomal (within chromosome, e.g., Scale>=100kb),
##   'I': Interchromosomal,
##   'B': Between Samples
##
## Second character= Nature of event
##   'D': Deletion,
##   'I': Insertion,
##   'X': Inversion Breakpoint,
##   'T': Translocation,
##   'C': Circularization
##
## Third Character= Nature of evidence
##   'J': Junction (an actual junction in one of the reads).   Indicator number is the overlap between the two reads at the junction 
##   'D': Discordant Pair (pair of reads that face away from rather than toward each other),  Indicator number is the start-to-start distance of the reads
##   'E': Extension (too long a distance between paired reads).  Indicator number is the start-to-start distance between reads
##   'C': Contraction (too short a distance between paired reads).  Indicator number is the overlap between reads (9 bases for a perfect nextera junction)  
##   'U': Unconnected pair (reads on different chromosomes or from different input DNAs).  Indicator number is the overlap at ends of the two reads   

Events0=['Circle_Discordant_1s2a', 'Circle_Discordant_1a2s', 'Circle_Contraction_1s2a', 'Circle_Contraction_1a2s',
         'Deletion_Contraction_1s2a', 'Deletion_Contraction_1a2s',
         'Circle_Junction_1s1s', 'Circle_Junction_2s2s', 'Circle_Junction_1a1a', 'Circle_Junction_2a2a',
         'Deletion_Junction_1s1s', 'Deletion_Junction_2s2s', 'Deletion_Junction_1a1a', 'Deletion_Junction_2a2a',
         'Insertion_Junction_Candidate_1s1s', 'Insertion_Junction_Candidate_2s2s', 'Insertion_Junction_Candidate_1a1a', 'Insertion_Junction_Candidate_2a2a']
if OtherSVCandidates1:
    Events0.extend(['Inversion_Candidate_1s2s', 'Inversion_Candidate_1a2a',
                   'Inversion_Candidate_1s1a', 'Inversion_Candidate_2s2a',
                   'Inversion_Candidate_1a1s', 'Inversion_Candidate_2a2s',
                   'Translocation_Candidate_1s2s', 'Translocation_Candidate_1s2a',
                   'Translocation_Candidate_1a2s', 'Translocation_Candidate_1a2a',
                   'Translocation_Candidate_1s1s', 'Translocation_Candidate_1s1a',
                   'Translocation_Candidate_1a1s', 'Translocation_Candidate_1a1a',
                   'Translocation_Candidate_2s2s', 'Translocation_Candidate_2s2a',
                   'Translocation_Candidate_2a2s', 'Translocation_Candidate_2a2a'])
Events1=[False] * len(Events0)
def strplus1(num1,lim1):
    if num1==lim1:
        return str(num1)+'+'
    else:
        return str(num1)
def CalculateFullCoverage1(Milestone1):
    FullCount1 = False
    if Milestone1.lower().endswith('final'):
        FullCount1 = True
    LogNote1("Starting Coverage Calculations for "+Milestone1,LogFile1)

    CoverageUnique1=[]
    BothStrandsUnique1=[]
    EitherStrandUnique1=[]

    CoverageLocalRpt1=[]
    BothStrandsLocal1=[]
    EitherStrandLocal1=[]

    CoverageChromosomalRpt1=[]
    BothStrandsChromosomal1=[]
    EitherStrandChromosomal1=[]

    CoverageDispersedRpt1=[]
    BothStrandsDispersed1=[]
    EitherStrandDispersed1=[]
    cknz1 = np.nonzero(CoverageK1)[0]
    for z in range(NumSeq1):
        cknzChr = cknz1[Sort_I1dedup[cknz1]==z]
        cknzChrP = cknzChr[Sort_P1dedup[cknzChr]>0]
        cknzChrM = cknzChr[Sort_P1dedup[cknzChr]<0]
        q1=cknzChrP[Multiplicities1[cknzChrP]==1]
        q2=cknzChrM[Multiplicities1[cknzChrM]==1]
        qx1= (Sort_P1dedup[q1]-1) % LD1[z]
        qx2= (-Sort_P1dedup[q2]-klen1) % LD1[z]
        CoverageZ1z=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1z[qx1]=np.minimum(CoverageK1[q1],MaxCoverHist1)
        CoverageZ1z[qx2]+=np.minimum(CoverageK1[q2],MaxCoverHist1)    
        CoverageUnique1.append(np.zeros(MaxCoverHist1+1,dtype=np.uint64))
        np.add.at(CoverageUnique1[-1],np.minimum(CoverageZ1z,MaxCoverHist1),1)
        CoverageZ1x=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1y=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1x[qx1]= (CoverageK1[q1]!=0)
        CoverageZ1y[qx2]= (CoverageK1[q2]!=0)   
        BothStrandsUnique1.append(np.count_nonzero(CoverageZ1x & CoverageZ1y))
        EitherStrandUnique1.append(np.sum(CoverageUnique1[-1][1:]))
        CoverageUnique1[-1][0] = TotalUnique1[z+1]-EitherStrandUnique1[-1]  ##11/22/20 Fixes bug in reporting numbers of uncovered k-mers
        if FullCount1 and (KmerBinCoverageUnique1 or KmerBinCoverageAll1):
            qxB1 = qx1//SeparationGranularity1
            qxB2 = qx2//SeparationGranularity1
            allchr = np.where(Sort_I1dedup==z)[0]
            allchrS = allchr[Sort_P1dedup[allchr]>0]
            allchrA = allchr[Sort_P1dedup[allchr]<0]
            allq1 = allchrS[Multiplicities1[allchrS]==1]
            allq2 = allchrA[Multiplicities1[allchrA]==1]
            allqx1 = (Sort_P1dedup[allq1]-1) % LD1[z]
            allqx2 = (-Sort_P1dedup[allq2]-klen1) % LD1[z]
            allqxB1 = allqx1//SeparationGranularity1
            allqxB2 = allqx2//SeparationGranularity1
            np.add.at(UniqueKNumber1[z],allqxB1,1)
            np.add.at(UniqueKObservedS1[z],qxB1,CoverageK1[q1]!=0)
            np.add.at(UniqueKObservedA1[z],qxB2,CoverageK1[q2]!=0)
            np.add.at(UniqueKCountS1[z],qxB1,CoverageK1[q1])
            np.add.at(UniqueKCountA1[z],qxB2,CoverageK1[q2])
            np.add.at(UniqueKCountS2[z],qxB1,CoverageK1[q1]**2)
            np.add.at(UniqueKCountA2[z],qxB2,CoverageK1[q2]**2)
            
        q1=cknzChrP[(Multiplicities1[cknzChrP] & Local1) > 0]
        q2=cknzChrM[(Multiplicities1[cknzChrM] & Local1) > 0]
        qx1= (Sort_P1dedup[q1]-1) % LD1[z]
        qx2= (-Sort_P1dedup[q2]-klen1) % LD1[z]
        CoverageZ1z=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1z[qx1]=np.minimum(CoverageK1[q1],MaxCoverHist1)
        CoverageZ1z[qx2]+=np.minimum(CoverageK1[q2],MaxCoverHist1)    
        CoverageLocalRpt1.append(np.zeros(MaxCoverHist1+1,dtype=np.uint64))
        np.add.at(CoverageLocalRpt1[-1],np.minimum(CoverageZ1z,MaxCoverHist1),1)
        CoverageZ1x=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1y=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1x[qx1]= (CoverageK1[q1]!=0)
        CoverageZ1y[qx2]= (CoverageK1[q2]!=0)   
        BothStrandsLocal1.append(np.count_nonzero(CoverageZ1x & CoverageZ1y))
        EitherStrandLocal1.append(np.sum(CoverageLocalRpt1[-1][1:]))
        CoverageLocalRpt1[-1][0] = TotalLocal1[z+1]-EitherStrandLocal1[-1]  ##11/22/20 Fixes bug in reporting numbers of uncovered k-mers
        if FullCount1 and KmerBinCoverageAll1:
            qxB1 = qx1//SeparationGranularity1
            qxB2 = qx2//SeparationGranularity1
            allq1 = allchrS[(Multiplicities1[allchrS] & Local1) > 0]
            allq2 = allchrA[(Multiplicities1[allchrA] & Local1) > 0]
            allqx1 = (Sort_P1dedup[allq1]-1) % LD1[z]
            allqx2 = (-Sort_P1dedup[allq2]-klen1) % LD1[z]
            allqxB1 = allqx1//SeparationGranularity1
            allqxB2 = allqx2//SeparationGranularity1
            np.add.at(LocalKNumber1[z],allqxB1,1)
            np.add.at(LocalKMultiplicity1[z],allqxB1,Multiplicities1[allq1] & Mult1)
            np.add.at(LocalKObservedS1[z],qxB1,CoverageK1[q1]!=0)
            np.add.at(LocalKObservedA1[z],qxB2,CoverageK1[q2]!=0)
            np.add.at(LocalKCountS1[z],qxB1,CoverageK1[q1])
            np.add.at(LocalKCountA1[z],qxB2,CoverageK1[q2])
            np.add.at(LocalKCountS2[z],qxB1,CoverageK1[q1]**2)
            np.add.at(LocalKCountA2[z],qxB2,CoverageK1[q2]**2)

        q1=cknzChrP[Multiplicities1[cknzChrP] > Chromosomal1]
        q2=cknzChrM[Multiplicities1[cknzChrM] > Chromosomal1]
        qx1= (Sort_P1dedup[q1]-1) % LD1[z]
        qx2= (-Sort_P1dedup[q2]-klen1) % LD1[z]
        CoverageZ1z=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1z[qx1]=np.minimum(CoverageK1[q1],MaxCoverHist1)
        CoverageZ1z[qx2]+=np.minimum(CoverageK1[q2],MaxCoverHist1)    
        CoverageChromosomalRpt1.append(np.zeros(MaxCoverHist1+1,dtype=np.uint64))
        np.add.at(CoverageChromosomalRpt1[-1],np.minimum(CoverageZ1z,MaxCoverHist1),1)
        CoverageZ1x=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1y=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1x[qx1]= (CoverageK1[q1]!=0)
        CoverageZ1y[qx2]= (CoverageK1[q2]!=0)   
        BothStrandsChromosomal1.append(np.count_nonzero(CoverageZ1x & CoverageZ1y))
        EitherStrandChromosomal1.append(np.sum(CoverageChromosomalRpt1[-1][1:]))
        CoverageChromosomalRpt1[-1][0] = TotalChromosomal1[z+1]-EitherStrandChromosomal1[-1]  ##11/22/20 Fixes bug in reporting numbers of uncovered k-mers
        if FullCount1 and KmerBinCoverageAll1:
            qxB1 = qx1//SeparationGranularity1
            qxB2 = qx2//SeparationGranularity1
            allq1 = allchrS[(Multiplicities1[allchrS] & Chromosomal1) > 0]
            allq2 = allchrA[(Multiplicities1[allchrA] & Chromosomal1) > 0]
            allqx1 = (Sort_P1dedup[allq1]-1) % LD1[z]
            allqx2 = (-Sort_P1dedup[allq2]-klen1) % LD1[z]
            allqxB1 = allqx1//SeparationGranularity1
            allqxB2 = allqx2//SeparationGranularity1
            np.add.at(ChromosomalKNumber1[z],allqxB1,1)
            np.add.at(ChromosomalKMultiplicity1[z],allqxB1,Multiplicities1[allq1] & Mult1)
            np.add.at(ChromosomalKObservedS1[z],qxB1,CoverageK1[q1]!=0)
            np.add.at(ChromosomalKObservedA1[z],qxB2,CoverageK1[q2]!=0)
            np.add.at(ChromosomalKCountS1[z],qxB1,CoverageK1[q1])
            np.add.at(ChromosomalKCountA1[z],qxB2,CoverageK1[q2])
            np.add.at(ChromosomalKCountS2[z],qxB1,CoverageK1[q1]**2)
            np.add.at(ChromosomalKCountA2[z],qxB2,CoverageK1[q2]**2)

        q1=cknzChrP[(Multiplicities1[cknzChrP] > 1) & (Multiplicities1[cknzChrP] < Local1)]
        q2=cknzChrM[(Multiplicities1[cknzChrM] > 1) & (Multiplicities1[cknzChrM] < Local1)]
        qx1 = (Sort_P1dedup[q1]-1) % LD1[z]
        qx2 = (-Sort_P1dedup[q2]-klen1) % LD1[z]
        CoverageZ1z=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1z[qx1]=np.minimum(CoverageK1[q1],MaxCoverHist1)
        CoverageZ1z[qx2]+=np.minimum(CoverageK1[q2],MaxCoverHist1)  ## This line is causing an unexpected bug with PyPy on Linux.      
        CoverageDispersedRpt1.append(np.zeros(MaxCoverHist1+1,dtype=np.uint64))
        np.add.at(CoverageDispersedRpt1[-1],np.minimum(CoverageZ1z,MaxCoverHist1),1)
        CoverageZ1x=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1y=np.zeros(LD1[z],dtype=np.uint8)
        CoverageZ1x[qx1]= (CoverageK1[q1]!=0)
        CoverageZ1y[qx2]= (CoverageK1[q2]!=0)   
        BothStrandsDispersed1.append(np.count_nonzero(CoverageZ1x & CoverageZ1y))
        EitherStrandDispersed1.append(np.sum(CoverageDispersedRpt1[-1][1:]))
        CoverageDispersedRpt1[-1][0] = TotalDispersed1[z+1]-EitherStrandDispersed1[-1]  ##11/22/20 Fixes bug in reporting numbers of uncovered k-mers
        if FullCount1 and KmerBinCoverageAll1:
            qxB1 = qx1//SeparationGranularity1
            qxB2 = qx2//SeparationGranularity1
            allq1 = allchrS[(Multiplicities1[allchrS] > 1) & (Multiplicities1[allchrS] < Local1)]
            allq2 = allchrA[(Multiplicities1[allchrA] > 1) & (Multiplicities1[allchrA] < Local1)]
            allqx1 = (Sort_P1dedup[allq1]-1) % LD1[z]
            allqx2 = (-Sort_P1dedup[allq2]-klen1) % LD1[z]
            allqxB1 = allqx1//SeparationGranularity1
            allqxB2 = allqx2//SeparationGranularity1
            np.add.at(DispersedKNumber1[z],allqxB1,1)
            np.add.at(DispersedKMultiplicity1[z],allqxB1,Multiplicities1[allq1] & Mult1)
            np.add.at(DispersedKObservedS1[z],qxB1,CoverageK1[q1]!=0)
            np.add.at(DispersedKObservedA1[z],qxB2,CoverageK1[q2]!=0)
            np.add.at(DispersedKCountS1[z],qxB1,CoverageK1[q1])
            np.add.at(DispersedKCountA1[z],qxB2,CoverageK1[q2])
            np.add.at(DispersedKCountS2[z],qxB1,CoverageK1[q1]**2)
            np.add.at(DispersedKCountA2[z],qxB2,CoverageK1[q2]**2)

    CoverageUnique1=prependtotal2D(CoverageUnique1)
    BothStrandsUnique1=prependtotal(BothStrandsUnique1)
    EitherStrandUnique1=prependtotal(EitherStrandUnique1)
    F5.write(F5Header1)
    F5.write('CumulativeCoverage_'+Milestone1+'\tTotalUnique\t'+'\t'.join(map(str,TotalUnique1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tCountEitherStrandUnique\t'+'\t'.join(map(str,EitherStrandUnique1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tUnique_'+str(klen1)+'-mers_PercentCoveredOnEitherStrand\t'+'\t'.join(per03(means1(TotalUnique1,EitherStrandUnique1)))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tCountBothStrandsUnique\t'+'\t'.join(map(str,BothStrandsUnique1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tUnique_'+str(klen1)+'-mers_PercentCoveredOnBothStrands\t'+'\t'.join(per03(means1(TotalUnique1,BothStrandsUnique1)))+'\t0\n')
    if Milestone1.lower().endswith('final'):
        for i in xrange(len((CoverageUnique1[0]))):
            F5.write('CoverageUnique\t'+strplus1(i,MaxCoverHist1)+'\t'+'\t'.join(map(str,[FLi[i] for FLi in CoverageUnique1]))+'\t0\n')  

    CoverageLocalRpt1=prependtotal2D(CoverageLocalRpt1)
    BothStrandsLocal1=prependtotal(BothStrandsLocal1)
    EitherStrandLocal1=prependtotal(EitherStrandLocal1)
    F5.write(F5Header1)        
    F5.write('CumulativeCoverage_'+Milestone1+'\tTotalLocalRpt\t'+'\t'.join(map(str,TotalLocal1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tCountEitherStrandLocalRpt\t'+'\t'.join(map(str,EitherStrandLocal1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tLocalRepeat_'+str(klen1)+'-mers_PercentCoveredOnEitherStrand\t'+'\t'.join(per03(means1(TotalLocal1,EitherStrandLocal1)))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tCountBothStrandsLocalRpt\t'+'\t'.join(map(str,BothStrandsLocal1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tLocalRepeat_'+str(klen1)+'-mers_PercentCoveredOnBothStrands\t'+'\t'.join(per03(means1(TotalLocal1,BothStrandsLocal1)))+'\t0\n')
    if Milestone1.lower().endswith('final'):
        for i in xrange(len((CoverageLocalRpt1[0]))):
            F5.write('CoverageLocalRpt\t'+strplus1(i,MaxCoverHist1)+'\t'+'\t'.join(map(str,[FLi[i] for FLi in CoverageLocalRpt1]))+'\t0\n')  

    CoverageChromosomalRpt1=prependtotal2D(CoverageChromosomalRpt1)
    BothStrandsChromosomal1=prependtotal(BothStrandsChromosomal1)
    EitherStrandChromosomal1=prependtotal(EitherStrandChromosomal1)
    F5.write(F5Header1)        
    F5.write('CumulativeCoverage_'+Milestone1+'\tTotalChromosomalRpt\t'+'\t'.join(map(str,TotalChromosomal1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tCountEitherStrandChromosomalRpt\t'+'\t'.join(map(str,EitherStrandChromosomal1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tChromosomeWideRepeat_'+str(klen1)+'-mers_PercentCoveredOnBothStrands\t'+'\t'.join(per03(means1(TotalChromosomal1,EitherStrandChromosomal1)))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tCountBothStrandsChromosomalRpt\t'+'\t'.join(map(str,BothStrandsChromosomal1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tChromosomeWideRepeat_'+str(klen1)+'-mers_PercentCoveredOnEitherStrand\t'+'\t'.join(per03(means1(TotalChromosomal1,BothStrandsChromosomal1)))+'\t0\n')
    if Milestone1.lower().endswith('final'):
        for i in xrange(len((CoverageChromosomalRpt1[0]))):
            F5.write('CoverageChromosomalRpt\t'+strplus1(i,MaxCoverHist1)+'\t'+'\t'.join(map(str,[FLi[i] for FLi in CoverageChromosomalRpt1]))+'\t0\n')  

    CoverageDispersedRpt1=prependtotal2D(CoverageDispersedRpt1)
    BothStrandsDispersed1=prependtotal(BothStrandsDispersed1)
    EitherStrandDispersed1=prependtotal(EitherStrandDispersed1)
    F5.write(F5Header1)        
    F5.write('CumulativeCoverage_'+Milestone1+'\tTotalDispersedRpt\t'+'\t'.join(map(str,TotalDispersed1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tCountEitherStrandDispersedRpt\t'+'\t'.join(map(str,EitherStrandDispersed1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tDispersedRepeat_'+str(klen1)+'-mers_PercentCoveredOnBothStrands\t'+'\t'.join(per03(means1(TotalDispersed1,EitherStrandDispersed1)))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tCountBothStrandsDispersedRpt\t'+'\t'.join(map(str,BothStrandsDispersed1))+'\t0\n')
    F5.write('CumulativeCoverage_'+Milestone1+'\tDispersedRepeat_'+str(klen1)+'-mers_PercentCoveredOnEitherStrand\t'+'\t'.join(per03(means1(TotalDispersed1,BothStrandsDispersed1)))+'\t0\n')
    if Milestone1.lower().endswith('final'):
        for i in xrange(len((CoverageDispersedRpt1[0]))):
            F5.write('CoverageDispersedRpt\t'+strplus1(i,MaxCoverHist1)+'\t'+'\t'.join(map(str,[FLi[i] for FLi in CoverageDispersedRpt1]))+'\t0\n')  
        
    LogNote1("Finished Coverage Calculations for "+Milestone1,LogFile1)

FiL1=[]
FiD1=[]
FiN1=[]
FiLR1=[]
FiLR2=[]
FiLN1=[]
if BaseByBase1:
    BfN1 = []
    ExtensionD1 = {}  ## Keys are (chromosome#,position), values are dictionaries with key being the extension nature and value being number of occurences
    if not(BaseByBaseRegionsFiles1) and not (BaseByBasePositions1):
        for c99 in NameARange1[:-1]:
            BfN1.append(NameA1[c99])  ## in this case BfN1 will simply be a copy of NameA1
            for Bfi1 in range(LD1[c99]):
                BaseByBaseD1[(c99,Bfi1)] = (c99,Bfi1)
    if BaseByBasePositions1:
        BfN1 = NameA1
        for pp0 in BaseByBasePositions1:
            pp1,pp2 = pp0.split(':')
            pp3=int(pp2.split('-')[0])
            pp4=int(pp2.split('-')[-1])
            ch = pp1
            if ch.lower().startswith('m') and not(ch in NameD1) and not('chr'+ch in NameD1):  ## very cumbersome temporary code trying to deal with different ways of naming MtDNA
                if ('M' in NameD1) or ('chrM' in NameD1):
                    ch = 'M'
            if not ch in NameD1:
                ch = 'chr'+ch
            if ch in NameD1:
                for pp5 in xrange(pp3-1,pp4):
                    BaseByBaseD1[(NameD1[ch],pp5)] = (NameD1[ch],pp5)
    if BaseByBaseRegionsFiles1:
         ## Mnemonic for each feature
        for CurRegionFile1 in BaseByBaseRegionsFiles1:
            if CurRegionFile1.endswith('gz'):
                if version.startswith('2.'):
                    bbFF1 = gzip.open(CurRegionFile1,mode='r')
                else:
                    bbFF1 = gzip.open(CurRegionFile1,mode='rt')
            else:
                if version.startswith('2.'):
                    bbFF1 = open(CurRegionFile1,mode='rU')
                else:
                    bbFF1 = open(CurRegionFile1,mode='r')
            LogOpeningFile1(CurRegionFile1)
            for nf0,L0 in enumerate(bbFF1):
                L1=L0.strip().split('\t')
                if len(L1)<5 or not(L1[3].isdigit()) or not(L1[4].isdigit()):
                    continue
                ch = L1[0]
                if ch.lower().startswith('m') and not(ch in NameD1) and not('chr'+ch in NameD1):  ## very cumbersome temporary code trying to deal with different ways of naming MtDNA
                    if ('M' in NameD1) or ('chrM' in NameD1):
                        ch = 'M'
                if not ch in NameD1:
                    ch = 'chr'+ch
                BfT1 = L1[2]
                if BaseByBaseCategories1 and not(BfT1 in BaseByBaseCategories1):
                    continue
                if BaseByBaseTags1:
                    Keep1 = False
                    for ft1 in BaseByBaseTags1:
                        if ft1.lower() in L0.lower():
                            Keep1 = True
                            break
                    if not(Keep1):
                        continue                       
                if (ch in NameD1):
                    BfN1.append( L1[6]+uidfinder(L0) )
                    BfC1 = NameD1[ch]
                    BfS1 = int(L1[3])
                    BfE1 = int(L1[4])
                    for i in range(BfS1-1,BfE1):
                        if not((BfC1,i)  in BaseByBaseD1):  ## In case of overlap between features (which is the rule for ENSEMBL datasets rather than the exception), Each base belongs to the first feature that contains it
                            if L1[6]=='-':
                                BaseByBaseD1[(BfC1,i)] = (len(BfN1)-1,i-BfE1) ## Keys are chromosome number, position, values are feature number, -position
                        else:
                            BaseByBaseD1[(BfC1,i)] = (len(BfN1)-1,i-BfS1+1) ## Keys are chromosome number, position, values are feature number, position
def versioner1(v88):
    s88 = ['0']
    for c1 in v88:
        if c1.isdigit():
            s88[-1]+=c1
        else:
            if s88[-1]!='0':
                s88.append('0')
    return(map(int,s88))
fqdVersion1 = ['0']
def Moment1(a):
    s1 = sum(a)
    s2 = sum(a[i]*i for i in range(len(a)))
    if s1==0:
        s1=1
    return float(s2)/s1

for f0,d0,n0 in zip(FiL0,FiD0,FiN0):
    if os.path.isfile(f0) or f0.endswith('.fastq') or f0.endswith('.fastq.gz') or f0.endswith('.fasta') or f0.endswith('.fasta.gz'):
        FiL1.append(f0)
        FiD1.append(d0)
        FiN1.append(n0)
    else:
        LogNote1(f0+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download",LogFile1)
        LogNote1("Preparing to download sequence read set "+f0+" from NCBI",LogFile1)
        try:
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,'-version'])
        except:
            LogNote1("Searching for a version of fastq-dump that will run; if this fails, you may need to redownload the program and unzip the archive, also add FastQDumpProgram=<path to program> to command line",LogFile1)
            os.environ["PATH"]=os.getenv("PATH")+':./:/opt/local/bin:/opt/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/X11/bin:/Applications:~/Downloads'
            for fqd1 in os.getenv("PATH").split(':'):
                if os.path.isdir(fqd1):
                    for fqd2 in os.listdir(fqd1):
                        if fqd2.startswith('sratoolkit'):
                            fqd3 = os.path.join(fqd1,fqd2)
                            if os.path.isdir(fqd3):
                                fqd4 = os.path.join(fqd3,'bin','fastq-dump')
                                if os.path.isfile(fqd4):
                                    if versioner1(fqd2) > versioner1(fqdVersion1):
                                        fqdVersion1 = fqd2
                                        FastQDumpProgram1 = fqd4
                                        TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,'-version'])
        LogNote1("Trying presumed fastq-dump program file located at "+FastQDumpProgram1,LogFile1)
        if LinesToProcess1<1:
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,
                                                 '--fasta',
                                                 '0',
                                                 '--origfmt',
                                                 '--split-files',
                                                 '--gzip',
                                                 '--outdir',
                                                 ScratchPath1,
                                                 f0])
        else:
            TryFastQDump1 = subprocess.check_output(['fastq-dump',
                                                 '--fasta',
                                                 '0',
                                                 '--origfmt',
                                                 '--split-files',
                                                 '--gzip',
                                                 '-X',
                                                 str(LinesToProcess1),
                                                 '--outdir',
                                                 ScratchPath1,
                                                 f0])
        LogNote1("Result of "+f0+" NCBI Download " + TryFastQDump1,LogFile1)
        PresumptiveRead1FilePath1 = os.path.join(ScratchPath1,f0+'_1.fasta.gz')
        PresumptiveRead1FilePath2 = os.path.join(ScratchPath1,f0+'_2.fasta.gz')
        if os.path.isfile(PresumptiveRead1FilePath1):
            SRRReadPath1 = os.path.join(ScratchPath1,f0.lower()+'R1_'+TempFileUID1+'.fasta.gz')
            SRRReadPath2 = os.path.join(ScratchPath1,f0.lower()+'R2_'+TempFileUID1+'.fasta.gz')
            os.rename(PresumptiveRead1FilePath1, SRRReadPath1)
            SRRTempList1.append(SRRReadPath1)
            FiL1.append(SRRReadPath1)
            FiD1.append('R1')
            FiN1.append(n0)
            if os.path.isfile(PresumptiveRead1FilePath2):
                os.rename(PresumptiveRead1FilePath2, SRRReadPath2)
                SRRTempList1.append(SRRReadPath2)
for f0,d0,n0 in zip(FiL1,FiD1,FiN1):
    if d0=='R1':
        rDom1 = 'R1'
    elif  d0=='R2':
        rDom1 = 'R2'
    else:
        rDom1='R1'
        if 'R1' in f0 and not('R2' in f0):
            rDom1='R1'
        elif 'R2' in f0 and not('R1' in f0) and not(('SRR2' in f0) and f0.count('R2')==1):
            rDom1='R2'
        elif 'R1' in f0 and 'R2' in f0:
            if f0.rfind('R1')>f0.rfind('R2'):
                rDom1='R1'
            else:
                rDom1='R2'
    if rDom1=='R1':
        f0r1=f0[:]
        f0r2=f0[::-1].replace('1R','2R',1)[::-1]
    else:
        f0r2=f0[:]
        f0r1=f0[::-1].replace('2R','1R',1)[::-1]
    if not(f0r1 in FiLR1):
        if os.path.isfile(f0r1):
            FiLR1.append(f0r1)
        else:
            FiLR1.append('')
        if (f0r1 != f0r2) and os.path.isfile(f0r2):
            FiLR2.append(f0r2)
        else:
            FiLR2.append('')
        FiLN1.append(n0)
        if FiLR1[-1]=='' and FiLR2[-1]=='':
            LogNote1("***Warning*** Input Files "+f0r1+' and '+f0r2+' not found.  Trying to proceed with any other files that are found',LogFile1)
if FiLR1==[]:
    LogNote1("NO READ INPUT FILES (e.g., fastq's) SPECIFIED- RUNNING IN DIAGNOSTIC MODE WITH A FEW TEST SEQUENCES",LogFile1)
    FiLR1.append(test1)
    FiLR2.append(test2)
    
if any(FiLR2):
    Read2Input = True
else:
    Read2Input = False
if any(FiLR1):
    Read1Input = True
else:
    Read1Input = False

TaskHeader1 = '<!--REVA_Task_Header: '+OutFileNameBase+'-->\n'+\
    '<!--Reference_File: '+RefSeqFile1+'-->\n'
for fin1,(fi1,fi2) in enumerate(zip(FiLR1,FiLR2)):
    if fi1:
        TaskHeader1+='<!--Read_1_File_'+str(fin1)+': '+FileInfo1(fi1)+'-->\n'
    if fi2:
        TaskHeader1+='<!--Read_2_File_'+str(fin1)+': '+FileInfo1(fi2)+'-->\n'
for fin1,GFF11 in enumerate(GFF1):
    TaskHeader1+='<!--Feature_File_'+str(fin1)+': '+FileInfo1(GFF11)+'-->\n'
for fin1,GFF11 in enumerate(BaseByBaseRegionsFiles1):
    TaskHeader1+='<!--BaseByBaseRegions_File_'+str(fin1)+': '+FileInfo1(GFF11)+'-->\n'
TaskHeader1 +=  '<!--DFAM_File: '+DFAM1+'-->\n'+\
    '<!--K-mer_Length: '+str(klen1)+'-->\n'+\
    '<!--Command_Line: '+' '.join(argv)+'-->\n'+\
    '<!--PythonVersion: '+','.join(version.splitlines())+'-->\n'+\
    '<!--REVA_Version: '+FileInfo1(argv[0])+'-->\n'
TaskHeader1+='\n'
AbbrevHeader1 = ''.join(TaskHeader1.splitlines()[:-1])+'<!--REVATableHeader-->'  ##ending with ':REVATableHeader' identifies a line as a row of table headers

NameA2 = [na1+'__'+RefAbbrev1 for na1 in NameA1]
F5Header1='Attribute\tValue\tAll\t'+'\t'.join(NameA2)+'\t'+AbbrevHeader1+'\n'

LogNote1("Starting to record Coverage Overhead",LogFile1)
CoverageU1=[]
for z in range(NumSeq1):
    CoverageU1.append(np.zeros(LD1[z],dtype=np.uint8))
F5.write('<!--REVA_Whole_Run_Statstics-->\n')
F5.write(TaskHeader1+'\n')
F5.write(HeaderTranspose(F5Header1)+'\n')
F5.write(F5Header1)        
F5.write('TotalUnique__'+RefAbbrev1+'\t'+str(klen1)+'-mers\t'+'\t'.join(map(str,TotalUnique1))+'\t0\n')
F5.write('TotalLocalRepeat__'+RefAbbrev1+'\t'+str(klen1)+'-mers\t'+'\t'.join(map(str,TotalLocal1))+'\t0\n')
F5.write('TotalChromosomalRepeat__'+RefAbbrev1+'\t'+str(klen1)+'-mers\t'+'\t'.join(map(str,TotalChromosomal1))+'\t0\n')
F5.write('TotalDispersedRepeat__'+RefAbbrev1+'\t'+str(klen1)+'-mers\t'+'\t'.join(map(str,TotalDispersed1))+'\t0\n')
LogNote1("Finishing recording of Coverage Overhead",LogFile1)
F5Lengths1='TotalBases__'+RefAbbrev1+'\tBases\t'+str(sum(LD1))+'\t'+'\t'.join(map(str,[LD1[z] for z in NameARange1]))+'\n'
F5.write(F5Lengths1)


if FindStructuralAnomalies1:
    OutputHeader0='\t'.join(['EventTypes',
            'EventChr__'+RefAbbrev1,
            'EventStartLo',
            'EventStartHi',
            'EventEndLo',
            'EventEndHi',
            'EventOvrelap',
            'UCSC_Link'])
    OutputHeader1='\t'+'\t'.join(['r1',
                'r1L_Chromosome',
                'r1L_StartBase',
                'r1L_HomologyLen',
                'r1L_ReferenceSeq',
                'r1R_Chromosome',
                'r1R_StartBase',
                'r1R_HomologyLen',
                'r1R_ReferenceSeq',
                'r1_RtoL_Distance_On_Genome',
                'r1_RtoL_Homology_Overlap'])
    OutputHeader2='\t'+'\t'.join(['r2',
                'r2L_Chromosome',
                'r2L_StartBase',
                'r2L_HomologyLen',
                'r2L_ReferenceSeq',
                'r2R_Chromosome',
                'r2R_StartBase',
                'r2R_HomologyLen',
                'r2R_ReferenceSeq',
                'r2_RtoL_Distance_On_Genome',
                'r2_RtoL_Homology_Overlap'])
    OutputHeader3='\tReadNumberInFile'
    F3.write('<!--REVA_StructuralAnomalyList-->\n')
    F3.write(TaskHeader1+'\n')
    if VerboseOutput1:
        if Read1Input:                                                   
            OutputHeader0 += OutputHeader1
        if Read2Input:                                                   
            OutputHeader0 += OutputHeader2
        OutputHeader0 += OutputHeader3
    if ReportBuffers1:
        if R1Buffer5:
            OutputHeader0 += '\tR1Buffer5'
        if R1Buffer3:
            OutputHeader0 += '\tR1Buffer3'
        if R2Buffer5:
            OutputHeader0 += '\tR2Buffer5'
        if R2Buffer3:
            OutputHeader0 += '\tR2Buffer3'
    F3.write(HeaderTranspose(OutputHeader0)+'\n')
    F3.write(OutputHeader0)
    F3.write('\t'+AbbrevHeader1)
    F3.write('\n')


EitherStrandUnique1=[0]*NumSeq1
XA1=[]
exptlinenum1=0
klen2 = klen1//2
if not(klen1d):
    klen1d = 2*klen1
klen1dm2 = 2*klen1-2
LastL0,LastM0 = '',''
MaxCoverage1=np.iinfo(CoverageDType1).max
SingleReadMode1 = 0  ## 0 means both R1 and R2, 1 for just R1, 2 for just R2

if FirstKTable1 or InterimKTable1:
    J0A  = [0]*(klen1d+1); J1A  = [0]*(klen1d+1); J2A  = [0]*(klen1d+1); J3A  = [0]*(klen1d+1)
    JS0A = [0]*(klen1d+1); JS1A = [0]*(klen1d+1); JS2A = [0]*(klen1d+1); JS3A = [0]*(klen1d+1)
    JA0A = [0]*(klen1d+1); JA1A = [0]*(klen1d+1); JA2A = [0]*(klen1d+1); JA3A = [0]*(klen1d+1)        
if FirstKTable1 and InterimKTable1:
    LogNote1('FirstKTable and InterimKTable options are incompatible, running InterimKTable Only.  Run without InterimKTable to obtain global KTable',LogFile1)        
    FirstKTable1 = False

for (f1,f2,n00) in zip(FiLR1,FiLR2,FiLN1):
    if R1Only:
        f2 = ''
    if R2Only:
        f1 = ''
    if f1:
        if f1.endswith('.gz'):
            if version.startswith('2.'):
                F1=gzip.open(f1,mode='r')
            else:
                F1=gzip.open(f1,mode='rt')
        else:
            if version.startswith('2.'):
                F1=open(f1,mode='rU')
            else:
                F1=open(f1,mode='r')                
        LogOpeningFile1(F1)
    else:
        F1 = iter(())
        SingleReadMode1 = 2
        if not(R2Only):
            LogNote1('no Read1 file for '+f2+' was found, running REVA in single read mode with just read file '+f1,LogFile1)
    if f2:
        if f2.endswith('.gz'):
            if version.startswith('2.'):
                F2=gzip.open(f2,mode='r')
            else:
                F2=gzip.open(f2,mode='rt')
        else:
            if version.startswith('2.'):
                F2=open(f2,mode='rU')
            else:
                F2=open(f2,mode='r')
        LogOpeningFile1(F2)
    else:
        F2 = iter(())
        SingleReadMode1 = 1
        if not(R1Only):
            LogNote1('no Read2 file for '+f1+' was found, running REVA in single read mode with just read file '+f1,LogFile1)
    filelinenum=0
    if f1:
        Mnemonic2=os.path.basename(f1)[::-1].split('1R',1)[-1][::-1].strip('_')  # file base for a given run
    else:
        Mnemonic2=os.path.basename(f2)[::-1].split('2R',1)[-1][::-1].strip('_')
    if '#' in n00:
        Mnemonic2 += 'ie'+n00.split('#',1)[1]
    FastAFile1 = False
    LineDensity1 = 4
    HotLine1 = 2
    if f1.lower().endswith('fasta') or  f1.lower().endswith('fasta.gz') or f2.lower().endswith('fasta') or  f2.lower().endswith('fasta.gz'):
        FastAFile1 = True
        LineDensity1 = 2
        HotLine1 = 0
    LinesToProcessA=LinesToProcess1*LineDensity1
    ReportIntervalA=LineDensity1*ReportInterval1
    FullCoverageIntervalA=LineDensity1*FullCoverageInterval1
    CoverageAccumulationIntervalA=LineDensity1*CoverageAccumulationInterval1
    PreFilters1 = 0
    PassFilters1 = 0
    if SeqIndexMode1>0:
        SeqIndexD1 = {}
        for rn1,(L0,M0) in enumerate(izip_longest(F1,F2,fillvalue='')):
            if rn1%LineDensity1 == 0:
                bcL1 = GetSeqIndexFromReadID1(L0)
                bcM1 = GetSeqIndexFromReadID1(M0)
                if not((bcL1,bcM1)) in SeqIndexD1:
                    SeqIndexD1[(bcL1,bcM1)] = 0
                SeqIndexD1[(bcL1,bcM1)] +=  1
        SeqIndexMax1 = max(SeqIndexD1.values())
        for key1 in list(SeqIndexD1):
            if not SeqIndexD1[key1]==SeqIndexMax1:
                del(SeqIndexD1[key1])
        if f1:
            F1.seek(0)
        if f2:
            F2.seek(0)
    if f1.lower().endswith('fasta') or  f1.lower().endswith('fasta.gz') or f2.lower().endswith('fasta') or  f2.lower().endswith('fasta.gz'):
        F1 = chain(F1,['>'])
        F2 = chain(F2,['>'])
    for L0,M0 in izip_longest(F1,F2,fillvalue=''):
        if (filelinenum & 3 == 0 or (FastAFile1 and (filelinenum & 3 == 2))) :
            GoodSeqIndex1 = True
            if SeqIndexD1:
                bcL1 = L0.strip().split(':')[-1]
                bcM1 = M0.strip().split(':')[-1]
                if not((bcL1,bcM1) in SeqIndexD1):
                    GoodSeqIndex1 = False
        filelinenum+=1
        if filelinenum % ReportIntervalA==0:
            LogNote1('finished '+Mnemonic2+
                     ' read: '+str(filelinenum//LineDensity1)+
                     ' PreFilters: '+str(PreFilters1)+
                     ' PassFilters: '+str(PassFilters1),LogFile1)
            PreFilters1 = 0
            PassFilters1 = 0
            if InterimKTable1:
                if SingleReadMode1==0:
                    LogNote1('    Average_FirstKValues_Cumulative: StartR1=' + '{0:.2f}'.format(Moment1(J0A))+
                                                 ', EndR1='   + '{0:.2f}'.format(Moment1(J1A))+
                                                 ', StartR2=' + '{0:.2f}'.format(Moment1(J2A))+
                                                 ', EndR2='   + '{0:.2f}'.format(Moment1(J3A)),LogFile1)
                elif SingleReadMode1==1:
                    LogNote1('    Average_FirstKValues_Cumulative: StartR1=' + '{0:.2f}'.format(Moment1(J0A))+
                                                 ', EndR1='   + '{0:.2f}'.format(Moment1(J1A)),LogFile1) 
                elif SingleReadMode1==2:
                    LogNote1('    Average_FirstKValues_Cumulative: StartR2=' + '{0:.2f}'.format(Moment1(J2A))+
                                                 ', EndR2='   + '{0:.2f}'.format(Moment1(J3A)),LogFile1)
                J0A  = [0]*(klen1d+1); J1A  = [0]*(klen1d+1); J2A  = [0]*(klen1d+1); J3A  = [0]*(klen1d+1)
                JS0A = [0]*(klen1d+1); JS1A = [0]*(klen1d+1); JS2A = [0]*(klen1d+1); JS3A = [0]*(klen1d+1)
                JA0A = [0]*(klen1d+1); JA1A = [0]*(klen1d+1); JA2A = [0]*(klen1d+1); JA3A = [0]*(klen1d+1)        

        if (BriefCoverage1 or FullCoverage1) and ((filelinenum % CoverageAccumulationIntervalA==0) or (filelinenum == LinesToProcessA)):
            LogNote1("Starting Interim Coverage Calculations for read "+str(exptlinenum1),LogFile1)
            XA2=np.unique(np.asarray(XA1,dtype=FastIndexDataType1),return_counts=True)
            if FullCoverage1:
                CoverageK1[XA2[0]]=np.minimum(CoverageK1[XA2[0]]+XA2[1],MaxCoverage1)
            XA1=[]
            if BriefCoverage1:
                IA1=Sort_I1dedup[XA2[0]]
                PAx1=Sort_P1dedup[XA2[0]]
                MA1=Multiplicities1[XA2[0]]
                for iA1 in range(NumSeq1):
                    PA2=PAx1[ (IA1==iA1) & (MA1==1)]
                    PA3=np.unique( ((np.abs(2*PA2+klen1-1)-klen1-1)>>1) % LD1[iA1] )
                    EitherStrandUnique1[iA1]+=len(PA3)-np.count_nonzero(CoverageU1[iA1][PA3])
                    CoverageU1[iA1][PA3]=1
                EitherStrandUnique2=prependtotal(EitherStrandUnique1)
                LabelID1="Interim"
                F5.write('UniqueCoverageCount'+LabelID1+'\t'+str(exptlinenum1)+'\t'+'\t'.join(map(str,EitherStrandUnique2))+'\t0\n')
                M1=per03(means1(TotalUnique1,EitherStrandUnique2))
                F5.write('UniqueCoveragePercent'+LabelID1+'\t'+str(exptlinenum1)+'\t'+'\t'.join(M1)+'\t0\n')
                if M1:
                    LogNote1('After '+str(exptlinenum1)+' read pairs, unique coverage stands at '+M1[0]+'% ',LogFile1)
                else:
                    LogNote1('After '+str(exptlinenum1)+' read pairs, unique coverage stands at '+'0.000'+'% ',LogFile1)
            LogNote1("Finished Interim Coverage Calculations for read "+str(exptlinenum1),LogFile1)
        if filelinenum == LinesToProcessA:
            LogNote1('Hit pre-specified line maximum of '+str(LinesToProcess1)+' reads.  Proceeding after '+Mnemonic2,LogFile1)
            break
        if FullCoverageInterval1>0 and filelinenum % FullCoverageIntervalA==0:
            MilestoneL1=Mnemonic2+'_'+str(filelinenum//LineDensity1)+'_'+str(exptlinenum1)
            if IndexConstructionBins1>1:
                npush('Sort_V1dedup') ## release this memory before calling CalculateFullCoverage
            CalculateFullCoverage1(MilestoneL1)
            if IndexConstructionBins>1:
                npull('Sort_V1dedup')
            LogNote1('recorded '+Mnemonic2+' coverage: '+str(filelinenum//LineDensity1),LogFile1)
        if (filelinenum & 3 == 2 or (FastAFile1 and (filelinenum & 3 == 0))) and ((GoodSeqIndex1 and SeqIndexMode1<2) or (not(GoodSeqIndex1) and SeqIndexMode1==2)):    ##only pay attention to the second line of each 4 (fastq file structure)
            PreFilters1 += 1
            if BarcodeLenR1>0 and L0:
                if BarcodeSeqR1 and BarcodeRequireR1 and L0[:BarcodeLenR1]!=BarcodeSeqR1:
                    L0 = ''
                elif L0[:BarcodeLenR1]==BarcodeSeqR1 or AnyBarcodeR1:
                    L0 = L0[:BarcodeLenR1]
            if BarcodeLenR2>0 and M0:
                if BarcodeSeqR2 and BarcodeRequireR2 and M0[:BarcodeLenR2]!=BarcodeSeqR2:
                    M0 = ''
                elif M0[:BarcodeLenR2]==BarcodeSeqR2 or AnyBarcodeR2:
                    M0 = M0[:BarcodeLenR2]
            if LinkerR1:
                LR1 = L0.rfind(LinkerR1)
                if LR1>=0:
                    L0 = L0[:LR1]
                elif LinkerRequireR1:
                    L0 = ''
            if LinkerR2:
                LR2 = M0.rfind(LinkerR2)
                if LR2>=0:
                    M0 = M0[:LR2]
                elif LinkerRequireR2:
                    M0 = ''
            exptlinenum1+=1
            TotalReadPairs1+=1
            L0=L0.strip()
            M0=M0.strip()
            if any((R1Buffer5,R1Buffer3,R2Buffer5,R2Buffer3)):
                L00=L0
                M00=M0
                if R1Buffer5>0:
                    L0 = L0[R1Buffer5:]
                if R1Buffer3>0:
                    L0 = L0[:-R1Buffer3]
                if R2Buffer5>0:
                    M0 = M0[R2Buffer5:]
                if R2Buffer3>0:
                    M0 = M0[:-R2Buffer3]
            if DeleteNs1:
                L0=L0.replace('N','')
                M0=M0.replace('N','')
            lL0=len(L0)
            lM0=len(M0)
            x0=-1; x1=-1; x2=-1; x3=-1 ## ordinal positions for first uniquely matched k-mers for beginning of R1, last k-mer for end of R1, first k-mer for beginning of read1, last k-mer for end of R2
            p0=0; p1=0; p2=0; p3=0  ## genome positions of the first matched base for x0-x3 
            i0=0; i1=lL0-klen1; i2=0; i3=lM0-klen1  ## i0,1,2,3 are the offsets between the beginning of L0 and M0 and the beginning of the first matched k-mer.  zero if the first uniquely matched k-mer is at the beginning of the sequence.  i1 and i3 are the distance between the end of the last uniquely matched k-mer and the end of the sequence 
            c0=NumSeq1; c1=NumSeq1; c2=NumSeq1; c3=NumSeq1  ## corresponding chromosomal assignments
            ix0=-1; ix1=-1; ix2=-1; ix3=-1  ##ipn, iin, and ixn are equivalent values but for the first match of any kind (not just unique)
            ip0=0; ip1=0; ip2=0; ip3=0  
            ii0=0; ii1=i1; ii2=0; ii3=i3  ## i0,1,2,3 are the offsets between the beginning of L0 and M0 and the beginning of the first matched k-mer.  All are zero if the first matched k-mer is at the beginning of the sequence
            ic0=NumSeq1; ic1=NumSeq1; ic2=NumSeq1; ic3=NumSeq1  ## corresponding chromosomal assignments
            X0=[]; X2=[]  ## Lists of hits
            HomologyLength1 = 0
            if lL0<klen1 and lM0<klen1:
                continue
            if lL0>=klen1:
                Lx0=seqnum(L0).astype(dtype1)
                Ly0=np.zeros(lL0-klen1+1,dtype=dtype1)
                z0=0
                for kin1 in eList32:
                    if klen1 & kin1:                   
                        if z0==0:
                            Ly0+=Lx0[klen1-kin1:]
                        else:
                            Ly0+=Lx0[klen1-kin1-z0:-z0]*e4[z0]
                        z0+=kin1
                    if klen1>z0:
                        Lx0=Lx0[:-kin1]*e4[kin1]+Lx0[kin1:]
                X0=FastFindB1(Sort_V1dedup,Ly0)
                while i0<=i1:
                    x0=X0.item(i0)
                    if x0<ulen1:
                        if ix0==-1:
                            ix0 = x0
                            ii0 = i0
                            ip0 = Sort_P1dedup.item(ix0)
                            ic0 = Sort_I1dedup.item(ix0)
                        if Multiplicities1.item(x0)==1:
                            p0 = Sort_P1dedup.item(x0)
                            c0 = Sort_I1dedup.item(x0)
                            break
                    i0 += 1
                while i1>=ii0:
                    x1=X0.item(i1)
                    if x1<ulen1:
                        if ix1==-1:
                            ix1 = x1
                            ii1 = i1
                            ip1 = Sort_P1dedup.item(ix1)
                            ic1 = Sort_I1dedup.item(ix1)
                        if Multiplicities1.item(x1)==1:
                            p1 = Sort_P1dedup.item(x1)
                            c1 = Sort_I1dedup.item(x1)
                            break
                    i1 -= 1
            if lM0>=klen1:
                Mx0=seqnum(M0).astype(dtype1)
                My0=np.zeros(lM0-klen1+1,dtype=dtype1)
                z0=0
                for kin1 in eList32:
                    if klen1 & kin1:                   
                        if z0==0:
                            My0+=Mx0[klen1-kin1:]
                        else:
                            My0+=Mx0[klen1-kin1-z0:-z0]*e4[z0]
                        z0+=kin1
                    if klen1>z0:
                        Mx0=Mx0[:-kin1]*e4[kin1]+Mx0[kin1:]
                X2=FastFindB1(Sort_V1dedup,My0)
                while i2<=i3:
                    x2=X2.item(i2)
                    if x2<ulen1:
                        if ix2==-1:
                            ix2 = x2
                            ii2 = i2
                            ip2 = Sort_P1dedup.item(ix2)
                            ic2 = Sort_I1dedup.item(ix2)
                        if Multiplicities1.item(x2)==1:
                            p2 = Sort_P1dedup.item(x2)
                            c2 = Sort_I1dedup.item(x2)
                            break
                    i2 += 1
                while i3>=ii2:
                    x3=X2.item(i3)
                    if x3<ulen1:
                        if ix3==-1:
                            ix3 = x3
                            ii3 = i3
                            ip3 = Sort_P1dedup.item(ix3)
                            ic3 = Sort_I1dedup.item(ix3)
                        if Multiplicities1.item(x3)==1:
                            p3 = Sort_P1dedup.item(x3)
                            c3 = Sort_I1dedup.item(x3)
                            break
                    i3 -= 1 
            ## Conditions for skipping a completely nonrelated read pair
            if ip0==0 and ip2==0: continue
            ## Calculate values indicating inferred starts and fragment lengths
            if p0!=0:
                Es0=p0-i0
            else:
                Es0=ip0-ii0
                c0=ic0
            if p1!=0:
                Es1=p1-i1
            else:
                Es1=ip1-ii1
                c1=ic1
            if p2!=0:
                Es2=p2-i2
            else:
                Es2=ip2-ii2
                c2=ic2
            if p3!=0:
                Es3=p3-i3
            else:
                Es3=ip3-ii3
                c3=ic3
            if SingleReadMode1==1:
                Es2 = -(Es1+lL0-1)
                Es3 = -(Es0+lL0-1)
                c2 = c1
                c3 = c0
            elif SingleReadMode1==2:
                Es0 = -(Es2+lM0-1)
                Es1 = -(Es3+lM0-1)
                c1 = c2
                c3 = c0
            Em0 = (abs(Es0-Es1)<=MatePairIndelLengthMax1) and c0==c1  ## sequences close enough to be considered coincident, with user defined slop value of MatePairIndelLengthMax1
            Em2 = (abs(Es2-Es3)<=MatePairIndelLengthMax1) and c2==c3
            ## User-specified conditions for skipping other read pairs
            VirtualFragLength1 = -Es0-Es2
            if CircD1[c0]:
                VirtualFragLength1 = VirtualFragLength1%LD1[c0]
            HomologyLength1 = ii1+klen1-ii0
            if HomologyLength1<MinHomology1 or HomologyLength1>MaxHomology1: continue
            if FivePrimeExtensionDisallowed1 and ii0>0: continue
            if StartHomology1 and not(StartHomology1==mySeq(SD1[ic0],AD1[ic0],CircD1[ic0],ip0,0,len(StartHomology1))): continue
            if AvoidR1Ambiguous1 and p0==0: continue
            if AvoidR2Ambiguous1 and p2==0: continue
            if AvoidBothAmbiguous1 and (p2==0 and p2==0): continue
            if SeparationMax1 or SeparationMin1:
                if c0!=c2 or p0*p2>=0: continue
                if SeparationMax1 and VirtualFragLength1>SeparationMax1: continue
                if SeparationMin1 and VirtualFragLength1<SeparationMin1: continue
            ## Now that conditions are met, record aspects of read pair
            PassFilters1 += 1
            if BriefCoverage1 or FullCoverage1:
                XA1.extend(X0)  
                XA1.extend(X2)             
            if ReportStarts1:
                if p0!=0:
                    CoverageS1[x0]= min(CoverageS1[x0]+1,CoverageSEMax1)  ## if there is a unique k-mer assign the start to the start of that k-mer; this results in a potential anomaly where the interface between repeated and unique sequence will serve as the effective start point for any read spanning the junction.  So caveat emptor!
                elif 0<=ix0<ulen1:
                    CoverageS1[ix0]= min(CoverageS1[ix0]+1,CoverageSEMax1) ## if there are no unique k-mers but there are repeated perfect matches, use the first such match to assign a provisional start
            if (SingleReadMode1==1) and ReportEnds1:
                if p1!=0:
                    CoverageE1[x1]= min(CoverageE1[x1]+1,CoverageSEMax1)
                elif 0<=ix1<ulen1:
                    CoverageE1[ix1]= min(CoverageE1[ix1]+1,CoverageSEMax1)
            if bbbT1 and ulen1>ix0>=0 and ii0>0 and (X0[-1]<ulen1 or KeepDoubleExtension1):
                if BaseByBase1 and ip0>0 and ((ic0,ip0) in BaseByBaseD1):
                    FivePrimeExtensionLength1 = min(ii0,MaxExtension1)
                    ProtoNubbin0 = mySeq(SD1[ic0],AD1[ic0],CircD1[ic0],ip0,-FivePrimeExtensionLength1,FivePrimeExtensionLength1)
                    Nubbin0,Match0,Mismatch0 = SimpleAlign1(L0[ii0-FivePrimeExtensionLength1:ii0],ProtoNubbin0)
                    FivePrimeExtension1 = "S5_"+Nubbin0
                    if Mismatch0>1 or Match0<SnpFilter1:
                        if not((ic0,ip0) in ExtensionD1):
                            ExtensionD1[(ic0,ip0)] = {}
                        if not FivePrimeExtension1 in ExtensionD1[(ic0,ip0)]:
                            ExtensionD1[(ic0,ip0)][FivePrimeExtension1] = 0
                        ExtensionD1[(ic0,ip0)][FivePrimeExtension1] += 1
                ip0A = -ip0 ##((-ip0-klen1) % LD1[ic1]) +1
                if BaseByBase1 and ip0<0 and ((ic0,ip0A) in BaseByBaseD1):                                
                    FivePrimeExtensionLength1 = min(ii0,MaxExtension1)
                    ProtoNubbin0 = mySeq(SD1[ic0],AD1[ic0],CircD1[ic0],-ip0A,-FivePrimeExtensionLength1,FivePrimeExtensionLength1)
                    Nubbin0,Match0,Mismatch0 = SimpleAlign1(L0[ii0-FivePrimeExtensionLength1:ii0],ProtoNubbin0)
                    if Mismatch0>1 or Match0<SnpFilter1:
                        FivePrimeExtension1 = "A5_"+Nubbin0                                     
                        if not((ic0,ip0A) in ExtensionD1):
                            ExtensionD1[(ic0,ip0A)] = {}
                        if not FivePrimeExtension1 in ExtensionD1[(ic0,ip0A)]:
                            ExtensionD1[(ic0,ip0A)][FivePrimeExtension1] = 0
                        ExtensionD1[(ic0,ip0A)][FivePrimeExtension1] += 1
            if bbbT1 and (SingleReadMode1==1) and ulen1>ix1>=0 and ii1<lL0-klen1 and (X0[0]<ulen1 or KeepDoubleExtension1):
                ip1S= ip1+klen1-1
                if BaseByBase1 and ip1>0 and (ic1,ip1S) in BaseByBaseD1:
                    ThreePrimeExtensionLength1 = min(lL0-klen1-ii1,MaxExtension1)
                    ProtoNubbin1 = mySeq(SD1[ic1],AD1[ic1],CircD1[ic1],ip1S,1,ThreePrimeExtensionLength1)
                    Nubbin1,Match1,Mismatch1 = SimpleAlign1(L0[ii1+klen1:ii1+klen1+ThreePrimeExtensionLength1],ProtoNubbin1)
                    if Mismatch1>1 or Match1<SnpFilter1:
                        ThreePrimeExtension1 = "S3_"+Nubbin1
                        if not((ic1,ip1S) in ExtensionD1):
                            ExtensionD1[(ic1,ip1S)] = {}
                        if not ThreePrimeExtension1 in ExtensionD1[(ic1,ip1S)]:
                            ExtensionD1[(ic1,ip1S)][ThreePrimeExtension1] = 0
                        ExtensionD1[(ic1,ip1S)][ThreePrimeExtension1] += 1
                ip1A = ((-ip1-klen1) %LD1[ic1])+1
                if BaseByBase1 and ip1<0 and (ic1,ip1A) in BaseByBaseD1:
                    ThreePrimeExtensionLength1 = min(lL0-klen1-ii1,MaxExtension1)
                    ProtoNubbin1 = mySeq(SD1[ic1],AD1[ic1],CircD1[ic1],-ip1A,1,ThreePrimeExtensionLength1)
                    Nubbin1,Match1,Mismatch1 = SimpleAlign1(L0[ii1+klen1:ii1+klen1+ThreePrimeExtensionLength1],ProtoNubbin1)
                    if Mismatch1>1 or Match1<SnpFilter1:
                        ThreePrimeExtension1 = "A3_"+Nubbin1
                        if not((ic1,ip1A) in ExtensionD1):
                            ExtensionD1[(ic1,ip1A)] = {}
                        if not ThreePrimeExtension1 in ExtensionD1[(ic1,ip1A)]:
                            ExtensionD1[(ic1,ip1A)][ThreePrimeExtension1] = 0
                        ExtensionD1[(ic1,ip1A)][ThreePrimeExtension1] += 1

            if (SingleReadMode1 == 0) and ReportEnds1:
                if p2!=0:
                    CoverageE1[x2]= min(CoverageE1[x2]+1,CoverageSEMax1)
                elif 0<=ix2<ulen1:
                    CoverageE1[ix2]= min(CoverageE1[ix2]+1,CoverageSEMax1)
            if bbbT1 and (SingleReadMode1 == 0) and ic0==ic2 and ip0>0 and abs(ip2+ip0)<ReadSeparationMax1:
                ip2S = -ip2  
                if BaseByBase1 and ip2<0 and ((ic2,ip2S) in BaseByBaseD1) and ii2>0:
                    ThreePrimeExtensionLength1 = min(ii2,MaxExtension1)
                    ProtoNubbin2 = antisense(mySeq(SD1[ic2],AD1[ic2],CircD1[ic2],ip2S,-ThreePrimeExtensionLength1,ThreePrimeExtensionLength1))
                    Nubbin2,Match2,Mismatch2 = SimpleAlign1(antisense(M0[ii2-MaxExtension1:ii2]),ProtoNubbin2)
                    if Mismatch2>1 or Match2<SnpFilter1:
                        ThreePrimeExtension1 = "S3_"+Nubbin2
                        if not((ic2,ip2) in ExtensionD1):
                            ExtensionD1[(ic2,ip2)] = {}
                        if not ThreePrimeExtension1 in ExtensionD1[(ic2,ip2)]:
                            ExtensionD1[(ic2,ip2)][ThreePrimeExtension1] = 0
                        ExtensionD1[(ic2,ip2)][ThreePrimeExtension1] += 1
                if BaseByBase1 and ip2>0 and ((ic2,ip2) in BaseByBaseD1):                                
                    ThreePrimeExtensionLength1 = min(ii2,MaxExtension1)
                    ProtoNubbin2 = antisense(mySeq(SD1[ic2],AD1[ic2],CircD1[ic2],ip2,-ThreePrimeExtensionLength1,ThreePrimeExtensionLength1))
                    Nubbin2,Match2,Mismatch2 = SimpleAlign1(antisense(M0[ii2-MaxExtension1:ii2]),ProtoNubbin2)
                    if Mismatch2>1 or Match2<SnpFilter1:
                        ThreePrimeExtension1 = "A3_"+Nubbin2
                        if not((ic2,ip2) in ExtensionD1):
                            ExtensionD1[(ic2,ip2)] = {}
                        if not ThreePrimeExtension1 in ExtensionD1[(ic2,ip2)]:
                            ExtensionD1[(ic2,ip2)][ThreePrimeExtension1] = 0
                        ExtensionD1[(ic2,ip2)][ThreePrimeExtension1] += 1
        ## does p0/p1 evidence a split read?
            ## The extrapolated start sites for each of the reads based on the first and last matching k-mer
            if FirstKTable1 or InterimKTable1:
                if ip0!=0:
                    J0A[min(klen1d,ii0)] += 1
                    J1A[min(klen1d,lL0-ii1-klen1)] += 1
                    if FirstKByStrand1:
                        if ip0>0:
                            JS0A[min(klen1d,ii0)] += 1
                        elif ip0<0:
                            JA0A[min(klen1d,ii0)] += 1
                        if ip1>0:
                            JS1A[min(klen1d,lL0-ii1-klen1)] += 1
                        elif ip1<0:
                            JA1A[min(klen1d,lL0-ii1-klen1)] += 1
                if ip2!=0:
                    J2A[min(klen1d,ii2)] += 1
                    J3A[min(klen1d,lM0-ii3-klen1)] += 1
                    if FirstKByStrand1:
                        if ip2<0:
                            JS2A[min(klen1d,ii2)] += 1
                        elif ip2>0:
                            JA2A[min(klen1d,ii2)] += 1
                        if ip3<0:
                            JS3A[min(klen1d,lM0-ii3-klen1)] += 1
                        elif ip3>0:
                            JA3A[min(klen1d,lM0-ii3-klen1)] += 1                
    ## First do some summary measurements of reads
            if (c0,c2,Es0,Es2) in D0:
                NewCombo1=False   ## True if this is a completely new combination
            else:
                NewCombo1=True
                D0[ (c0,c2,Es0,Es2) ] = [0,0]
            D0[ (c0,c2,Es0,Es2) ][0] += 1
            if p0!=0 or p2!=0: LocatableReadPairs1+=1
            if p0!=0 and p2!=0: DoubleLocatableReadPairs1+=1
            if p0!=0 and p2!=0 and Em0 and Em2: DefinitiveDoubleLocatableReadPairs1+=1
            if p0!=0:
                pIndex = min(abs(Es0)//SeparationGranularity1,len(PA1[0][c0])-1)
                PA1[0][c0][pIndex]+=1
                if Em0:
                    PA1[10][c0][pIndex]+=1
                if NewCombo1:
                    PA1[4][c0][pIndex]+=1                    
            if p2!=0:
                pIndex = min(abs(Es2)//SeparationGranularity1,len(PA1[0][c2])-1)
                PA1[0][c2][pIndex]+=1
                if Em2:
                    PA1[10][c2][pIndex]+=1
                if NewCombo1:
                    PA1[4][c2][pIndex]+=1                    
            if DFAMCount1:
                if ip0!=0:
                    ##v0a=FastFindR1(sc0,abs(rp0))  ## Slower than the explicit array but more memory parsimonious
                    v0=DfamQ1[ic0][abs(ip0)+klen2-1]
                    if v0>0:
                        DfamC1[v0-1]+=1
                if ip2!=0:
                    ##v2=FastFindR1(sc2,abs(rp2)) ## see comment above
                    v2=DfamQ1[ic2][abs(ip2)+klen2-1]
                    if v2>0:
                        DfamC1[v2-1]+=1
            if c0==c2 and p0!=0 and p2!=0 and -Tn5DupMax1<=(-Es0-Es2)<ReadSeparationMax1:
                SeparationArray1[c0][-Es0-Es2+Tn5DupMax1]+=1
            elif SingleReadMode1 == 1  and HomologyLength1<ReadSeparationMax1:
                SeparationArray1[c0][HomologyLength1]+=1
            if c0==c2 and 7<=-Es0-Es2<=9:
                p20=(abs(Es0)+abs(Es2))//2
                pIndex = min(p20//SeparationGranularity1,len(PA1[0][c0])-1)
                if p0!=0 and p2!=0:
                    PA1[18][c0][pIndex]+=1
                    if NewCombo1:
                        PA1[19][c0][pIndex]+=1
                PA1[20][c0][pIndex]+=1
                if NewCombo1:
                    PA1[21][c0][pIndex]+=1
            if ix0>=0:
                im0 = Multiplicities1[ix0] & Mult1
                if not im0 in MultiplicityD1[ic0]:
                    MultiplicityD1[ic0][im0]=0
                MultiplicityD1[ic0][im0]+=1
            if ix2>=0:
                im2 = Multiplicities1[ix2] & Mult1
                if not im2 in MultiplicityD1[ic2]:
                    MultiplicityD1[ic2][im2]=0
                MultiplicityD1[ic2][im2]+=1
            if p0==0 and ip0!=0:
                pIndex = min(abs(ip0)//SeparationGranularity1,len(PA1[0][ic0])-1)
                if Multiplicities1[ix0] & Local1:
                    PA1[1][ic0][pIndex]+=1
                    if NewCombo1:
                        PA1[5][ic0][pIndex]+=1
                elif Multiplicities1[ix0] & Chromosomal1:
                    PA1[2][ic0][pIndex]+=1
                    if NewCombo1:
                        PA1[6][ic0][pIndex]+=1
                else:
                    PA1[3][ic0][pIndex]+=1                
                    if NewCombo1:
                        PA1[7][ic0][pIndex]+=1
            if p2==0 and ip2!=0:
                pIndex = min(abs(ip2)//SeparationGranularity1,len(PA1[0][ic2])-1)
                if Multiplicities1[ix2] & Local1:
                    PA1[1][ic2][pIndex]+=1
                    if NewCombo1:
                        PA1[5][ic2][pIndex]+=1
                elif Multiplicities1[ix2] & Chromosomal1:
                    PA1[2][ic2][pIndex]+=1                
                    if NewCombo1:
                        PA1[6][ic2][pIndex]+=1
                else:
                    PA1[3][ic2][pIndex]+=1                
                    if NewCombo1:
                        PA1[7][ic2][pIndex]+=1
            ## Look for "Normal" read pairs and place them in a dictionary-- goal here is to get a rough idea of redundancy in the sequencing dataset
            if ((p0<0 and p2>0) or (p0>0 and p2<0)) and c0==c2 and Tn5DupMax1<=-Es0-Es2<=ReadSeparationMax1:
                pIndex = min(((abs(Es0)+abs(Es2))//2)//SeparationGranularity1,len(PA1[0][c0])-1)
                PA1[8][c0][pIndex]+=1
                if Em0 and Em2: 
                    D0[ (c0,c2,Es0,Es2) ][1] +=1
                if NewCombo1:
                    PA1[9][c0][pIndex]+=1
            if FindStructuralAnomalies1:
                MPLS11 = MPLS1 or (p1>=p0+klen1 and p3>=p2+klen1) ## if MPLS1 is False, use high stringency matching (two nonoverlapping unique k-mers miust match to call junction
                MinCirDiff1 = max(lL0-Tn5DupMax1,InsertionMin1)
                MinCirDiff2 = max(lM0-Tn5DupMax1,InsertionMin1)
                Events1[0] = p0>0 and p2<0 and 0<Es0+Es2<CircleMax1 and Em0 and Em2 and c0==c2 and MPLS11
                Events1[1] = p0<0 and p2>0 and 0<Es0+Es2<CircleMax1 and Em0 and Em2 and c0==c2 and MPLS11  
                Events1[2] = p0>0 and p2<0 and Tn5DupMin1<=-Es0-Es2<=Tn5DupMax1 and Em0 and Em2 and c0==c2 and MPLS11        
                Events1[3] = p0<0 and p2>0 and Tn5DupMin1<=-Es0-Es2<=Tn5DupMax1 and Em0 and Em2 and c0==c2 and MPLS11
                Events1[4] = p0>0 and p2<0 and DeletionMax1>=-Es0-Es2>ReadSeparationMax1 and Em0 and Em2 and c0==c2 and MPLS11
                Events1[5] = p0<0 and p2>0 and DeletionMax1>=-Es0-Es2>ReadSeparationMax1 and Em0 and Em2 and c0==c2 and MPLS11
                Events1[6] = p0>0 and p1>0 and c0==c1 and CircleMax1>=Es0-Es1>=MinCirDiff1
                Events1[7] = p2>0 and p3>0 and c2==c3 and CircleMax1>=Es2-Es3>=MinCirDiff2
                Events1[8] = p0<0 and p1<0 and c0==c1 and CircleMax1>=Es0-Es1>=MinCirDiff1
                Events1[9] = p2<0 and p3<0 and c2==c3 and CircleMax1>=Es2-Es3>=MinCirDiff2
                Events1[10] = p0>0 and p1>0 and c0==c1 and DeletionMax1>=Es1-Es0>=MinCirDiff1
                Events1[11] = p2>0 and p3>0 and c2==c3 and DeletionMax1>=Es3-Es2>=MinCirDiff2
                Events1[12] = p0<0 and p1<0 and c0==c1 and DeletionMax1>=Es1-Es0>=MinCirDiff1
                Events1[13] = p2<0 and p3<0 and c2==c3 and DeletionMax1>=Es3-Es2>=MinCirDiff2
                Events1[14] = p0>0 and p1>0 and c0==c1 and MinCirDiff1>Es0-Es1>=InsertionMin1
                Events1[15] = p2>0 and p3>0 and c2==c3 and MinCirDiff2>Es2-Es3>=InsertionMin1
                Events1[16] = p0<0 and p1<0 and c0==c1 and MinCirDiff1>Es0-Es1>=InsertionMin1
                Events1[17] = p2<0 and p3<0 and c2==c3 and MinCirDiff2>Es2-Es3>=InsertionMin1
                if OtherSVCandidates1:
                    Events1[18] = p0>0 and p2>0 and Em0 and Em2 and c0==c2 and MPLS11  ## 'Inversion_Candidate_1s2s'
                    Events1[19] = p0<0 and p2<0 and Em0 and Em2 and c0==c2 and MPLS11  ## 'Inversion_Candidate_1a2a'
                    Events1[20] = p0>0 and p1<0 and c0==c1 ## 'Inversion_Candidate_1s1a'
                    Events1[21] = p2>0 and p3<0 and c2==c3 ## 'Inversion_Candidate_2s2a',
                    Events1[22] = p0<0 and p1>0 and c0==c1 ## 'Inversion_Candidate_1a1s'
                    Events1[23] = p2<0 and p3>0 and c2==c3 ## 'Inversion_Candidate_2a2s',
                    Events1[24] = p0>0 and p2>0 and Em0 and Em2 and c0!=c2 and MPLS11 ## 'Translocation_Candidate_1s2s'
                    Events1[25] = p0>0 and p2<0 and Em0 and Em2 and c0!=c2 and MPLS11 ## 'Translocation_Candidate_1s2a'
                    Events1[26] = p0<0 and p2>0 and Em0 and Em2 and c0!=c2 and MPLS11 ## 'Translocation_Candidate_1a2s'
                    Events1[27] = p0<0 and p2<0 and Em0 and Em2 and c0!=c2 and MPLS11 ## 'Translocation_Candidate_1a2a'
                    Events1[28] = p0>0 and p1>0 and c0!=c1 ## 'Translocation_Candidate_1s1s'
                    Events1[29] = p0>0 and p1<0 and c0!=c1 ## 'Translocation_Candidate_1s1a'
                    Events1[30] = p0<0 and p1>0 and c0!=c1 ## 'Translocation_Candidate_1a1s'
                    Events1[31] = p0<0 and p1<0 and c0!=c1 ## 'Translocation_Candidate_1a1a'
                    Events1[32] = p2>0 and p3>0 and c2!=c3 ## 'Translocation_Candidate_2s2s'
                    Events1[33] = p2>0 and p3<0 and c2!=c3 ## 'Translocation_Candidate_2s2a'
                    Events1[34] = p2<0 and p3>0 and c2!=c3 ## 'Translocation_Candidate_2a2s'
                    Events1[35] = p2<0 and p3<0 and c2!=c3 ## 'Translocation_Candidate_2a2a'

                if any(Events1[0:6]) or (OtherSVCandidates1 and any(Events1[18:20]+Events1[24:28])):
                    ESeq0=ESeq1=mySeq(SD1[c0],AD1[c0],CircD1[c0],p0,-i0,lL0)
                    ESeq2=ESeq3=mySeq(SD1[c2],AD1[c2],CircD1[c2],p2,-i2,lM0)##corrected with 6.20.17 version of MySeq call (do not remove this comment)
                    ScoreL1,MatchL1,MisMatchL1,IndelS1,IndelE1 = nw1(L0,ESeq0,span=MatePairIndelLengthMax1)
                    ScoreL2,MatchL2,MisMatchL2,IndelS2,IndelE2 = nw1(M0,ESeq2,span=MatePairIndelLengthMax1)
                    if (len(ESeq0)<lL0-MatePairIndelLengthMax1 or
                        MisMatchL1>MatePairSubstitutionMax1 or
                        IndelS1>MatePairIndelCountMax1 or
                        IndelE1>MatePairIndelLengthMax1 or
                        len(ESeq2)<lM0-MatePairIndelLengthMax1 or
                        MisMatchL2>MatePairSubstitutionMax1 or
                        IndelS2>MatePairIndelCountMax1 or
                        IndelE2>MatePairIndelLengthMax1):
                        Events1[:6]=[False]*6
                    Hl0=Hl1=MatchL1
                    Hl2=Hl3=MatchL2
                elif any(Events1[6:18]) or (OtherSVCandidates1 and any(Events1[20:24]+Events1[28:36])):
                    R1Events1 = list(range(6,17,2))
                    R2Events1 = list(range(7,18,2))
                    if OtherSVCandidates1:
                        R1Events1 += list(range(20,23,2))+list(range(28,35,2))
                        R2Events1 += list(range(21,24,2))+list(range(29,36,2))
                    ESeq0=mySeq(SD1[c0],AD1[c0],CircD1[c0],p0,-i0,lL0)
                    ESeq1=mySeq(SD1[c1],AD1[c1],CircD1[c1],p1,-i1,lL0)
                    ESeq2=mySeq(SD1[c2],AD1[c2],CircD1[c2],p2,-i2,lM0)
                    ESeq3=mySeq(SD1[c3],AD1[c3],CircD1[c3],p3,-i3,lM0) ## Extrapolated sequences (what the full read should have looked like) for each unique match##corrected with 6.20.17 version of MySeq call (do not remove this comment)
                    if RequireKMerOnly1:
                        Hl0 = i0+klen1; Hl1 = lL0-i1; Hl2 = i2+klen1; Hl3 = lM0-i3
                    else:
                        Hl0=0; Hl1=0; Hl2=0; Hl3=0  ## Homology Lengths
                    if Hl0<lL0:
                        while lL0 and (L0[Hl0]==ESeq0[Hl0]):
                            Hl0+=1
                            if Hl0==lL0:
                                break
                    if Hl1<lL0:
                        while lL0 and (L0[-(Hl1+1)]==ESeq1[-(Hl1+1)]):
                            Hl1+=1            
                            if Hl1==lL0:
                                break
                    if Hl2<lM0:
                        while lM0 and (M0[Hl2]==ESeq2[Hl2]):
                            Hl2+=1
                            if Hl2==lM0:
                                break
                    if Hl3<lM0:
                        while lM0 and (M0[-(Hl3+1)]==ESeq3[-(Hl3+1)]):
                            Hl3+=1
                            if Hl3==lM0:
                                break
                    if RequireKMerOnly1:
                        if any(Events1[6:17:2]) or (OtherSVCandidates1 and (any(Events1[20:23:2]+Events1[28:35:2]))):
                            nwANum0 = nw1(L0[:i0+klen1],ESeq0[:i0+klen1],span=SplitReadIndelLengthMax1)
                            nwANum1 = nw1(L0[i1:],ESeq1[i1:],span=SplitReadIndelLengthMax1)
                            if nwANum0[2]>SplitReadSubstitutionMax1 or nwANum0[3]>SplitReadIndelCountMax1 or nwANum1[2]>SplitReadSubstitutionMax1 or nwANum1[3]>SplitReadIndelCountMax1:
                                Hl0=0; Hl1=0
                        if any(Events1[7:18:2]) or (OtherSVCandidates1 and (any(Events1[21:24:2]+Events1[29:36:2]))):
                            nwANum2 = nw1(M0[:i2+klen1],ESeq2[:i2+klen1],span=SplitReadIndelLengthMax1)
                            nwANum3 = nw1(M0[i3:],ESeq3[i3:],span=SplitReadIndelLengthMax1)
                            if nwANum2[2]>SplitReadSubstitutionMax1 or nwANum2[3]>SplitReadIndelCountMax1 or nwANum3[2]>SplitReadSubstitutionMax1 or nwANum3[3]>SplitReadIndelCountMax1:
                                Hl2=0; Hl3=0
                    for g11 in R1Events1:
                        if Events1[g11] and ((Hl0<i0+klen1) or (Hl1<lL0-i1)):
                            Events1[g11]=False
                    for g11 in R2Events1:
                        if Events1[g11] and ((Hl2<i2+klen1) or (Hl3<lM0-i3)):
                            Events1[g11]=False
                else:
                    continue
                if any(Events1):
        ##         ['Circle_Discordant_1s2a', 'Circle_Discordant_1a2s', 'Circle_Contraction_1s2a', 'Circle_Contraction_1a2s',
        ##         'Deletion_Contraction_1s2a', 'Deletion_Contraction_1a2s',
        ##         'Circle_Junction_1s1s', 'Circle_Junction_2s2s', 'Circle_Junction_1a1a', 'Circle_Junction_2a2a',
        ##         'Deletion_Junction_1s1s', 'Deletion_Junction_2s2s', 'Deletion_Junction_1a1a', 'Deletion_Junction_2a2a',
        ##         'Insertion_Junction_Candidate_1s1s', 'Insertion_Junction_Candidate_2s2s', 'Insertion_Junction_Candidate_1a1a', 'Insertion_Junction_Candidate_2a2a']


        ##Events1[18] = p0>0 and p2>0 and Em0 and Em2 and c0==c2 and MPLS11  ## 'Inversion_Candidate_1s2s'
        ##Events1[19] = p0<0 and p2<0 and Em0 and Em2 and c0==c2 and MPLS11  ## 'Inversion_Candidate_1a2a'
        ##Events1[20] = p0>0 and p1<0 and c0==c1 ## 'Inversion_Candidate_1s1a'
        ##Events1[21] = p2>0 and p3<0 and c2==c3 ## 'Inversion_Candidate_2s2a',
        ##Events1[22] = p0<0 and p1>0 and c0==c1 ## 'Inversion_Candidate_1a1s'
        ##Events1[23] = p2<0 and p3>0 and c2==c3 ## 'Inversion_Candidate_2a2s',
        ##Events1[24] = p0>0 and p2>0 and Em0 and Em2 and c0!=c2 and MPLS11 ## 'Translocation_Candidate_1s2s'
        ##Events1[25] = p0>0 and p2<0 and Em0 and Em2 and c0!=c2 and MPLS11 ## 'Translocation_Candidate_1s2a'
        ##Events1[26] = p0<0 and p2>0 and Em0 and Em2 and c0!=c2 and MPLS11 ## 'Translocation_Candidate_1a2s'
        ##Events1[27] = p0<0 and p2<0 and Em0 and Em2 and c0!=c2 and MPLS11 ## 'Translocation_Candidate_1a2a'
        ##Events1[28] = p0>0 and p1>0 and c0!=c1 ## 'Translocation_Candidate_1s1s'
        ##Events1[29] = p0>0 and p1<0 and c0!=c1 ## 'Translocation_Candidate_1s1a'
        ##Events1[30] = p0<0 and p1>0 and c0!=c1 ## 'Translocation_Candidate_1a1s'
        ##Events1[31] = p0<0 and p1<0 and c0!=c1 ## 'Translocation_Candidate_1a1a'
        ##Events1[32] = p2>0 and p3>0 and c2!=c3 ## 'Translocation_Candidate_2s2s'
        ##Events1[33] = p2>0 and p3<0 and c2!=c3 ## 'Translocation_Candidate_2s2a'
        ##Events1[34] = p2<0 and p3>0 and c2!=c3 ## 'Translocation_Candidate_2a2s'
        ##Events1[35] = p2<0 and p3<0 and c2!=c3 ## 'Translocation_Candidate_2a2a'
                    if OtherSVCandidates1:
                        if Events1[18]:
                            EventChr=c0
                            EventStartLo=Es0+lL0
                            EventStartHi=Es0+ReadSeparationMax1
                            EventEndLo=Es2+lM0
                            EventEndHi=Es2+ReadSeparationMax1
                            EventOverlap=EventEndLo-EventStartLo
                        if Events1[19]:
                            EventChr=c0
                            EventStartLo=-Es0-ReadSeparationMax1
                            EventStartHi=-Es0-lL0
                            EventEndLo=-Es2-ReadSeparationMax1
                            EventEndHi=-Es2-lM0
                            EventOverlap=EventEndLo-EventStartLo
                        if Events1[24]:
                            EventChr=(c0,c2)
                            EventStartLo=Es0+lL0
                            EventStartHi=Es0+ReadSeparationMax1
                            EventEndLo=Es2+lM0
                            EventEndHi=Es2+ReadSeparationMax1
                            EventOverlap=0
                        if Events1[25]:
                            EventChr=(c0,c2)
                            EventStartLo=Es0+lL0
                            EventStartHi=Es0+ReadSeparationMax1
                            EventEndLo=-Es2-ReadSeparationMax1
                            EventEndHi=-Es2-lM0
                            EventOverlap=0
                        if Events1[26]:
                            EventChr=(c0,c2)
                            EventStartLo=-Es0-ReadSeparationMax1
                            EventStartHi=-Es0-lL0
                            EventEndLo=Es2+lM0
                            EventEndHi=Es2+ReadSeparationMax1
                            EventOverlap=0
                        if Events1[27]:
                            EventChr=(c0,c2)
                            EventStartLo=-Es0-ReadSeparationMax1
                            EventStartHi=-Es0-lL0
                            EventEndLo=-Es2-ReadSeparationMax1
                            EventEndHi=-Es2-lM0
                            EventOverlap=0
                        if Events1[20]:
                            EventChr=c0
                            EventStartHi=Es0+Hl0
                            EventStartLo=min(EventStartHi, EventStartHi-(Hl0+Hl1-lL0))
                            EventEndLo=-Es1-Hl1
                            EventEndHi=max(EventEndLo, EventEndLo+(Hl0+Hl1-lL0))
                            EventOverlap=Hl0+Hl1-lL0
                        if Events1[21]:
                            EventChr=c2
                            EventStartHi=Es2+Hl2
                            EventStartLo=min(EventStartHi, EventStartHi-(Hl2+Hl3-lM0))
                            EventEndLo=-Es3-Hl3
                            EventEndHi=max(EventEndLo, EventEndLo+(Hl2+Hl3-lM0))
                            EventOverlap=Hl2+Hl3-lM0
                        if Events1[22]:
                            EventChr=c0
                            EventStartLo=-Es0-Hl0
                            EventStartHi=max(EventStartLo, EventStartLo+(Hl0+Hl1-lL0))
                            EventEndHi=Es1+Hl1
                            EventEndLo=min(EventEndHi, EventEndHi-(Hl0+Hl1-lL0))
                            EventOverlap=Hl0+Hl1-lL0
                        if Events1[23]:
                            EventChr=c2
                            EventStartLo=-Es2-Hl2
                            EventStartHi=max(EventStartLo, EventStartLo+(Hl2+Hl3-lM0))
                            EventEndHi=Es3+Hl3
                            EventEndLo=min(EventEndHi, EventEndHi-(Hl2+Hl3-lM0))
                            EventOverlap=Hl2+Hl3-lM0
                        if Events1[28]:
                            EventChr=(c0,c1)
                            EventStartHi=Es0+Hl0
                            EventStartLo=min(EventStartHi, EventStartHi-(Hl0+Hl1-lL0))
                            EventEndHi=Es1+Hl1
                            EventEndLo=min(EventEndHi, EventEndHi-(Hl0+Hl1-lL0))
                            EventOverlap=Hl0+Hl1-lL0
                        if Events1[29]:
                            EventChr=(c0,c1)
                            EventStartHi=Es0+Hl0
                            EventStartLo=min(EventStartHi, EventStartHi-(Hl0+Hl1-lL0))
                            EventEndLo=-Es1-Hl1
                            EventEndHi=max(EventEndLo, EventEndLo+(Hl0+Hl1-lL0))
                            EventOverlap=Hl0+Hl1-lL0
                        if Events1[30]:
                            EventChr=(c0,c1)
                            EventStartLo=-Es0-Hl0
                            EventStartHi=max(EventStartLo, EventStartLo+(Hl0+Hl1-lL0))
                            EventEndHi=Es1+Hl1
                            EventEndLo=min(EventEndHi, EventEndHi-(Hl0+Hl1-lL0))
                            EventOverlap=Hl0+Hl1-lL0
                        if Events1[31]:
                            EventChr=(c0,c1)
                            EventStartLo=-Es0-Hl0
                            EventStartHi=max(EventStartLo, EventStartLo+(Hl0+Hl1-lL0))
                            EventEndLo=-Es1-Hl1
                            EventEndHi=max(EventEndLo, EventEndLo+(Hl0+Hl1-lL0))
                            EventOverlap=Hl0+Hl1-lL0
                        if Events1[32]:
                            EventChr=(c2,c3)
                            EventStartHi=Es2+Hl2
                            EventStartLo=min(EventStartHi, EventStartHi-(Hl2+Hl3-lM0))
                            EventEndHi=Es3+Hl3
                            EventEndLo=min(EventEndHi, EventEndHi-(Hl2+Hl3-lM0))
                            EventOverlap=Hl2+Hl3-lM0
                        if Events1[33]:
                            EventChr=(c2,c3)
                            EventStartHi=Es2+Hl2
                            EventStartLo=min(EventStartHi, EventStartHi-(Hl2+Hl3-lM0))
                            EventEndLo=-Es3-Hl3
                            EventEndHi=max(EventEndLo, EventEndLo+(Hl2+Hl3-lM0))
                            EventOverlap=Hl2+Hl3-lM0
                        if Events1[34]:
                            EventChr=(c2,c3)
                            EventStartLo=-Es2-Hl2
                            EventStartHi=max(EventStartLo, EventStartLo+(Hl2+Hl3-lM0))
                            EventEndHi=Es3+Hl3
                            EventEndLo=min(EventEndHi, EventEndHi-(Hl2+Hl3-lM0))
                            EventOverlap=Hl2+Hl3-lM0
                        if Events1[35]:
                            EventChr=(c2,c3)
                            EventStartLo=-Es2-Hl2
                            EventStartHi=max(EventStartLo, EventStartLo+(Hl2+Hl3-lM0))
                            EventEndLo=-Es3-Hl3
                            EventEndHi=max(EventEndLo, EventEndLo+(Hl2+Hl3-lM0))
                            EventOverlap=Hl2+Hl3-lM0                       
                    if Events1[0] or Events1[2]:
                        EventChr=c0
                        EventEndLo=Es0+lL0
                        EventEndHi=Es0+ReadSeparationMax1
                        EventStartHi=-Es2-lM0
                        EventStartLo=-Es2-ReadSeparationMax1
                        EventOverlap=-Es2-Es0
                    if Events1[1] or Events1[3]:
                        EventChr=c0
                        EventEndLo=Es2+lM0
                        EventEndHi=Es2+ReadSeparationMax1
                        EventStartHi=-Es0-lL0
                        EventStartLo=-Es0-ReadSeparationMax1
                        EventOverlap=-Es0-Es2
                    if Events1[4]:
                        EventChr=c0
                        EventStartLo=Es0+lL0
                        EventEndHi=-Es2-lM0
                        EventStartHi=min(Es0+ReadSeparationMax1,EventEndHi)
                        EventEndLo=max(-Es2-ReadSeparationMax1,EventStartLo)
                        EventOverlap=-Es0-Es2
                    if Events1[5]:
                        EventChr=c0
                        EventStartLo=Es2+lM0
                        EventEndHi=-Es0-lL0
                        EventStartHi=min(Es2+ReadSeparationMax1,EventEndHi)
                        EventEndLo=max(-Es0-ReadSeparationMax1,EventStartLo)
                        EventOverlap=-Es0-Es2
                    if Events1[6] or Events1[14]:
                        EventChr=c0
                        EventEndHi=Es0+Hl0
                        EventStartLo=Es1+lL0-Hl1
                        EventEndLo=min(EventEndHi,EventEndHi-(Hl0+Hl1-lL0))
                        EventStartHi=max(EventStartLo,EventStartLo+(Hl0+Hl1-lL0))
                        EventOverlap=Hl0+Hl1-lL0
                    if Events1[7] or Events1[15]:
                        EventChr=c2
                        EventEndHi=Es2+Hl2
                        EventStartLo=Es3+lM0-Hl3
                        EventEndLo=min(EventEndHi,EventEndHi-(Hl2+Hl3-lM0))
                        EventStartHi=max(EventStartLo,EventStartLo+(Hl2+Hl3-lM0))
                        EventOverlap=Hl2+Hl3-lM0
                    if Events1[8] or Events1[16]:
                        EventChr=c0
                        EventStartLo=-Es0-Hl0
                        EventEndHi=-Es1-lL0+Hl1
                        EventStartHi=max(EventStartLo,EventStartLo+(lL0-Hl0-Hl1))
                        EventEndLo=min(EventEndHi,EventEndHi-(lL0-Hl0-Hl1))
                        EventOverlap=Hl0+Hl1-lL0
                    if Events1[9] or Events1[17]:
                        EventChr=c2
                        EventStartLo=-Es2-Hl2
                        EventEndHi=-Es3-lM0+Hl3
                        EventStartHi=max(EventStartLo,EventStartLo+(lM0-Hl2-Hl3))
                        EventEndLo=min(EventEndHi,EventEndHi-(lM0-Hl2-Hl3))
                        EventOverlap=Hl2+Hl3-lM0
                    if Events1[10]:
                        EventChr=c0
                        EventStartHi=Es0+Hl0
                        EventEndLo=Es1+lL0-Hl1
                        EventStartLo=min(EventStartHi,EventStartHi-(Hl0+Hl1-lL0))
                        EventEndHi=max(EventEndLo,EventEndLo+(Hl0+Hl1-lL0))
                        EventOverlap=Hl0+Hl1-lL0
                    if Events1[11]:
                        EventChr=c2
                        EventStartHi=Es2+Hl2
                        EventEndLo=Es3+lM0-Hl3
                        EventStartLo=min(EventStartHi,EventStartHi-(Hl2+Hl3-lM0))
                        EventEndHi=max(EventEndLo,EventEndLo+(Hl2+Hl3-lM0))
                        EventOverlap=Hl2+Hl3-lM0
                    if Events1[12]:
                        EventChr=c0
                        EventStartHi=-Es1-lL0+Hl1  
                        EventEndLo=-Es0-Hl0
                        EventStartLo=min(EventStartHi,EventStartHi-(lL0-Hl0-Hl1))
                        EventEndHi=max(EventEndLo,EventEndLo+(lL0-Hl0-Hl1))
                        EventOverlap=Hl0+Hl1-lL0
                    if Events1[13]:
                        EventChr=c2
                        EventStartHi=-Es3-lM0+Hl3  
                        EventEndLo=-Es2-Hl2
                        EventStartLo=min(EventStartHi,EventStartHi-(lM0-Hl2-Hl3))
                        EventEndHi=max(EventEndLo,EventEndLo+(lM0-Hl2-Hl3))
                        EventOverlap=Hl2+Hl3-lM0
        ##        Events0=['Circle_Discordant_1s2a', 'Circle_Discordant_1a2s', 'Circle_Contraction_1s2a', 'Circle_Contraction_1a2s',
        ##                 'Deletion_Contraction_1s2a', 'Deletion_Contraction_1a2s',
        ##                 'Circle_Junction_1s1s', 'Circle_Junction_2s2s', 'Circle_Junction_1a1a', 'Circle_Junction_2a2a',
        ##                 'Deletion_Junction_1s1s', 'Deletion_Junction_2s2s', 'Deletion_Junction_1a1a', 'Deletion_Junction_2a2a'
        ##                 'Insertion_Junction_1s1s', 'Insertion_Junction_2s2s', 'Insertion_Junction_1a1a', 'Insertion_Junction_2a2a']
                    AllEvents0=[]
                    for ev0,ev1 in enumerate(Events1):
                        if ev1:
                            AllEvents0.append(Events0[ev0])
                    AllEvents1=','.join(AllEvents0)
                    if type(EventChr)==tuple:
                        Event11 = NameA1[EventChr[0]]+'_to_'+NameA1[EventChr[1]]
                        EventChr = EventChr[0]
                    else:
                        Event11 = NameA1[EventChr]
                    BasicInfo1=[AllEvents1,
                                Event11,
                                EventStartLo,
                                EventStartHi,
                                EventEndLo,
                                EventEndHi,
                                EventOverlap,
                                UCSClink1(NameA1[EventChr],EventStartLo-UCSCBuffer1,EventEndHi+UCSCBuffer1)]
                    F3.write('\t'.join(map(str,BasicInfo1)))
                    if VerboseOutput1:
                        MoreInfo1=[]
                        if Read1Input:
                            MoreInfo1 += [L0,
                                    NameA1[c0],
                                    Es0,
                                    Hl0,
                                    ESeq0,
                                    NameA1[c1],
                                    Es1,
                                    Hl1,
                                    ESeq1,
                                    Es1-Es0,
                                    -lL0+Hl0+Hl1]
                        if Read2Input:
                            MoreInfo1 += [M0,
                                NameA1[c2],
                                Es2,
                                Hl2,
                                ESeq2,
                                NameA1[c3],
                                Es3,
                                Hl3,
                                ESeq3,
                                Es3-Es2,
                                -lM0+Hl2+Hl3]
                        if LastL0.startswith('>'):
                            MoreInfo1.append(LastL0.strip()[1:])
                        elif LastM0.startswith('>'):
                            MoreInfo1.append(LastM0.strip()[1:])
                        else:
                            MoreInfo1.append(Mnemonic2+':'+str(filelinenum//LineDensity1))
                        F3.write('\t'+'\t'.join(map(str,MoreInfo1)))
                    if ReportBuffers1:
                        if R1Buffer5:
                            F3.write('\t'+L00.strip()[:R1Buffer5])
                        if R1Buffer3:
                            F3.write('\t'+L00.strip()[-R1Buffer3:])
                        if R2Buffer5:
                            F3.write('\t'+M00.strip()[:R2Buffer5])
                        if R2Buffer3:
                            F3.write('\t'+M00.strip()[-R2Buffer3:])
                    F3.write('\n')
        LastL0,LastM0 = L0,M0
    if ( XA1 != [] ) and (BriefCoverage1 or FullCoverage1):
        LogNote1("Starting End-of-file Coverage Calculations for read "+str(exptlinenum1),LogFile1)
        XA2=np.unique(np.asarray(XA1,dtype=FastIndexDataType1),return_counts=True)
        if FullCoverage1:
            CoverageK1[XA2[0]]=np.minimum(CoverageK1[XA2[0]]+XA2[1],MaxCoverage1)
        XA1=[]
        if BriefCoverage1:
            IA1=Sort_I1dedup[XA2[0]]
            PAx1=Sort_P1dedup[XA2[0]]
            MA1=Multiplicities1[XA2[0]]
            for iA1 in range(NumSeq1):
                PA2=PAx1[ (IA1==iA1) & (MA1==1)]
                PA3=np.unique( ((np.abs(2*PA2+klen1-1)-klen1-1)>>1) % LD1[iA1] )
                EitherStrandUnique1[iA1]+=len(PA3)-np.count_nonzero(CoverageU1[iA1][PA3])
                CoverageU1[iA1][PA3]=1
            EitherStrandUnique2=prependtotal(EitherStrandUnique1)
            LabelID1='EOF'
            if f1==FiLR1[-1] and f2==FiLR2[-1]:
                LabelID1='Complete'
            F5.write('UniqueCoverageCount'+LabelID1+'\t'+str(exptlinenum1)+'\t'+'\t'.join(map(str,EitherStrandUnique2))+'\t0\n')
            M1=per03(means1(TotalUnique1,EitherStrandUnique2))
            F5.write('UniqueCoveragePercent'+LabelID1+'\t'+str(exptlinenum1)+'\t'+'\t'.join(M1)+'\t0\n')
            if M1:
                LogNote1('After '+str(exptlinenum1)+' read pairs, unique coverage stands at '+M1[0]+'% ',LogFile1)
            else:
                LogNote1('After '+str(exptlinenum1)+' read pairs, unique coverage stands at '+'0.000'+'% ',LogFile1)
        LogNote1("Finished Inerim Coverage Calculations for read "+str(exptlinenum1),LogFile1)
    if FullCoverageByFile1 and f1!=FiLR1[-1] and f2!=FiLR2[-1]:
        MilestoneL1=Mnemonic2+'_'+str(filelinenum//LineDensity1)+'_'+str(exptlinenum1)
        if IndexConstructionBins1>1:
            npush('Sort_V1dedup') ## release this memory before calling CalculateFullCoverage
        CalculateFullCoverage1(MilestoneL1)
        if IndexConstructionBins1>1:
            npull('Sort_V1dedup')
    elif FullCoverage1 and f1==FiLR1[-1] and f2==FiLR2[-1]:
        MilestoneL1=Mnemonic1+'_'+str(filelinenum//LineDensity1)+'_'+str(exptlinenum1)+'_final'
        if IndexConstructionBins1>1:
            npush('Sort_V1dedup')
        CalculateFullCoverage1(MilestoneL1)
    if 'file' in str(type(F1)).lower():
        F1.close()
        LogClosingFile1(F1)
    if 'file' in str(type(F2)).lower():
        F2.close()
        LogClosingFile1(F2)         

                    

LogNote1("Total Read Pairs: "+str(TotalReadPairs1),LogFile1)
LogNote1("Locatable Read Pairs: "+str(LocatableReadPairs1),LogFile1)
LogNote1("Double Locatable Read Pairs: "+str(DoubleLocatableReadPairs1),LogFile1)
LogNote1("Definitive Double Locatable Read Pairs: "+str(DefinitiveDoubleLocatableReadPairs1),LogFile1)
LogNote1("Starting Summary Calculations",LogFile1)

if CombinationList1:
    F4=open("PairedEndCombinations_"+OutFileNameBase+'.tdv',mode='w')
    LogOpeningFile1(F4)
    F4.write('<!--REVA Counts of paired end r1:r2 start combinations for -->\n')  ## write a list of 'chromosome' (reference sequence strings) as number (0 based) and name
    F4.write('\n')
    F4.write('\tChromosomeNumber\tChromosomeName\tChromosomeLength'+'\n')  ## write a list of 'chromosome' (reference sequence strings) as number (0 based) and name
    for iu1,nu1 in enumerate(NameA1):
        F4.write('\t'+str(iu1)+'\t'+nu1+'\t'+str(LD1[iu1])+'\n')  ## write a list of 'chromosome' (reference sequence strings) as number (0 based) and name
    OutputHeaderT4='\t'.join(['r1_Chromosome__'+RefAbbrev1,
                'r1_StartBase__'+RefAbbrev1,
                'r2_Chromosome__'+RefAbbrev1,
                'r2_StartBase__'+RefAbbrev1,
                'CombinationCount__'+Mnemonic1,
                AbbrevHeader1])
    F4.write(TaskHeader1)
    F4.write(HeaderTranspose(OutputHeaderT4)+'\n')
    F4.write(OutputHeaderT4)
    F4.write('\n')

def PositionSortKey1(myList):
    return tuple((abs(2*pskn1-1) for pskn1 in myList))
for (c0,c2,e0,e2) in sorted(list(D0.keys()),key=PositionSortKey1):
    r0=D0[(c0,c2,e0,e2)][1]
    if c0==c2 and r0>0:
        if not(r0 in RecurrenceD1[c0]):
            RecurrenceD0[c0][r0]=0
            RecurrenceD1[c0][r0]=0
        RecurrenceD0[c0][r0]+=1
        RecurrenceD1[c0][r0]+=r0
        p20=(abs(e0)+abs(e2))//2
        pIndex = min(p20//SeparationGranularity1,len(PA1[0][c0])-1)
        sep1=-e0-e2
        PA1[11][c0][pIndex]+=r0
        PA1[12][c0][pIndex]+=1
        PA1[13][c0][pIndex]+=sep1
        PA1[14][c0][pIndex]+=sep1**2
        if 0<sep1<=Short1:
            PA1[15][c0][pIndex]+=1
        elif Short1<sep1<Long1:
            PA1[16][c0][pIndex]+=1
        elif sep1>Long1:
            PA1[17][c0][pIndex]+=1
    if CombinationList1==1 and c0==c2:
        F4.write('\t'.join(map(str,[c0,e0,c2,e2,D0[(c0,c2,e0,e2)][0]]))+'\n')
    elif CombinationList1==2 and (((c0,abs(e0)) in BaseByBaseD1) or ((c2,abs(e2)) in BaseByBaseD1)) and c0==c2:
        F4.write('\t'.join(map(str,[c0,e0,c2,e2,D0[(c0,c2,e0,e2)][0]]))+'\n')
    elif CombinationList1==3 and 7<=-e0-e2<=9 and c0==c2:
        F4.write('\t'.join(map(str,[c0,e0,c2,e2,D0[(c0,c2,e0,e2)][0]]))+'\n')
       
if CombinationList1:
    F4.close()
    LogClosingFile1(F4)

    
RecurrenceNumbers1=sorted(list(set.union(*[set(list(x.keys())) for x in RecurrenceD1])))
MultiplicityNumbers1=sorted(list(set.union(*[set(list(x.keys())) for x in MultiplicityD1])))
F5RTotals0=[sum( (RecurrenceD1[c11][i]/i for i in RecurrenceD1[c11]) ) for c11 in NameARange1]
F5RTotals1=[sum( (RecurrenceD1[c11][i] for i in RecurrenceD1[c11]) ) for c11 in NameARange1]
F5RTotals2=[sum( (RecurrenceD1[c11][i]*i for i in RecurrenceD1[c11]) ) for c11 in NameARange1]
F5RTotals0.insert(0,sum(F5RTotals0))
F5RTotals1.insert(0,sum(F5RTotals1))
F5RTotals2.insert(0,sum(F5RTotals2))
F5.write('RecurrenceSummary\tSpecies\t'+'\t'.join(map(str,F5RTotals0))+'\n')
F5.write('RecurrenceSummary\tMeanInstancesPerSpecies\t'+'\t'.join(str02(means1(F5RTotals0,F5RTotals1)))+'\n')
F5.write('RecurrenceSummary\tStandardDeviation\t'+'\t'.join(str02(sdeviations1(F5RTotals0,F5RTotals1,F5RTotals2)) )+'\n')
for i in RecurrenceNumbers1:
    row1=[]
    for c11 in NameARange1:
        if i in RecurrenceD1[c11]:
            row1.append(RecurrenceD1[c11][i])
        else:
            row1.append(0)             
    sum1=sum(row1)
    F5.write('Recurrence\t'+str(i)+'\t'+str(sum1)+'\t'+'\t'.join(map(str,row1))+'\n')
F5.write(F5Header1)
F5MTotals0=[sum(MultiplicityD1[c11].values()) for c11 in NameARange1]
F5MTotals1=[sum( (MultiplicityD1[c11][i]*i for i in MultiplicityD1[c11]) ) for c11 in NameARange1]
F5MTotals2=[sum( (MultiplicityD1[c11][i]*i*i for i in MultiplicityD1[c11]) ) for c11 in NameARange1]
F5MTotals0.insert(0,sum(F5MTotals0))
F5MTotals1.insert(0,sum(F5MTotals1))
F5MTotals2.insert(0,sum(F5MTotals2))
F5.write('MultiplicitySummary\tTotalInstances\t'+'\t'.join(map(str,F5MTotals0))+'\n')
F5.write('MultiplicitySummary\tMeanHitsPerInitial_k-mer_w_cap='+str(MaxMultiplicityReported1)+'\t'+'\t'.join(str02(means1(F5MTotals0,F5MTotals1)))+'\n')
F5.write('MultiplicitySummary\tStandardDeviation_w_cap='+str(MaxMultiplicityReported1)+'\t'+'\t'.join(str02(sdeviations1(F5MTotals0,F5MTotals1,F5MTotals2)))+'\n')
for i in MultiplicityNumbers1:
    row1=[]
    for c11 in NameARange1:
        if i in MultiplicityD1[c11]:
            row1.append(MultiplicityD1[c11][i])
        else:
            row1.append(0)             
    sum1=sum(row1)
    F5.write('Multiplicity\t'+str(i)+'\t'+str(sum1)+'\t'+'\t'.join(map(str,row1))+'\n')

F5.write(F5Header1)
F5STotals0=[sum(SeparationArray1[c11]) for c11 in NameARange1]
F5STotals1=[sum( (SeparationArray1[c11][i]*(i-Tn5DupMax1) for i in xrange(len(SeparationArray1[0])))) for c11 in NameARange1]
F5STotals2=[sum( (SeparationArray1[c11][i]*(i-Tn5DupMax1)**2 for i in xrange(len(SeparationArray1[0])))) for c11 in NameARange1]
F5STotals0.insert(0,sum(F5STotals0))
F5STotals1.insert(0,sum(F5STotals1))
F5STotals2.insert(0,sum(F5STotals2))
F5.write('SeparationSummary\tTotalInstances\t'+'\t'.join(map(str,F5STotals0))+'\n')
F5.write('SeparationSummary\tMeanSeparation_w_cap='+str(ReadSeparationMax1)+'\t'+'\t'.join(str02(means1(F5STotals0,F5STotals1)))+'\n')
F5.write('SeparationSummary\tStandardDeviation_w_cap='+str(ReadSeparationMax1)+'\t'+'\t'.join(str02(sdeviations1(F5STotals0,F5STotals1,F5STotals2)))+'\n')
for i in xrange(len(SeparationArray1[0])):
    sep1=i-Tn5DupMax1
    row1=[SeparationArray1[c11][i] for c11 in NameARange1]
    sum1=sum(row1)
    F5.write('Separation\t'+str(sep1)+'\t'+str(sum1)+'\t'+'\t'.join(map(str,row1))+'\n')
F5.write(F5Header1)
## PA1[] -  column in output - Description
## 0 c3 Non-deduplicated reads with at least one well-positioned k-mer
## 1 c4 Non-deduplicated reads mapping to focal repeats in a given region
## 2 c5 Non-deduplicated reads mapping to chromosomal repeats in a given region
## 3 c6 Non-deduplicated reads mapping to dispersed repeats in a given region
## 4 c8 Deduplicated reads with at least one well-positioned k-mer
## 5 c9 Deduplicated reads mapping to focal repeats in a given region
## 6 c10 Deduplicated reads mapping to chromosomal repeats in a given region
## 7 c11 Deduplicated reads mapping to dispersed repeats in a given region
## 8 c7 Non-deduplicated read pairs, each with at least one well-positioned k-mer
## 9 c12 Deduplicated read pairs, each with at least one well-positioned k-mer
## 10 c13 Non-deduplicated reads with at least two well-positioned k-mers
## 11 c14 Non-deduplicated read pairs, each with at least two well-positioned k-mers
## 12 c15 Deduplicated read pairs, each with at least two well-positioned k-mers
## 13 c16 Total length of unique well-positioned read pair species (means in final output)
## 14 c17 Total length**2 of unique well-positioned read pair species (std dev in final output)
## 15 c18 Total unique reads mapping to a given region with a "short" span (<=Short1
## 16 c19 Total unique reads mapping to a given region with a "medium" span (>Short1, Long1
## 17 c20 Total unique reads mapping to a given region with a "short" span (>=Long
## 18 c21 Non-deduplicated, Uniquely mapped dual-positioned circle-single-tagmentation candidates
## 19 c22 Deduplicated uniquely mapped dual-positioned circle-single-tagmentation candidates
## 20 c23 Non-deduplicated, Unique+repetative circle-single-tagmentation candidates
## 21 c24 Deduplicated, Unique+repetative circle-single-tagmentation candidates
## 22 c25 G Bases
## 23 c26 A Bases
## 24 c27 T Bases
## 25 c28 C Bases
## 26 c29 N Bases

F6Headers1='\t'.join(['Chromosome__'+RefAbbrev1,
                      'BinStart__'+RefAbbrev1,
                      'BinLength__'+RefAbbrev1,
                      'MappableSingleReads__'+Mnemonic1,
                      'FocalRepeatSingleReads__'+Mnemonic1,
                      'ChromosomeLimitedRepeatSingleReads__'+Mnemonic1,
                      'DispersedRepeatSingleReads__'+Mnemonic1,
                      'MappableSingleReadSpecies(dedup)__'+Mnemonic1,
                      'FocalRepeatSingleReadSpecies(dedup)__'+Mnemonic1,
                      'ChromosomeLimitedRepeatSingleReadSpecies(dedup)__'+Mnemonic1,
                      'DispersedRepeatSingleReadSpecies(dedup)__'+Mnemonic1,
                      'MappableReadPairs__'+Mnemonic1,
                      'MappableReadPairSpecies(dedup)__'+Mnemonic1,
                      'DMR_DualMappableReads(Dual-kmer-mapped)__'+Mnemonic1,
                      'DMRP_DualMappableReadPairs(Dual-kmer-mapped)__'+Mnemonic1,
                      'DMRP_DualMappableReadPairSpecies(Dual-kmer-mapped_dedup)__'+Mnemonic1,
                      'MeanSepDMRP__'+Mnemonic1,
                      'StdDSepDMRP__'+Mnemonic1,
                      'UniqueShortFragReadPairSpecies_Len<='+str(Short1)+'__'+Mnemonic1,
                      'UniqueMedFragReadPairSpecies_Len<'+str(Short1)+'_Len>'+str(Long1)+'__'+Mnemonic1,
                      'UniqueLongFragReadPairSpecies_Len>='+str(Long1)+'__'+Mnemonic1,
                      'UniqSingleTagmentationCandidateReads'+'__'+Mnemonic1,
                      'UniqSingleTagmentationCandidateSpecies'+'__'+Mnemonic1,
                      'AllSingleTagmentationCandidateReads'+'__'+Mnemonic1,
                      'AllSingleTagmentationCandidateSpecies'+'__'+Mnemonic1,
                      'G_Bases__'+RefAbbrev1,
                      'A_Bases__'+RefAbbrev1,
                      'T_Bases__'+RefAbbrev1,
                      'C_Bases__'+RefAbbrev1,
                      'N_Bases__'+RefAbbrev1])
if KmerBinCoverageUnique1 or KmerBinCoverageAll1:
    F6Headers1 += '\t'+'\t'.join(('UniqueKNumber__'+RefAbbrev1,
                                  'UniqueKObservedS'+'__'+Mnemonic1,
                                  'UniqueKObservedA'+'__'+Mnemonic1,
                                  'UniqueKCountS'+'__'+Mnemonic1,
                                  'UniqueKCountA'+'__'+Mnemonic1,
                                  'UniqueKCountS2'+'__'+Mnemonic1,
                                  'UniqueKCountA2'+'__'+Mnemonic1))
if KmerBinCoverageAll1:
    F6Headers1 += '\t'+'\t'.join(('LocalKNumber__'+RefAbbrev1,
                                  'LocalKMultiplicity__'+RefAbbrev1,
                                  'LocalKObservedS'+'__'+Mnemonic1,
                                  'LocalKObservedA'+'__'+Mnemonic1,
                                  'LocalKCountS'+'__'+Mnemonic1,
                                  'LocalKCountA'+'__'+Mnemonic1,
                                  'LocalKCountS2'+'__'+Mnemonic1,
                                  'LocalKCountA2'+'__'+Mnemonic1,
                                  'ChromosomalKNumber__'+RefAbbrev1,
                                  'ChromosomalKMultiplicity__'+RefAbbrev1,
                                  'ChromosomalKObservedS'+'__'+Mnemonic1,
                                  'ChromosomalKObservedA'+'__'+Mnemonic1,
                                  'ChromosomalKCountS'+'__'+Mnemonic1,
                                  'ChromosomalKCountA'+'__'+Mnemonic1,
                                  'ChromosomalKCountS2'+'__'+Mnemonic1,
                                  'ChromosomalKCountA2'+'__'+Mnemonic1,
                                  'DispersedKNumber__'+RefAbbrev1,
                                  'DispersedKMultiplicity__'+RefAbbrev1,
                                  'DispersedKObservedS'+'__'+Mnemonic1,
                                  'DispersedKObservedA'+'__'+Mnemonic1,
                                  'DispersedKCountS'+'__'+Mnemonic1,
                                  'DispersedKCountA'+'__'+Mnemonic1,
                                  'DispersedKCountS2'+'__'+Mnemonic1,
                                  'DispersedKCountA2'+'__'+Mnemonic1))
F6Headers1 += '\tUCSCLink__'+RefAbbrev1+'\t'+AbbrevHeader1
PositionArray2_M1=[[] for z in range(NumSeq1)]
PositionArray2_S1=[[] for z in range(NumSeq1)]
for c99 in NameARange1[:-1]:
    PositionArray2_M1[c99]=means1(PA1[12][c99],PA1[13][c99])
    PositionArray2_S1[c99]=sdeviations1(PA1[12][c99],PA1[13][c99],PA1[14][c99])
    for j in xrange(len(PA1[0][c99])):
        SeqStart1 = min(LD1[c99],j*SeparationGranularity1)
        SeqEnd1 = min(LD1[c99],(j+1)*SeparationGranularity1)
        CurSeq99 = SD1[c99][SeqStart1:SeqEnd1]
        PA1[22][c99][j] = CurSeq99.count('G')
        PA1[23][c99][j] = CurSeq99.count('A')
        PA1[24][c99][j] = CurSeq99.count('T')
        PA1[25][c99][j] = CurSeq99.count('C')
        PA1[26][c99][j] = CurSeq99.count('N')
        
Totals1=[]
F6.write('<!--REVA_BinByBinEventSummary-->\n')
F6.write(TaskHeader1+'\n')
F6.write(HeaderTranspose(F6Headers1)+'\n')
F6.write(F6Headers1+'\n')
for c99 in NameARange1[:-1]:
    Totals1.append(list(map(sum,(PA1[ii][c99] for ii in range(len(PA1))))))
Totals2=[]
for i in range(len(Totals1[0])):
    Totals2.append(sum( [Totals1[j][i] for j in NameARange1[:-1] ] ))
if Totals2[12]==0:
    M2='0.00'
    S2='0.00'
else:
    M2="{0:.2f}".format((Totals2[13]*1.0)/Totals2[12])
    S2="{0:.2f}".format(((Totals2[14]*1.0)/Totals2[12]-((Totals2[13]**2*1.0)/Totals2[12]**2))**0.5)
Totals2[13]=M2
Totals2[14]=S2

F6.write('All\t1\t'+str(sum(LD1))+'\t'+'\t'.join(map(str,Totals2))+'\n')


for c99 in NameARange1[:-1]:
    nanolist0=Totals1[c99][:]
    if nanolist0[12]==0:
        M0='0.00'
        S0='0.00'
    else:
        M0="{0:.2f}".format((nanolist0[13]*1.0)/nanolist0[12])
        S0="{0:.2f}".format(((nanolist0[14]*1.0)/nanolist0[12]-((nanolist0[13]**2*1.0)/nanolist0[12]**2))**0.5)
    nanolist0[13]=M0
    nanolist0[14]=S0
    F6.write(NameA1[c99]+'\t1\t'+str(LD1[c99])+'\t'+'\t'.join(map(str,nanolist0))+'\n')         
    PA1[13][c99]=str02(PositionArray2_M1[c99])
    PA1[14][c99]=str02(PositionArray2_S1[c99])

for c99 in NameARange1[:-1]:
    i99=0
    for nanolist1 in zip(*[PA1[i][c99] for i in range(len(PA1))]):        
        if LD1[c99]>=i99+SeparationGranularity1:
            BinLength1=SeparationGranularity1
        else:
            BinLength1=LD1[c99]-i99
        UCSC1=UCSClink1(c99,i99+1,i99+BinLength1)
        F6.write(NameA1[c99]+'\t'+str(i99+1)+'\t'+str(BinLength1)+'\t'+'\t'.join(map(str,nanolist1))+'\t'+UCSC1+'\n')
        i99+=SeparationGranularity1
OutputHeaders1 = []
## Make generalized header list for bin-by-bin or feature-by feature output summaries
if ReportUnique1:
    OutputHeaders1.append('Reference_SingleCopyKMerCount___'+RefAbbrev1)
if ReportRepeats1:
    if ReportSense1:
        OutputHeaders1.append('Reference_LocalRepeatKMerCount___'+RefAbbrev1)
    if ReportSense1:
        OutputHeaders1.append('Reference_ChromosomalRepeatKMerCount___'+RefAbbrev1)
    if ReportSense1:
        OutputHeaders1.append('Reference_DispersedRepeatKMerCount___'+RefAbbrev1)
MultiplicityTypes1 = []
if ReportUnique1:
    MultiplicityTypes1.append('SingleCopyKmers')
if ReportRepeats1:
    MultiplicityTypes1.append('LocalRepeatKmers')
if ReportRepeats1:
    MultiplicityTypes1.append('ChromosomalRepeatKmers')
if ReportRepeats1:
    MultiplicityTypes1.append('DispersedRepeatKmers')
EventTypes1 = []
if ReportCoverage1:
    EventTypes1.append('Covering')
if ReportStarts1:
    EventTypes1.append('Start')
if ReportEnds1:
    EventTypes1.append('End')
Strands1 = []
if ReportSense1:
    Strands1.append('Sense')
if ReportAntisense1:
    Strands1.append('Antisense')
if ReportBothStrands1:
    Strands1.append('')
DuplicityTypes1 = []
if ReportPositions1:
    DuplicityTypes1.append('Positions')
if ReportReads1:
    DuplicityTypes1.append('Counts')
for ri3 in MultiplicityTypes1:
    for ri1 in EventTypes1:
        for ri4 in Strands1:
            for ri2 in DuplicityTypes1:
                CurHead1 = ''.join((ri4,ri1,ri2,'_',ri3))
                CurHead1 = CurHead1.replace("CoveringPositions","CoveredPositions")
                CurHead1 = CurHead1.replace("CoveringCounts","SummedKMerCoverage")
                if PreREVA1:
                    CurHead1 = CurHead1+"_PreREVA"
                OutputHeaders1.append(CurHead1+'__'+Mnemonic1)
if UCSCLinks1:
    OutputHeaders1.append('UCSC_Link__'+RefAbbrev1)
OutputHeaders1 = '\t'.join(OutputHeaders1)+'\t'+AbbrevHeader1
if GFF1:
    F7Buff1 = []
    Dft1 = [] # list of feature types
    Dfl1 = [] # list of stripped lines in GFF file(s)
    Dfc1 = [] # list of ordinal chromosome numbers for gffs
    Dfs1 = [] # list of GFF feature starts
    Dfe1 = [] # list of GFF feature ends
    Dfo1 = [] # orientation of GFF feature (+/-/.)  ## not used for now but could be used
    HeadWritten1 = False
    for GFF11 in GFF1:
        LogOpeningFile1(GFF11)
        if GFF11.endswith('gz'):
            if version.startswith('2.'):
                GFFile1=gzip.open(GFF11,mode='r')
            else:
                GFFile1=gzip.open(GFF11,mode='rt')
        else:
            if version.startswith('2.'):
                GFFile1=open(GFF11,mode='rU')
            else:
                GFFile1=open(GFF11,mode='r')
        for L0 in GFFile1:
            L1=L0.strip().split('\t')
            if len(L1)<5 and not(HeadWritten1):
                F7Buff1.append(L0.strip())
                continue
            if not(HeadWritten1) and not(L1[3].isdigit() and L1[4].isdigit()):
                F7Buff1.append(L0.strip()+'\t'+OutputHeaders1+'\t'+AbbrevHeader1+'\n')
                HeadWritten1 = True
                continue
            if not(HeadWritten1):
                GFFHead0 = ['Seqname','Source','Feature','Start','End','Score','Strand','Frame','Attribute']
                if len(L1)>len(GFFHead0):
                    GFFHead0.extend(['GFF_ColumnUnlabeled']*(len(L1)-len(GFFHead0)))
                F7Buff1.append('\t'.join([GFFt1+'__'+GFF11 for GFFt1 in GFFHead0[:len(L1)]])+'\t'+OutputHeaders1+'\t'+AbbrevHeader1+'\n')
                HeadWritten1 = True
            ch = L1[0]
            if ch.lower().startswith('m') and not(ch in NameD1) and not('chr'+ch in NameD1):  ## very cumbersome temporary code trying to deal with different ways of naming MtDNA
                if ('M' in NameD1) or ('chrM' in NameD1):
                    ch = 'M'
            if not ch in NameD1:
                ch = 'chr'+ch
            if (ch in NameD1) and ((L1[2] in Feature1D) or Feature1All):
                if FeatureTags1:
                    Keep1 = False
                    for ft1 in FeatureTags1:
                        if ft1 in L0.lower():
                            Keep1 = True
                            break
                    if not(Keep1):
                        continue                       
                Dfc1.append(NameD1[ch])
                Dft1.append(L1[2])
                Dfs1.append(int(L1[3]))
                Dfe1.append(int(L1[4]))
                Dfl1.append(L0.strip())
                Dfo1.append(L1[6])
        GFFile1.close()
    F7=open("FeatureSummary_"+OutFileNameBase+'.tdv',mode='w')
    LogOpeningFile1(F7)
    F7.write('<!--REVA_FeatureByFeatureCounts-->\n')
    F7.write(TaskHeader1+'\n')
    F7.write(HeaderTranspose(F7Buff1[-1])+'\n')
    F7.write('\n'.join(F7Buff1))

if BaseByBase1:
    F9 = open(BaseByBaseOutFile1,mode='w')
    LogOpeningFile1(F9)
    F9.write('<!--REVA_BaseByBaseCounts-->\n')
    F9.write(TaskHeader1+'\n')
    F9.write(HeaderTranspose(BaseByBaseHeader1)+'\n')
    F9.write(BaseByBaseHeader1+'\t'+AbbrevHeader1+'\n')


if BinByBin1:
    BinHead1  = '\t'.join(['Chromosome__'+RefAbbrev1,
                          'BinStart__'+RefAbbrev1,
                          'BinLength__'+RefAbbrev1])+'\t'+OutputHeaders1

    F10=open("BinByBinReadCountSummary_"+OutFileNameBase+'.tdv',mode='w')
    LogOpeningFile1(F10)
    F10.write('<!--REVA_BinByBinReadCountSummary-->\n')
    F10.write(TaskHeader1+'\n')
    F10.write(HeaderTranspose(BinHead1)+'\n')
    F10.write(BinHead1+'\n')

if ChromosomeByChromosome1:
    ChromosomeHead1  = '\t'.join(['Chromosome__'+RefAbbrev1,
                          'BinStart__'+RefAbbrev1,
                          'BinLength__'+RefAbbrev1])+'\t'+OutputHeaders1

    F11=open("ChromosomeByChromosomeReadCountSummary_"+OutFileNameBase+'.tdv',mode='w')
    LogOpeningFile1(F11)
    F11.write('<!--REVA_ChromosomeByChromosomeReadCountSummary-->\n')
    F11.write(TaskHeader1+'\n')
    F11.write(HeaderTranspose(ChromosomeHead1)+'\n')
    F11.write(BinHead1+'\n')


if GFF1 or BaseByBase1 or BinByBin1 or ChromosomeByChromosome1:
    for z in range(NumSeq1):        
        PA7s = []
        PA7a = []
        PA7k = []
        allchr = np.where(Sort_I1dedup==z)[0]
        allchrS = allchr[Sort_P1dedup[allchr]>0]
        allchrA = allchr[Sort_P1dedup[allchr]<0]
        allqS = (Sort_P1dedup[allchrS]-1) % LD1[z]
        allqA = (-Sort_P1dedup[allchrA]-klen1) % LD1[z]
        MultS = np.zeros(LD1[z],Multiplicities1.dtype)
        MultA = np.zeros(LD1[z],Multiplicities1.dtype)
        MultS[allqS] = Multiplicities1[allchrS]
        MultA[allqA] = Multiplicities1[allchrA]
        CovKzS = np.zeros(LD1[z],CoverageK1.dtype)
        CovKzA = np.zeros(LD1[z],CoverageK1.dtype)
        CovKzS[allqS] = CoverageK1[allchrS]
        CovKzA[allqA] = CoverageK1[allchrA]
        if ReportStarts1:
            StartzS = np.zeros(LD1[z],CoverageS1.dtype)
            StartzA = np.zeros(LD1[z],CoverageS1.dtype)
            StartzS[allqS] = CoverageS1[allchrS]
            StartzA[allqA] = CoverageS1[allchrA]
        if ReportEnds1:
            EndzS = np.zeros(LD1[z],CoverageE1.dtype)
            EndzA = np.zeros(LD1[z],CoverageE1.dtype)
            EndzS[allqS] = CoverageE1[allchrS]
            EndzA[allqA] = CoverageE1[allchrA]
        if BaseByBase1:
            for bbi1 in xrange(LD1[z]):
                if (z,bbi1) in BaseByBaseD1:
                    BaseByBaseLine1 = []
                    bbN1,bbP1 = BaseByBaseD1[(z,bbi1)]
                    if bbbC1:
                        BaseByBaseLine1.append(NameA1[z])
                    if bbbP1:
                        BaseByBaseLine1.append(str(bbi1+1))
                    if bbbF1:
                        BaseByBaseLine1.append(str(BfN1[bbN1]))
                    if bbbO1:
                        BaseByBaseLine1.append(str(bbP1+1))
                    if bbbB1:
                        BaseByBaseLine1.append(SD1[z][bbi1])
                    if bbbR1:
                        bbm1 = MultS[(bbi1-klen2)%LD1[z]]
                        if bbm1==1:
                            bbr1 = 'U'
                        elif bbm1 & Local1:
                            bbr1 = 'L'
                        elif bbm1 & Chromosomal1:
                            bbr1 = 'C'
                        elif bbm1>1:
                            bbr1 = 'D'
                        else: 
                            bbr1 = 'E'                          
                        BaseByBaseLine1.append(bbr1)
                    if bbbM1:
                        BaseByBaseLine1.append(str(Mult1&MultS[(bbi1-klen2)%LD1[z]]))
                    if bbbK1:
                        BaseByBaseLine1.append(str(CovKzS[(bbi1-klen2)%LD1[z]]))
                        BaseByBaseLine1.append(str(CovKzA[(bbi1-klen2)%LD1[z]]))
                    if bbbS1:
                        BaseByBaseLine1.append(str(StartzS[bbi1]))
                        BaseByBaseLine1.append(str(StartzA[(bbi1-klen1+1) % LD1[z]]))
                    if bbbE1:
                        BaseByBaseLine1.append(str(EndzS[(bbi1-klen1+1) % LD1[z]]))
                        BaseByBaseLine1.append(str(EndzA[bbi1]))
                    if bbbT1:
                        if (z,bbi1+1) in ExtensionD1:
                            BaseByBaseLine1.append(str(ExtensionD1[(z,bbi1+1)]).replace("'",""))
                        else:
                            BaseByBaseLine1.append('')
                    BaseByBaseLine1 = '\t'.join(BaseByBaseLine1)
                    F9.write(BaseByBaseLine1+'\n')            
        if ReportUnique1:
            UniqueS1 = np.where( MultS==1 )[0]
            UniqueA1 = np.where( MultA==1 )[0]
            KzSU = np.zeros(LD1[z],np.uint8)
            KzSU[(UniqueS1+klen2)%LD1[z]] = 1
            PA7k.append(KzSU)
            if ReportCoverage1:
                CovzSU = np.zeros(LD1[z],CoverageK1.dtype)
                CovzAU = np.zeros(LD1[z],CoverageK1.dtype)
                CovzSU[(UniqueS1+klen2)%LD1[z]] = CovKzS[UniqueS1]
                CovzAU[(UniqueA1+klen2)%LD1[z]] = CovKzA[UniqueA1]
                PA7s.append(CovzSU)
                PA7a.append(CovzAU)
            if ReportStarts1:
                StartzSU = np.zeros(LD1[z],CoverageS1.dtype)
                StartzAU = np.zeros(LD1[z],CoverageS1.dtype)
                StartzSU[UniqueS1] = StartzS[UniqueS1]
                StartzAU[(UniqueA1+klen1-1)%LD1[z]] = StartzA[UniqueA1]
                PA7s.append(StartzSU)
                PA7a.append(StartzAU)
            if ReportEnds1:
                EndzSU = np.zeros(LD1[z],CoverageS1.dtype)
                EndzAU = np.zeros(LD1[z],CoverageS1.dtype)
                EndzSU[UniqueS1] = EndzS[UniqueS1]
                EndzAU[(UniqueA1+klen1-1)%LD1[z]] = EndzA[UniqueA1]
                PA7s.append(EndzSU)
                PA7a.append(EndzAU)
        if ReportRepeats1:
            LocalS1 = np.where( (MultS & Local1)>0 )[0]
            LocalA1 = np.where( (MultA & Local1)>0 )[0]
            KzSL = np.zeros(LD1[z],np.uint8)
            KzSL[(LocalS1+klen2)%LD1[z]] = 1
            PA7k.append(KzSL)
            if ReportCoverage1:
                CovzSL = np.zeros(LD1[z],CoverageK1.dtype)
                CovzAL = np.zeros(LD1[z],CoverageK1.dtype)
                CovzSL[(LocalS1+klen2)%LD1[z]] = CovKzS[LocalS1]
                CovzAL[(LocalA1+klen2)%LD1[z]] = CovKzA[LocalA1]
                PA7s.append(CovzSL)
                PA7a.append(CovzAL)
            if ReportStarts1:
                StartzSL = np.zeros(LD1[z],CoverageS1.dtype)
                StartzAL = np.zeros(LD1[z],CoverageS1.dtype)
                StartzSL[LocalS1] = StartzS[LocalS1]
                StartzAL[(LocalA1+klen1-1)%LD1[z]] = StartzA[LocalA1]
                PA7s.append(StartzSL)
                PA7a.append(StartzAL)
            if ReportEnds1:
                EndzSL = np.zeros(LD1[z],CoverageS1.dtype)
                EndzAL = np.zeros(LD1[z],CoverageS1.dtype)
                EndzSL[LocalS1] = EndzS[LocalS1]
                EndzAL[(LocalA1+klen1-1)%LD1[z]] = EndzA[LocalA1]
                PA7s.append(EndzSL)
                PA7a.append(EndzAL)
            ChromosomalS1 = np.where( (MultS & Chromosomal1)>0 )[0]
            ChromosomalA1 = np.where( (MultA & Chromosomal1)>0 )[0]
            KzSC = np.zeros(LD1[z],np.uint8)
            KzSC[(ChromosomalS1+klen2)%LD1[z]] = 1
            PA7k.append(KzSC)
            if ReportCoverage1:
                CovzSC = np.zeros(LD1[z],CoverageK1.dtype)
                CovzAC = np.zeros(LD1[z],CoverageK1.dtype)
                CovzSC[(ChromosomalS1+klen2)%LD1[z]] = CovKzS[ChromosomalS1]
                CovzAC[(ChromosomalA1+klen2)%LD1[z]] = CovKzA[ChromosomalA1]
                PA7s.append(CovzSC)
                PA7a.append(CovzAC)
            if ReportStarts1:
                StartzSC = np.zeros(LD1[z],CoverageS1.dtype)
                StartzAC = np.zeros(LD1[z],CoverageS1.dtype)
                StartzSC[ChromosomalS1] = StartzS[ChromosomalS1]
                StartzAC[(ChromosomalA1+klen1-1)%LD1[z]] = StartzA[ChromosomalA1]
                PA7s.append(StartzSC)
                PA7a.append(StartzAC)
            if ReportEnds1:
                EndzSC = np.zeros(LD1[z],CoverageS1.dtype)
                EndzAC = np.zeros(LD1[z],CoverageS1.dtype)
                EndzSC[ChromosomalS1] = EndzS[ChromosomalS1]
                EndzAC[(ChromosomalA1+klen1-1)%LD1[z]] = EndzA[ChromosomalA1]
                PA7s.append(EndzSC)
                PA7a.append(EndzAC)
            DispersedS1 = np.where( (MultS < Local1) & (MultS > 1) )[0]
            DispersedA1 = np.where( (MultA < Local1) & (MultA > 1) )[0]
            KzSD = np.zeros(LD1[z],np.uint8)
            KzSD[(DispersedS1+klen2)%LD1[z]] = 1
            PA7k.append(KzSD)
            if ReportCoverage1:
                CovzSD = np.zeros(LD1[z],CoverageK1.dtype)
                CovzAD = np.zeros(LD1[z],CoverageK1.dtype)
                CovzSD[(DispersedS1+klen2)%LD1[z]] = CovKzS[DispersedS1]
                CovzAD[(DispersedA1+klen2)%LD1[z]] = CovKzA[DispersedA1]
                PA7s.append(CovzSD)
                PA7a.append(CovzAD)
            if ReportStarts1:
                StartzSD = np.zeros(LD1[z],CoverageS1.dtype)
                StartzAD = np.zeros(LD1[z],CoverageS1.dtype)
                StartzSD[DispersedS1] = StartzS[DispersedS1]
                StartzAD[(DispersedA1+klen1-1)%LD1[z]] = StartzA[DispersedA1]
                PA7s.append(StartzSD)
                PA7a.append(StartzAD)
            if ReportEnds1:
                EndzSD = np.zeros(LD1[z],CoverageS1.dtype)
                EndzAD = np.zeros(LD1[z],CoverageS1.dtype)
                EndzSD[DispersedS1] = EndzS[DispersedS1]
                EndzAD[(DispersedA1+klen1-1)%LD1[z]] = EndzA[DispersedA1]
                PA7s.append(EndzSD)
                PA7a.append(EndzAD)
        if GFF1:
            for t7,l7,c7,s7,e7,o7 in izip(Dft1,Dfl1,Dfc1,Dfs1,Dfe1,Dfo1):
                if c7==z:
                    if s7<=0:
                        s7=1
                    if e7>LD1[z]:
                        e7=LD1[z]
                    F7Buffer1= [l7,]
                    for pak7 in PA7k:
                        F7Buffer1.append(np.sum(pak7[s7-1:e7]))
                    for pa7s,pa7a in izip(PA7s,PA7a):
                        if o7=='-':
                            pa7s,pa7a = pa7a,pa7s
                        if ReportSense1:
                            if ReportPositions1:
                                F7Buffer1.append(np.count_nonzero(pa7s[s7-1:e7]))
                            if ReportReads1:
                                F7Buffer1.append(np.sum(pa7s[s7-1:e7]))
                        if ReportAntisense1:
                            if ReportPositions1:
                                F7Buffer1.append(np.count_nonzero(pa7a[s7-1:e7]))
                            if ReportReads1:
                                F7Buffer1.append(np.sum(pa7a[s7-1:e7]))
                        if ReportBothStrands1:
                            if ReportPositions1:
                                F7Buffer1.append(np.count_nonzero(pa7a[s7-1:e7])+np.count_nonzero(pa7s[s7-1:e7]))
                            if ReportReads1:
                                F7Buffer1.append(np.sum(pa7a[s7-1:e7])+np.sum(pa7s[s7-1:e7]))
                    if UCSCLinks1:
                        try:
                            F7Buffer1.append(UCSClink1(c7,s7-UCSCBuffer1,e7+UCSCBuffer1))
                        except:
                            pass
                    ##Wormbaselink1(n7),
                    F7.write('\t'.join(map(str,F7Buffer1))+'\n')
        if BinByBin1:
            for bi1 in xrange(0,LD1[z],SeparationGranularity1):
                end1 = min(bi1+SeparationGranularity1,LD1[z])
                F10Buffer1= [NameA1[z],]
                F10Buffer1.append(bi1+1)
                F10Buffer1.append(end1-bi1)
                for pak10 in PA7k:
                    F10Buffer1.append(np.sum(pak10[bi1:end1]))
                for pa7s,pa7a in izip(PA7s,PA7a):
                    if ReportSense1:
                        if ReportPositions1:
                            F10Buffer1.append(np.count_nonzero(pa7s[bi1:end1]))
                        if ReportReads1:
                            F10Buffer1.append(np.sum(pa7s[bi1:end1]))
                    if ReportAntisense1:
                        if ReportPositions1:
                            F10Buffer1.append(np.count_nonzero(pa7a[bi1:end1]))
                        if ReportReads1:
                            F10Buffer1.append(np.sum(pa7a[bi1:end1]))
                    if ReportBothStrands1:
                        if ReportPositions1:
                            F10Buffer1.append(np.count_nonzero(pa7a[bi1:end1])+np.count_nonzero(pa7s[bi1:end1]))
                        if ReportReads1:
                            F10Buffer1.append(np.sum(pa7a[bi1:end1])+np.sum(pa7s[bi1:end1]))
                F10.write('\t'.join(map(str,F10Buffer1))+'\n')
        if ChromosomeByChromosome1:
            F11Buffer1= [NameA1[z],]
            F11Buffer1.append(1)
            F11Buffer1.append(LD1[z])
            for pak11 in PA7k:
                F11Buffer1.append(np.sum(pak11))
            for pa7s,pa7a in izip(PA7s,PA7a):
                if ReportSense1:
                    if ReportPositions1:
                        F11Buffer1.append(np.count_nonzero(pa7s))
                    if ReportReads1:
                        F11Buffer1.append(np.sum(pa7s))
                if ReportAntisense1:
                    if ReportPositions1:
                        F11Buffer1.append(np.count_nonzero(pa7a))
                    if ReportReads1:
                        F11Buffer1.append(np.sum(pa7a))
                if ReportBothStrands1:
                    if ReportPositions1:
                        F11Buffer1.append(np.count_nonzero(pa7a)+np.count_nonzero(pa7s))
                    if ReportReads1:
                        F11Buffer1.append(np.sum(pa7a)+np.sum(pa7s))
            F11.write('\t'.join(map(str,F11Buffer1))+'\n')
    if GFF1:
        F7.close()
        LogClosingFile1(F7)
    if BaseByBase1:
        F9.close()
        LogClosingFile1(F9)
    if BinByBin1:
        F10.close()
        LogClosingFile1(F10)
    if ChromosomeByChromosome1:
        F11.close()
        LogClosingFile1(F11)
        


if DFAMCount1:
    F8=open("DFAM_RepeatSummary_"+OutFileNameBase+'.tdv',mode='w')
    LogOpeningFile1(F8)
    F8Headers1='\t'.join(['REVATableDFAM_Repeat_Name__'+DFAM1,
                      'Repeat_Length__'+DFAM1,
                      'FeatureID__'+DFAM1,                         
                      'Repeat_Description__'+DFAM1,
                      'DFAMLink__'+DFAM1,
                      'Hit_Count__'+Mnemonic1,
                      'PercentTotalReads__'+Mnemonic1,
                      'RPKM__'+Mnemonic1,
                      'Aligned_Bases__'+Mnemonic1,
                      'Identified_Instances__'+Mnemonic1])
    F8.write('<!--REVA_DFAM_Repeat_Summary-->\n')
    F8.write(TaskHeader1+'\n')
    F3.write(HeaderTranspose(F8Headers1)+'\n')
    F8.write(F8Headers1+'\t'+AbbrevHeader1+'\n')
    for j in xrange(len(DfamX1)):
        rpkm1=(1000000000.0 * DfamC1[j]) / (2*TotalReadPairs1 * DfamB1[j])
        ptr1=(100.0 * DfamC1[j]) / (2*TotalReadPairs1)
        F8.write('\t'.join((DfamX1[j],
                               str(DfamL1[j]),
                               DfamW1[j],
                               DfamD1[j],
                               DfamW2[j],
                               str(DfamC1[j]),
                               '{0:.4f}'.format(ptr1),
                               '{0:.2f}'.format(rpkm1),
                               str(DfamB1[j]),
                               str(DfamT1[j])+'\n')))        
    F8.close()
    LogClosingFile1(F8)


if FindStructuralAnomalies1:
    F3.close()
    LogClosingFile1(F3)
F5.close(); F6.close() 
LogClosingFile1(F5); LogClosingFile1(F6)         


if FirstKTable1:
    J0As = float(max(1,sum(J0A)))
    J1As = float(max(1,sum(J1A)))
    J2As = float(max(1,sum(J2A)))
    J3As = float(max(1,sum(J3A)))
    if FirstKByStrand1:
        JS0As = float(max(1,sum(JS0A)))
        JS1As = float(max(1,sum(JS1A)))
        JS2As = float(max(1,sum(JS2A)))
        JS3As = float(max(1,sum(JS3A)))
        JA0As = float(max(1,sum(JA0A)))
        JA1As = float(max(1,sum(JA1A)))
        JA2As = float(max(1,sum(JA2A)))
        JA3As = float(max(1,sum(JA3A)))
    FirstKHead1 = 'FirstKOffset__'+RefAbbrev1
    if SingleReadMode1==0 or (J0As>0 and J2As>0):
        AveNote1 = '    Average_FirstKValues_Cumulative: StartR1=' + '{0:.2f}'.format(Moment1(J0A))+\
                                     ', EndR1='   + '{0:.2f}'.format(Moment1(J1A))+\
                                     ', StartR2=' + '{0:.2f}'.format(Moment1(J2A))+\
                                     ', EndR2='   + '{0:.2f}'.format(Moment1(J3A))
        FirstKHead1 += '\tR1Start_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR1End_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR2Start_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR2End_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR1Start_FirstKOffsetFraction__'+RefAbbrev1
        FirstKHead1 += '\tR1End_FirstKOffsetFraction__'+RefAbbrev1
        FirstKHead1 += '\tR2Start_FirstKOffsetFraction__'+RefAbbrev1
        FirstKHead1 += '\tR2End_FirstKOffsetFraction__'+RefAbbrev1
        if FirstKByStrand1:
            FirstKHead1 += '\tR1StartSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1EndSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2StartSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2EndSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1StartSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1EndSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2StartSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2EndSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1StartAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1EndAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2StartAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2EndAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1StartAntiSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1EndAntiSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2StartAntiSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2EndAntiSense_FirstKOffsetFraction__'+RefAbbrev1
    elif SingleReadMode1==1 or J0As>0:
        AveNote1 = '    Average_FirstKValues_Cumulative: StartR1=' + '{0:.2f}'.format(Moment1(J0A))+\
                                     ', EndR1='   + '{0:.2f}'.format(Moment1(J1A)) 
        FirstKHead1 += '\tR1Start_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR1End_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR1Start_FirstKOffsetFraction__'+RefAbbrev1
        FirstKHead1 += '\tR1End_FirstKOffsetFraction__'+RefAbbrev1
        if FirstKByStrand1:
            FirstKHead1 += '\tR1StartSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1EndSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1StartSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1EndSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1StartAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1EndAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1StartAntiSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1EndAntiSense_FirstKOffsetFraction__'+RefAbbrev1
    elif SingleReadMode1==2 or J2As>0:
        AveNote1 = '    Average_FirstKValues_Cumulative: StartR2=' + '{0:.2f}'.format(Moment1(J2A))+\
                                     ', EndR2='   + '{0:.2f}'.format(Moment1(J3A))
        FirstKHead1+= '\tR2Start_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1+= '\tR2End_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR2Start_FirstKOffsetFraction__'+RefAbbrev1
        FirstKHead1 += '\tR2End_FirstKOffsetFraction__'+RefAbbrev1
        if FirstKByStrand1:
            FirstKHead1 += '\tR2StartSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2EndSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2StartSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2EndSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2StartAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2EndAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2StartAntiSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2EndAntiSense_FirstKOffsetFraction__'+RefAbbrev1
    CapNote1 = "Note that all offset values above k-len cutoff (in this case cutoff=" +str(len(J0A)-1)+ ") have been assigned offset value at cutoff"    
    LogNote1(AveNote1,LogFile1)
    LogNote1(CapNote1,LogFile1)
    F13=open("FirstKTable_"+OutFileNameBase+'.tdv',mode='w')
    LogOpeningFile1(F13)
    F13.write('<!--REVA_FirstMatchedKMerSummary-->\n')
    F13.write(TaskHeader1+'\n')
    FirstKHead1 += '\t'+AbbrevHeader1
    F13.write(HeaderTranspose(FirstKHead1)+'\n')
    F13.write('<!--'+AveNote1+'\n')
    F13.write('<!--'+CapNote1+'\n')
    F13.write('\n')
    F13.write(FirstKHead1+'\n')
    for i in range(len(J0A)):
        F13.write(str(i)+'\t')
        if SingleReadMode1==0 or (J0As>0 and J2As>0):
            F13.write(str(J0A[i])+'\t')
            F13.write(str(J1A[i])+'\t')
            F13.write(str(J2A[i])+'\t')
            F13.write(str(J3A[i])+'\t')
            F13.write('{0:.6f}'.format(J0A[i]/J0As)+'\t')
            F13.write('{0:.6f}'.format(J1A[i]/J1As)+'\t')
            F13.write('{0:.6f}'.format(J2A[i]/J2As)+'\t')
            F13.write('{0:.6f}'.format(J3A[i]/J3As)+'\t')
            if FirstKByStrand1:
                F13.write(str(JS0A[i])+'\t')
                F13.write(str(JS1A[i])+'\t')
                F13.write(str(JS2A[i])+'\t')
                F13.write(str(JS3A[i])+'\t')
                F13.write('{0:.6f}'.format(JS0A[i]/JS0As)+'\t')
                F13.write('{0:.6f}'.format(JS1A[i]/JS1As)+'\t')
                F13.write('{0:.6f}'.format(JS2A[i]/JS2As)+'\t')
                F13.write('{0:.6f}'.format(JS3A[i]/JS3As)+'\t')
                F13.write(str(JA0A[i])+'\t')
                F13.write(str(JA1A[i])+'\t')
                F13.write(str(JA2A[i])+'\t')
                F13.write(str(JA3A[i])+'\t')
                F13.write('{0:.6f}'.format(JA0A[i]/JA0As)+'\t')
                F13.write('{0:.6f}'.format(JA1A[i]/JA1As)+'\t')
                F13.write('{0:.6f}'.format(JA2A[i]/JA2As)+'\t')
                F13.write('{0:.6f}'.format(JA3A[i]/JA3As)+'\t')

        elif SingleReadMode1==1 or J0As>0:
            F13.write(str(J0A[i])+'\t')
            F13.write(str(J1A[i])+'\t')
            F13.write('{0:.6f}'.format(J0A[i]/J0As)+'\t')
            F13.write('{0:.6f}'.format(J1A[i]/J1As)+'\t')
            if FirstKByStrand1:
                F13.write(str(JS0A[i])+'\t')
                F13.write(str(JS1A[i])+'\t')
                F13.write('{0:.6f}'.format(JS0A[i]/JS0As)+'\t')
                F13.write('{0:.6f}'.format(JS1A[i]/JS1As)+'\t')
                F13.write(str(JA0A[i])+'\t')
                F13.write(str(JA1A[i])+'\t')
                F13.write('{0:.6f}'.format(JA0A[i]/JA0As)+'\t')
                F13.write('{0:.6f}'.format(JA1A[i]/JA1As)+'\t')
        elif SingleReadMode1==2 or J2As>0:
            F13.write(str(J2A[i])+'\t')
            F13.write(str(J3A[i])+'\t')
            F13.write('{0:.6f}'.format(J2A[i]/J2As)+'\t')
            F13.write('{0:.6f}'.format(J3A[i]/J3As)+'\t')
            if FirstKByStrand1:
                F13.write(str(JS2A[i])+'\t')
                F13.write(str(JS3A[i])+'\t')
                F13.write('{0:.6f}'.format(JS2A[i]/JS2As)+'\t')
                F13.write('{0:.6f}'.format(JS3A[i]/JS3As)+'\t')
                F13.write(str(JA2A[i])+'\t')
                F13.write(str(JA3A[i])+'\t')
                F13.write('{0:.6f}'.format(JA2A[i]/JA2As)+'\t')
                F13.write('{0:.6f}'.format(JA3A[i]/JA3As)+'\t')
        F13.write('\n')
    F13.close()
    LogClosingFile1(F13)     


if 'darwin' in platform:
    try:
        Coffee_process1.kill()
    except:
        LogNote1("Couldn't kill 'caffeinate' process, you may need to manually reset your energy saving preferences (System Preferences, Energy)",LogFile1)
PurgeREVATempFiles(ScratchPath1)
LogNote1('Finished running '+' '.join(argv),LogFile1)
try:
    LogNote1(open(MyCode1,mode='r').read(),LogFile1)
except:
    pass


