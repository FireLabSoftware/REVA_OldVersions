#!/usr/bin/env python -i

## REVA Output JunctionParser version xy11 01-21-20

## Goals: Parse possible split-read junctions from Illumina data and identify
## most likely candidates for auithentic structural variants

##        Single line display of Junction read versus reference sequence regions
##        Here s1 is the junction read, s2 and s3 are the inferred left and right ends of the reference 
##            A a single parsed string that describes the alignment
##                # mismatches in s2 and s3 relative to s1
##                M,m mismatches relative to s1 in only s2 or only s3 respectively 
##                v Deletion in s2 and s3 relative to s1
##                D,d deletions relative to s1 in only s2 or only s3 respectively 
##                + Insertion in s2 and s3 relative to s1
##                I,i insertions relative to s1 in only s2 or only s3 respectively 
##                () possible region of junction
##                - both sequenes match
##                * different variants on the two strands
##                . padding base-- terminal sequence not part of alignment (optional... no effect on quality)

##    Public Parameters:
##        s1 = sequence 1 (also an __init__ parameter)
##        s2 = sequence 2 (not reversed for ReverseAlign) (also an __init__ parameter)
##        ReverseAlign = True if the sequences are aligned end first (also an __init__ parameter)
##        TwoStep = True if alignment is to be done in two steps-- first finding the best alignment score and then re-choosing all downstream points to obtain best alignments that go through the maximum
##        MakeCigA = True if we're making an alignment-describing string called cigA (see below)
##        span = The maximum slippage allowed between the two sequences (also an __init__ parameter, default = 12)
##        l1 = len(s1)
##        l2 = len(s2)
##        s1Max = Position in one-based Sequence s1 corresponding to the maximum alignment score
##        s2Max = Position in one-based Sequence s2 corresponding to the maximum alignment score
##        vMax = Best alignment score
##        mV = array of best alignment values for each position in s1.  mV[0] is zero and mv(x) is the best alignment score for s1 position x-1
##        rV = s2position-s1position offset for the end of best alignment at position x+1 in s1.  Like mV, rV is offset by one, so rV[x] is the s2-s1 offset at position x-1
##        vB = Best alignment score for whole sequence... such alignments allow trailing or leading bases to be cut off with no penalty, but otherwise require the entire participation of s1 and s2 in the alignment
##        s1B = position in one-based s1 at the end of the best whole sequence alignment ( this will be len(s1) unless there are trailing bases from s1 once aligned with s2)
##        s2B = position in one-based s2 at the end of the best whole sequence alignment ( this will be len(s2) unless there are trailing bases from s1 once aligned with s2)

## cigA language
    ## M is a match between string 1 and string 2
    ## m is a mismatch between string 1 and string 2
    ## d is a base that is present in string 2 but missing in string 1
    ## D is a base that is present in string 1 but missing in string 2
    ## e is a terminal "filler" to be added to string 1 to allow alignment
    ## E is a terminal "filler" to be added to string 2 to allow alignment

## Fenceposting
## For a sequence S and Base B, let's say the ordinal (starting with 1) position of B is n
## then the index of B   in the sequence array S is n-1
## In antisense(S), Antisense(B) is at index (Len(S)-1)-n-1 = len(S)-n
## In antisense array, antisense(S), Antisense(B) is at len(S)-n-1
## d1 and e1 arrays are one-based, so the value of aligning positions 1..n in a sequence is d1[n]
## The mode of addition for base B onto the alignment (match/mismatch, deleted, insertion after) is e1[n]
## A perfectly-matched alignment to base n should give a score and mismatch count of n.
## For a forward alignment, mV[n] is the best alignment up to base n.  mV[0] is alignment to base zero and by definition zero
## rV[n] is the relative position in sequence S2 for alignment to base n of sequence S1.  So rV[n]==0 means that alignment to base B gives the highest score if base B is aligned with the n'th base of S2,
## rV[n]==-1 means that based n (B) is aligned with base n-1 of S2, rV[n]==1 means that base n of S1 is aligned with base n+1 of S2
## l1 is the length of sequence 1, l2 is the length of sequence 2. Iterating through for alignment, one needs a starter row, so the number of rows for iteration is l1+1, which is defined as X.
## the position in sequence 2 is allowed to vary in the alignment by a distance of +/-span of the position in S1,
## So there are a total of 2*span+1 positions allowed for alignment at any position in s1.
## Value j1 gives tha the currently-considered offset between S1 and S2 (positive if the current base in S2 has a greater index than the current base in S1)
## j1 varies between -span and +span inclusive
## i2 is a joint index that consists of the absolute position of the currently-considered alignment pair
## for an ordinal position i1 in S1 and offset j1 for the position in S2, i2 = (i1+1)*Y+(j1+span+1)
## sMax1 and sMax2 are the one-based (ordinal) positions in S1 and S2 where the most favorable alignment terminates.  So sMax1+1 and sMax2+1 are the zero-based indices of the last bases in those alignments
## s1B and s2B are simlarly one-based positions of optimal alignments ending at the terminal base on either s1 or s2 (so s1B should be l1 or s2B should be l2)
## for reverse alignment mV[n] is the value of aligning bases n+1 to the end of S1 (base l1), inclusive.  so mV[0] is the value of aligning base 1 to the end, mV[1] the value of aligning base 2 to the end, etc.  mV[l1] is the value of aligning  nothing, so zero
## for reverse alignment rV[n] is the relative value (offset S1-S2) for base n+1 in an optimal alignment from base n+1 to the end.  So rV[n]==3 means that base n+1 of S1 and base n+1+3 of sequence S2 are paired in the optimal alignment from base n+1 of S1 to the end of each
## for reverse alignment, sMax1 is the first index position participating in optimal alignment.  So sMax1==0 means that every base participates, sMax2=l2 means that no base participates 
## for reverse alignment, sB1 is the first index position participating in the maximum-extend alignment.  So either sB1==0 or sB2==0 should be true in all alignments
## Bond indices-- bond zero is the bond before the initial base in the sequence, bond 1 after the first base, bond len(s) after the last base.  Bond x+1 is the bond after base with index x.

## self.s1JMin is the last base 1-based in optimal alignment, also the first bond that is not in optimal alignment
## thus first bond that might be junction is bond index s1Min
## self.s1JMax is the first base 1-based in optimal alignment on right, says bond index s1JMax is last possible bond as junction
## junction is between s1JMin and s1JMax bond index inclusive
## 

##        If MakeTextDisplay
##            A a single parsed string that describes the alignment
##                # mismatches in s2 and s3 relative to s1
##                M,m mismatches relative to s1 in only s2 or only s3 respectively 
##                v Deletion in s2 and s3 relative to s1
##                D,d deletions relative to s1 in only s2 or only s3 respectively 
##                + Insertion in s2 and s3 relative to s1
##                I,i insertions relative to s1 in only s2 or only s3 respectively 
##                () possible region of junction
##                - both sequenes match
##                * different variants on the two strands
##                . padding base-- terminal sequence not part of alignment (optional... no effect on quality)
##            B1, B2, and B3 are space-delimited versions of s1,s2,s3 designed to be printed on top of A

            
    

print "Welcome to REVAJunctionParser A FireLab product designed to rid your junction feed of fake news"
from time import localtime, strftime, time
from sys import argv
from glob import glob
t0 = time()
now1=strftime('_Date-%m-%d-%y_Time-%H-%M-%S_',localtime())

ReportDeletions1 = True ## Report Putative Deletions?
ReportCircles1 = False ## Report Putative Circles?
MetaName1 = ''  ## Mnemonic to use for to identify output files -- default is a list of input files with extensions trimmed
IgnoreNum1 = 0 ## Number of bases to ignore at the beginning of each read

FL1 = []  ## List of REVA Junction Files to start with

while ai1<len(argv):
    a1=argv[ai1]
    ai1+=1
    a1=a1.replace('"','').replace("'","")
    if not('=') in a1 and not(a1.startswith('-')):
        if len(MetaName1)<100:
            MetaName1 += a1.strip('JunctionEvents_REVA_').replace('*','').split('.')[0]
        else:
            if len(MetaName1) <108:
                MetaName1 += '_'
        if '*' in a1:
            FL1.extend(glob(a1))
        else:
            FL1.append(a1)
    if a1[0]=='-':
        a11=a1.strip()[1:].lower()
        a22=argv[ai1].strip().replace('"','').replace("'","")
        ai1+=1
    else:
        a11=a1.split('=')[0].strip().lower()
        a22=a1.split('=')[-1].strip()
    if a11.startswith('ignore'):
        IgnoreNum1=int(a22)

if nof(FL1):
    print("File not found or no file specified.  Use call syntax syntax:")
    print("python REVAJunctionParser##.py REVAJunctionFileName (IgnoreBases=##))."
    exit()
FGood1 = open(MetaName1.split('.')[0]+now1+'_GoodJunctions.tdv', mode='w')
FBad1 = open(MetaName1.split('.')[0]+now1+'_QuestionableJunctions.tdv', mode='w')
def antisense(s):
    return s.replace('G','c').replace('C','g').replace('A','t').replace('T','a')[::-1].upper()


def unpackB(v0,b=12):
    'output parsed list (score,match,mismatch,indelS,indelE) from input of raw alignment value,  used by nw1 to deliver integral values of score, matches, mismatches, etc'
    ## This is designed to work with sequences up to 1kb with up to 99 distinct insertions/deletions
    score = v0 >> (b*4)
    v1 = v0-(score<<(b*4))
    match = v1>>(b*3)
    v2 = v1-(match<<(b*3))
    mismatch = v2>>(b*2)
    v3 = v2-(mismatch<<(b*2))
    indelE = v3>>b
    indelS = v3-(indelE<<b)
    return [score,match,mismatch,indelS,indelE]
def unpackD(v0,v1,b=12):
    'output [dif-score,dif-match,dif-mismatch,dif-indelS,dif-indelE ] from two alignment raw scores as inputs '
    return [a-b for (a,b) in zip(unpackB(v0),unpackB(v1))]

class A2:
    def __init__(self,s1,s2,span=24,matchV=1,mismatchV=-1,indelSV=-3,indelEV=-1,ReverseAlign=True, LateFavor=True, TwoStep=True, MakeCigA=True):  ## s1 and s2 are sequence strings or byte arrays.  No case checking done so all upper or lower to shart
    ## Two way sequence aligner: Input s1,s2: two sequences. Also optional parameters that we do not recommend changing
    ## penalties/rewards are offset very slightly to take advantage of later 'unpacking' of individual components of the score
    ## TwoStep alignment does a first alignment, gets the best score, then re-aligns sequences downstream of that best score position
    ## MakeCigA makes a string that describes the alignment (self.cigA) with parameters as below
        self.ReverseAlign = ReverseAlign
        self.TwoStep = TwoStep
        self.s1 = s1
        self.s2 = s2
        self.span=span
        if ReverseAlign:
            s1 = s1[::-1]
            s2 = s2[::-1]
        b1 = 12  ## bit length for individual counts of events.  2**b1 is the maximum sequence length for analysis
                 ## for anything like 'current' technology (2018) the implicit maximum of 4096 bases per read seems quite adequate
        maxL = 2**b1
        scoreI = 2**(b1*4)
        matchI = matchV*scoreI+2**(b1*3)
        mismatchI = mismatchV*scoreI+2**(b1*2)
        indelEI = indelEV*scoreI+2**b1
        indelSI = indelSV*scoreI+1
        rd1 = matchI-mismatchI
        pgd1 = indelEI-indelSI
        defaultStart = -2**(b1*5-1)
        self.l1 = min(maxL,len(s1))
        self.l2 = min(maxL,len(s2))
        X = self.l1+1
        Y = 2*span+1
        d1 = [0]*(X*Y) ## d1[(1+n1)*Y+d1] is the highest alignment value for bases zero to n1 of S1 with bases zero to (1+n1+(d1-span-1) of S2.  It is also the value such an alignment contributes to a junction score from the S1 side (thorugh position n1) if the junctino point is set just after position n1 in S1 (zero based)
        e1 = [0]*(X*Y)  ## programmer's note.. you are never going to believe me, but this is really faster with a native Python data structure and direct comparisons (rather than max statements) than with Numpy.  It really is.  You can recode it with Numpy and it will run, but just 20x slower.   
        self.mV = [0] * X
        self.nV = [defaultStart] * X
        self.rV = [0] * X
        self.vMax = defaultStart
        self.s1Max = 0
        self.s2Max = 0
        self.s1B = 0
        self.s2B = 0
        self.vB = defaultStart
        j1 = -span-1
        i1 = 0
        for i2 in xrange(Y,X*Y):  ##i1+1 is the position in the sequence, i1 is the position in array
            j1 += 1
            if j1 > span:
                j1 = -span
                i1 += 1
            if j1+i1<0 or j1+i1>=self.l2:
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
            if (f1>f2 and f1>f3) or (LateFavor and (f1>=f2 and f1>=f3)):
                d1[i2] = f1
                e1[i2] = 1                                
            elif f2>f3 or (LateFavor and f2==f3):
                d1[i2] = f2
                e1[i2] = 2                
            else:
                d1[i2] = f3
                e1[i2] = 3
            if d1[i2]>self.mV[i1+1]: 
                self.mV[i1+1] = d1[i2]
                self.rV[i1+1] = j1
                if d1[i2]>self.vMax or (LateFavor and d1[i2]==self.vMax):
                    self.s1Max = i1+1
                    self.s2Max = i1+j1+1
                    self.vMax = d1[i2]
                    iMax = i2
        if TwoStep:
            for px1 in xrange(iMax-Y+1,iMax):
                d1[px1] = defaultStart
            if (self.s1Max==self.l1 or self.s2Max==self.l2):
                p0 = iMax
                self.vB = d1[iMax]
                self.s1B = self.s1Max
                self.s2B = self.s2Max
            i1 = self.s1Max-1
            j1 = self.s2Max-i1-1
            for i2 in xrange(iMax+1,X*Y):  ##i1+1 is the position in the sequence, i1 is the position in array
                j1 += 1
                if j1 > span:
                    j1 = -span
                    i1 += 1
                if j1+i1<0 or j1+i1>=self.l2:
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
                if (f1>f2 and f1>f3) or (LateFavor and (f1>=f2 and f1>=f3)):
                    d1[i2] = f1
                    e1[i2] = 1                                
                elif f2>f3 or (LateFavor and f2>=f3):
                    d1[i2] = f2
                    e1[i2] = 2                
                else:
                    d1[i2] = f3
                    e1[i2] = 3
                if (d1[i2]>self.vB or (LateFavor and d1[i2]==self.vB)) and (i1==self.l1-1 or i1+j1==self.l2-1):
                    p0 = i2
                    self.vB = d1[i2]
                    self.s1B = i1+1
                    self.s2B = j1+i1+1
        if MakeCigA:
            self.cigA = ['E']*(self.l1-self.s1B)
            self.cigA.extend(['e']*(self.l2-self.s2B))        
            while p0>=Y:
                i1 = p0//Y-1
                j1 = p0-(i1+1)*Y-span
                if d1[p0] > self.nV[i1+1]:
                    self.nV[i1+1] = d1[p0]
                if e1[p0]==1:
                    if i1+j1<self.l2 and s1[i1]==s2[i1+j1]:
                        self.cigA.append('M')
                    else:
                        self.cigA.append('m')
                    p0-=Y
                    self.rV[i1+1] = j1
                elif e1[p0] == 2:
                    self.cigA.append('d')
                    p0-=1
                elif e1[p0] == 3:
                    self.cigA.append('D')
                    p0 -= Y-1
                else:
                    self.cigA.extend('E'*(p0//Y))
                    p0=span
                    break
            self.cigA.extend(['e']*(p0-span))
            self.cigA = ''.join(self.cigA[::-1])
        if ReverseAlign:
            self.mV = self.mV[::-1]
            self.rV = [self.l1-self.l2-x for x in self.rV[::-1]]
            self.nV = self.nV[::-1]
            self.s1Max = self.l1-self.s1Max
            self.s2Max = self.l2-self.s2Max
            if TwoStep:
                self.s1B = self.l1-self.s1B
                self.s2B = self.l2-self.s2B
            if MakeCigA:
                self.cigA = self.cigA[::-1]
            


class A3:
    def __init__(self,s1,s2,s3,JLen = 25,MakeTextDisplay=True,SenseStrand=True):
        self.aF12 = A2(s1,s2,ReverseAlign=False,LateFavor=not(SenseStrand))  ## forward alignment, sequence s2
        self.aR13 = A2(s1,s3,ReverseAlign=True,LateFavor=SenseStrand)
        self.JLen = JLen
        self.pivot = self.aR13.s1Max-self.aF12.s1Max ## number of bases to be inserted between best alighment of s2 and best alignment of s3.  zero for a flush alignment
        if self.pivot>=0:
            self.s1JMin = self.aF12.s1Max ## first bond index for possible junction.  index of first S1 base that is not in the left optimal alignement.  The bond before t
            self.s1JMax = self.aR13.s1Max ## index of bond that is the last bond that is not connecting bases in the optimum right alignment.  S1 base that is not in optimal right alignment
            self.s2JMin = self.aF12.s2Max
            self.s2JMax = self.aF12.s2Max
            self.s3JMin = self.aR13.s2Max
            self.s3JMax = self.aR13.s2Max
            self.JIns = s1[self.s1JMin:self.s1JMax]
            self.JDup = ''
            self.s2J5p = s2[max(0,self.s2JMax-JLen):self.s2JMax]
            self.s2J3p = s2[self.s2JMax:self.s2JMax+JLen]
            self.s3J5p = s3[max(0,self.s3JMax-JLen):self.s3JMax]
            self.s3J3p = s3[self.s3JMax:self.s3JMax+JLen]
            self.qual12 = unpackB(self.aF12.vB)
            self.qual13 = unpackB(self.aR13.vB)
            self.qMatchJF12 = unpackB(self.aF12.vMax)
            self.qMatchJR13 = unpackB(self.aR13.vMax)
            self.qMisMatchJF12 = unpackD(self.aF12.vB,self.aF12.vMax)
            self.qMisMatchJR13 = unpackD(self.aR13.vB,self.aR13.vMax)
        else:
            self.JScore = map(sum,zip(self.aF12.nV,self.aR13.nV))  ## differential value of including position j in left arm
            self.JScoreMax = max(self.JScore)
            self.s1JMin = self.JScore.index(self.JScoreMax)
            self.s2JMin = self.s1JMin+self.aF12.rV[self.s1JMin]
            self.s3JMin = self.s1JMin+self.aR13.rV[self.s1JMin]
            self.s1JMax = len(s1)-self.JScore[::-1].index(self.JScoreMax)
            self.s2JMax = self.s1JMax+self.aF12.rV[self.s1JMax]
            self.s3JMax = self.s1JMax+self.aR13.rV[self.s1JMax]
            self.JIns = ''
            self.JDup = s1[self.s1JMin:self.s1JMax]
            self.s2J5p = s2[max(0,self.s2JMin-JLen):self.s2JMin]
            self.s2J3p = s2[self.s2JMax:self.s2JMax+JLen]
            self.s3J5p = s3[max(0,self.s3JMin-JLen):self.s3JMin]
            self.s3J3p = s3[self.s3JMax:self.s3JMax+JLen]
            self.qual12 = unpackB(self.aF12.vB)
            self.qual13 = unpackB(self.aR13.vB)
            self.qMatchJF12 = unpackB(self.aF12.mV[self.s1JMax])
            self.qMatchJR13 = unpackB(self.aR13.mV[self.s1JMin])
            self.qMisMatchJF12 = unpackD(self.aF12.vB,self.aF12.mV[self.s1JMax])
            self.qMisMatchJR13 = unpackD(self.aR13.vB,self.aR13.mV[self.s1JMin])
        if MakeTextDisplay:
            i1=0
            i2=0
            i3=0
            i12=0
            i13=0
            A = []
            B1 = []
            B2 = []
            B3 = []
            D1 = {('M','M'):('-',1,1,1,1,1),
                  ('M','m'):('m',1,1,1,1,1),
                  ('m','M'):('M',1,1,1,1,1),
                  ('m','m'):('#',1,1,1,1,1),
                  ('d','d'):('v',1,1,0,1,1),
                  ('D','D'):('+',1,1,1,0,0),
                  ('d','M'):('I',1,0,0,1,0),
                  ('d','m'):('I',1,0,0,1,0),
                  ('d','D'):('I',1,0,0,1,0),
                  ('M','d'):('i',0,1,0,0,1),
                  ('m','d'):('i',0,1,0,0,1),
                  ('D','d'):('i',0,1,0,0,1),
                  ('D','M'):('D',1,1,1,0,1),
                  ('D','m'):('*',1,1,1,0,1),
                  ('M','D'):('d',1,1,1,1,0),
                  ('m','D'):('*',1,1,1,1,0)}
            cig12 = self.aF12.cigA.replace('e','d').replace('E','D')
            cig13 = self.aR13.cigA.replace('e','d').replace('E','D')
            fP1 = False ## Have we put in forward parenthesis yet to indicate first possible junction point 
            rP1 = False ## Have we put in reverse parenthesis yet to indicate last possible junction point 
            while True:
                if i12>=len(cig12) and i13>=len(cig13):
                    break
                if i12<len(cig12):
                    c12 = cig12[i12]
                else:
                    c12 = 'D'
                if i13<len(cig13):
                    c13 = cig13[i13]
                else:
                    c13 = 'D'
                a1,inc12,inc13,inc1,inc2,inc3=D1[(c12,c13)]
                if self.pivot<=0 and i1>=self.s1JMin and i2>=self.s2JMin and i3>=self.s3JMin and not(fP1):
                    A+='['                        
                    B1+='|'
                    B2+='|'
                    B3+='|'
                    fP1 = True
                if self.pivot>0 and i1>=self.s1JMin and i2>=self.s2JMin and not(fP1):
                    A+='('                        
                    B1+='|'
                    B2+='|'
                    B3+='|'
                    fP1 = True
                if inc1==1 and i1<len(s1):
                    if c12=='m' and c13=='m':
                        B1+=s1[i1].upper()
                    else:
                        B1+=s1[i1].lower()            
                else:
                    B1+=' '
                if inc2==1 and i2<len(s2):
                    if c12=='m' and c13=='M':
                        B2+=s2[i2].upper()
                    else:
                        B2+=s2[i2].lower()
                else:
                    B2+=' '
                if inc3==1 and i3<len(s3):
                    if c12=='M' and c13=='m':
                        B3+=s3[i3].upper()
                    else:
                        B3+=s3[i3].lower()
                else:
                    B3+=' '
                A += a1
                i12+=inc12
                i13+=inc13
                i1+=inc1
                i2+=inc2
                i3+=inc3
                if self.pivot<=0 and (i1>self.s1JMax or i2>self.s2JMax or i3>self.s3JMax) and not(rP1):
                    A.insert(-1,']')                        
                    B1.insert(-1,'|')
                    B2.insert(-1,'|')
                    B3.insert(-1,'|')
                    rP1 = True
                if self.pivot>0 and (i1>self.s1JMax or i3>self.s3JMax) and not(rP1):
                    A.insert(-1,')')                        
                    B1.insert(-1,'|')
                    B2.insert(-1,'|')
                    B3.insert(-1,'|')
                    rP1 = True

            for x,c in enumerate(A):
                if c in 'dDiI+v':
                    A[x]='.'
                else:
                    break
            for x,c in enumerate(A[::-1]):
                if c in 'dDiI+v':
                    A[len(A)-1-x]='.'
                else:
                    break
            self.A = ''.join(A)
            self.B1 = ''.join(B1)
            self.B2 = ''.join(B2)
            self.B3 = ''.join(B3)
##            print s1[:self.s1JMin]
##            print ' '*self.s1JMin+s1[self.s1JMin:self.s1JMax]
##            print ' '*self.s1JMax+s1[self.s1JMax:]
##            print self.B1
##            print self.B2
##            print self.B3
##            print self.A
##        public parameters
##            pivot: ## number of bases to be inserted between best alighment of s2 and best alignment of s3.  zero for a flush alignment
##            s1JMin,s2JMin,s3JMin: Earliest optimal junction positions (in s1,s2,s3)
##            s1JMax,s2JMax,s3JMax: Latest  optimal junction positions (in s1,s2,s3)
##            s1JMax: Latest position in S1 where junction fits 
##            JIns: Reconstruction of sequences inserted at junction from neither s2 nor s3 
##            JDup: Sequences shared between s2 and s3 that are ducplicated a the junction 
##            s2J5p: Sequences in S2 5' to earliest optimal junction 
##            s2J3p: Sequences in S2 3' to latest optimal junction
##            s3J5p: Sequences in S3 5' to earliest optimal junction  
##            s3J3p: Sequences in S3 5' to earliest optimal junction  
##            qual12: Quality scores (Score, Matches, Mismatches, Indel Count, Indel Total Length) for whole-length match of s1 and s2
##            qual13: Quality scores (Score, Matches, Mismatches, Indel Count, Indel Total Length) for whole-length match of s1 and s3 
##            qMatchJF12 Quality scores (Score, Matches, Mismatches, Indel Count, Indel Total Length) for best foreward match of s1 and s2
##            qMatchJR13 Quality scores (Score, Matches, Mismatches, Indel Count, Indel Total Length) for best reverse match of s1 and s3
##            qMisMatchJF12 Differential quality score-- how much does whole sequence lose from match quality of best match (ie how convincing is it that we actually need a junction as contrasted to just errors in sequence)
##            qMisMatchJR13
##        If MakeTextDisplay
##            A a single parsed string that describes the alignment
##                # mismatches in s2 and s3 relative to s1
##                M,m mismatches relative to s1 in only s2 or only s3 respectively 
##                v Deletion in s2 and s3 relative to s1
##                D,d deletions relative to s1 in only s2 or only s3 respectively 
##                + Insertion in s2 and s3 relative to s1
##                I,i insertions relative to s1 in only s2 or only s3 respectively 
##                () possible region of junction
##                - both sequenes match
##                * different variants on the two strands
##                . padding base-- terminal sequence not part of alignment (optional... no effect on quality)
##            B1, B2, and B3 are space-delimited versions of s1,s2,s3 designed to be printed on top of A

orderingstring = '#*-+IiDdMm()[]vGATC'  ## reverse order of preference in reconciling alignments if a "tie" is detected
## note that bases beyond the alignment aren't counted at all '.'
def mostcommonchar(S):
    ## S is a list or string. Return the most frequently used character
    ## For reproducibility the list is sorted first and the last element is returne in case of ties
    D = {}
    for c in S:
        if not(c in D):
            D[c]=0
        D[c]+=1
    mV = 0
    mC = '.'
    for c in orderingstring:
        if c in D and D[c]>mV:
            mC = c
            mV = D[c]
    sstr = ''.join(S).replace('.','').replace(',','')
    ValList = [D[x] for x in D if x!='.' and x!=',']
    ValSum = sum(ValList)
    CoincidenceOpportunities = ValSum**2-ValSum
    Coincidences = sum([x**2-x for x in ValList])    
    return mC,CoincidenceOpportunities,Coincidences
def alignreverse(aL):
    return aL[::-1].swapcase().replace('[','1').replace('(','2').replace(']','[').replace(')','(').replace('1',']').replace('2',')')
def alignreconcile(aL):
    ## Calculate a best case alignment and also provide a rough value to assess concordance between different alignment assessments.
    CoincidenceOpportunitiesSum = 0
    CoincidencesSum = 0
    if len(aL)==0:
        return '',0,0
    L = len(aL)
    nL = []
    dL = ''
    for a1 in aL:
        nL.append(max(a1.find('['),a1.find('(')))
    pL = nL[:]
    while True:
        nC = ['.']  ## list of next character
        for i in xrange(L):
            if len(aL[i])>pL[i]:
                nC.append(aL[i][pL[i]])
        cN,coO,co = mostcommonchar(nC)
        if not(cN in '[]().,'):
            CoincidenceOpportunitiesSum += coO
            CoincidencesSum += co
        dL += cN
        if cN == '.':
            break
        if cN in 'dDv':
            for i in range(L):
                if len(aL[i])>pL[i] and aL[i][pL[i]] in 'dDv':
                    pL[i]+=1
        else:
            for i in range(L):
                while len(aL[i])>pL[i] and aL[i][pL[i]] in 'dDv':
                    pL[i]+=1
                pL[i]+=1
    pL = nL[:]
    uL=''
    while True:
        nC = ['.']  ## list of next character
        for i in xrange(L):
            if pL[i]>=0:
                nC.append(aL[i][pL[i]])
        cN,coO,co = mostcommonchar(nC)
        if not(cN in '[]().,'):
            CoincidenceOpportunitiesSum += coO
            CoincidencesSum += co
        uL += cN
        if cN == '.':
            break
        if cN in 'dDv':
            for i in range(L):
                if pL[i]>=0 and aL[i][pL[i]] in 'dDv':
                    pL[i]-=1
        else:
            for i in range(L):
                while pL[i]>=0 and aL[i][pL[i]] in 'dDv':
                    pL[i]-=1
                pL[i]-=1
    CoincidenceScore = CoincidencesSum/(0.0000001+CoincidenceOpportunitiesSum)
    return (uL[1:][::-1]+dL).strip('.'),CoincidenceScore

def consensus(a1):  ## assumes all sequences are the same length
    ## Calculate a best case alignment and also provide a rough value to assess concordance between different alignment assessments.
    CoincidenceOpportunitiesSum = 0
    CoincidencesSum = 0
    if len(a1)==0:
        return '',0.0
    L = len(a1)
    dL = ''
    k = 0
    while True:
        nC = []  ## list of next character
        for i in xrange(L):
            if len(a1[i])>k:
                nC.append(a1[i][k])
        if not(nC):
            break
        cN,coO,co = mostcommonchar(nC)
        CoincidenceOpportunitiesSum += coO
        CoincidencesSum += co
        if len(nC)>L//2:
            dL += cN
        k+=1
    CoincidenceScore = CoincidencesSum/(0.0000001+CoincidenceOpportunitiesSum)
    return dL,CoincidenceScore        
        

def alignquality(a1):
    ## A very simple view of match quality between putative aggregate junstion and actual read
    a1 = a1.strip('.')[1:-1].strip('dDiI+v*#mM')
    leadingN1 = max(a1.find('['),a1.find('('))
    laggingN1 = max(a1.find(']'),a1.find(')'))
    Lead1 = a1[:leadingN1]
    Lag1 = a1[laggingN1:]
    if Lag1 and Lag1[0]==']' and laggingN1>leadingN1-1:
        Overlap1 = a1[laggingN1:leadingN1]
    else:
        Overlap1 =''
    
    if Lag1 and Lag1[0]==')' and laggingN1>leadingN1-1:
        Insertion1 = a1[laggingN1:leadingN1]
    else:
        Insertion1 =''
    Mismatch1 = (Lead1.count('M')+
                 Lead1.count('I')+
                 Lead1.count('D')+
                 Lag1.count('m')+
                 Lag1.count('i')+
                 Lag1.count('d')+
                 a1.count('+')+
                 a1.count('v')+
                 a1.count('*')+
                 a1.count('#')+
                 Overlap1.count('M')+
                 Overlap1.count('m')+
                 Overlap1.count('D')+
                 Overlap1.count('d')+
                 Overlap1.count('I')+
                 Overlap1.count('i'))
    Match1 = a1.count('-')
    LeftInconsistent1 = (Lead1.count('M')+
                 Lead1.count('I')+
                 Lead1.count('D')+
                 Overlap1.count('M')+
                 Overlap1.count('I')+
                 Overlap1.count('D'))
    RightInconsistent1 = (Lag1.count('m')+
                 Lag1.count('i')+
                 Lag1.count('d')+
                 Overlap1.count('m')+
                 Overlap1.count('i')+
                 Overlap1.count('d'))
    LeftSupportive1 = (Lead1.count('m')+
                 Lead1.count('i')+
                 Lead1.count('d'))    
    RightSupportive1 = (Lag1.count('M')+
                 Lag1.count('I')+
                 Lag1.count('D'))    
    return Match1, Mismatch1,LeftInconsistent1,RightInconsistent1,LeftSupportive1,RightSupportive1
                 
    
x0=0
y0=0
y1=0
y2=0
D1 = {}
GoodyHeader1 = '\t'.join(['JLeftLow',
                    'JLeftHigh',
                    'JRightLow',
                    'JRightHigh',
                    'Insertion',
                    'Duplication',
                    'PivotDistance',
                    'Left_5p',
                    'Left_3p',
                    'Right_5p',
                    'Right_3p',
                    'MatchString'])

MetaD1 = {}  ## Dictionary are unique rearrangments indexed by <chr, start, end>
             ## Values are lists of rearrangement 
MetaD2 = {}  ## Dictionary indices are rearrangments
             ## Values are lists of alignments 
MetaD3 = {}  ## Dictionary indices are rearrangments
             ## Values are lists of flanking sequences (Ups 5', Ups 3', Dwn 5', Dwn 3', ) 
MetaD4 = {}  ## Dictionary indices are rearrangments
             ## Values are lists of duplications
MetaD5 = {}  ## Dictionary indices are rearrangments
             ## Values are lists of pivot distances
MetaD6 = {}  ## Dictionary indices are rearrangments
             ## Values are ucsc links

HeadersWritten1 = False
for Fn1 in FL1:
    F1 = open(Fn1,mode='rU')
    print("Opening File "+Fn1+'. Time='+"{0:.2f}".format(time()-t0)+' sec')
    SourceFile1 = F1.name
    if len(SourceFile1.split('_'))>2:
        SourceFile1 = SourceFile1.split('_')[2]
    dc00 = 0
    for g1,L0 in enumerate(F1):
        L1 = L0.strip().split('\t')
        if not(HeadersWritten1) and len(L1)>20:
            FGood1.write(L0.strip()+'\t')
            FGood1.write(GoodyHeader1+'\r')
            FBad1.write(L0.strip()+'\t')
            FBad1.write(GoodyHeader1+'\r')
            HeadersWritten1 = True            
        SkipThisLine1 = True
        if len(L1)>30 and ((ReportDeletions1 and ('Deletion_Junction' in L1[0])) or (ReportCircles1 and ('Circle_Junction' in L1[0]))):
            SkipThisLine1 = False
        if SkipThisLine1:
            continue
        OutFile=FBad1
        dc00 += 1
        if dc00%10000==0:
            print("Completing Candidate "+str(dc00)+' from '+SourceFile1+' . Time='+"{0:.2f}".format(time()-t0)+' sec')
        if len(L1)<31 or not(('#' in L1[30]) and (':' in L1[30])):
            SampleID1 = SourceFile1
        else:
            SampleID1 = L1[30].split('#')[0].split(':')[-1]        
        if ('1s1s' in L1[0]) or ('1a1a' in L1[0]):
            s1 = L1[8][IgnoreNum1:]
            s2 = L1[12][IgnoreNum1:]
            s3 = L1[16][IgnoreNum1:]
            a123 = A3(s1,s2,s3,SenseStrand=('1s1s' in L1[0]))
            J12Inconsistent = a123.qMatchJF12[2]+a123.qMatchJF12[4]
            J12Supportive = a123.qMisMatchJF12[2]+a123.qMatchJF12[4]
            J13Inconsistent = a123.qMatchJR13[2]+a123.qMatchJR13[4]
            J13Supportive = a123.qMisMatchJR13[2]+a123.qMatchJR13[4]
            if not(J12Supportive<2 or
               J13Supportive<2 or
               J12Inconsistent>2 or
               J13Inconsistent>2 or
               J12Supportive<=J13Inconsistent+2 or
               J13Supportive<=J12Inconsistent+2 or
               a123.qMatchJF12[2]>3 or
               a123.qMatchJR13[2]>3):
                OutFile=FGood1
            if '1s1s' in L1[0]:
                PreciseEventStartLo1 = int(L1[10])+a123.s2JMin
                PreciseEventStartHi1 = int(L1[10])+a123.s2JMax
                PreciseEventEndLo1 = int(L1[14])+a123.s3JMin
                PreciseEventEndHi1 = int(L1[14])+a123.s3JMax
                DupSeq = a123.JDup
                InsSeq = a123.JIns
                JU5p = a123.s2J5p
                JU3p = a123.s2J3p
                JD5p = a123.s3J5p
                JD3p = a123.s3J3p
                if OutFile==FGood1:
                    Index1 = (L1[9],PreciseEventStartLo1,PreciseEventStartHi1,L1[13],PreciseEventEndLo1,PreciseEventEndHi1,InsSeq)
                    if not Index1 in MetaD1:
                        MetaD1[Index1] = {}
                        MetaD2[Index1] = []
                        MetaD3[Index1] = ['','','','']
                        MetaD4[Index1] = []
                        MetaD5[Index1] = []
                        MetaD6[Index1] = L1[7]
                    MetaD4[Index1].append(DupSeq)
                    MetaD5[Index1].append(a123.pivot)
                    Index2 = '_'.join(('1s',str(a123.s2JMin),SampleID1))
                    if not(Index2 in MetaD1[Index1]):
                        MetaD1[Index1][Index2] = 0
                    MetaD1[Index1][Index2] += 1
                    MetaD2[Index1].append(a123.A)
                    if len(JU5p)>len(MetaD3[Index1][0]):
                        MetaD3[Index1][0]=JU5p
                    if len(JU3p)>len(MetaD3[Index1][1]):
                        MetaD3[Index1][1]=JU3p
                    if len(JD5p)>len(MetaD3[Index1][2]):
                        MetaD3[Index1][2]=JD5p
                    if len(JD3p)>len(MetaD3[Index1][3]):
                        MetaD3[Index1][3]=JD3p
            elif '1a1a' in L1[0]:
                PreciseEventEndLo1 = -int(L1[10])-a123.s2JMax+1
                PreciseEventEndHi1 = -int(L1[10])-a123.s2JMin+1
                PreciseEventStartLo1 = -int(L1[14])-a123.s3JMax+1
                PreciseEventStartHi1 = -int(L1[14])-a123.s3JMin+1
                DupSeq = antisense(a123.JDup)
                InsSeq = antisense(a123.JIns)
                JU5p = antisense(a123.s3J3p)
                JU3p = antisense(a123.s3J5p)
                JD5p = antisense(a123.s2J3p)
                JD3p = antisense(a123.s2J5p)
                if OutFile==FGood1:
                    Index1 = (L1[13],PreciseEventStartLo1,PreciseEventStartHi1,L1[9],PreciseEventEndLo1,PreciseEventEndHi1,InsSeq)##,JU5p,JU3p,JD5p,JD3p)
                    if not Index1 in MetaD1:
                        MetaD1[Index1] = {}
                        MetaD2[Index1] = []
                        MetaD3[Index1] = ['','','','']
                        MetaD4[Index1] = []
                        MetaD5[Index1] = []
                        MetaD6[Index1] = L1[7]
                    MetaD4[Index1].append(DupSeq)
                    MetaD5[Index1].append(a123.pivot)
                    Index2 = '_'.join(('1a',str(a123.s2JMin),SampleID1))
                    if not(Index2 in MetaD1[Index1]):
                        MetaD1[Index1][Index2] = 0
                    MetaD1[Index1][Index2] += 1
                    MetaD2[Index1].append(alignreverse(a123.A))
                    if len(JU5p)>len(MetaD3[Index1][0]):
                        MetaD3[Index1][0]=JU5p
                    if len(JU3p)>len(MetaD3[Index1][1]):
                        MetaD3[Index1][1]=JU3p
                    if len(JD5p)>len(MetaD3[Index1][2]):
                        MetaD3[Index1][2]=JD5p
                    if len(JD3p)>len(MetaD3[Index1][3]):
                        MetaD3[Index1][3]=JD3p
            Goodies = [ PreciseEventStartLo1,
                        PreciseEventStartHi1,
                        PreciseEventEndLo1,
                        PreciseEventEndHi1,
                        InsSeq,
                        DupSeq,
                        -a123.pivot,
                        JU5p,
                        JU3p,
                        JD5p,
                        JD3p,
                        a123.A]
            L0 = L0.strip() + '\t' + '\t'.join(map(str,Goodies))+'\r'
        if ('2s2s' in L1[0]) or ('2a2a' in L1[0]):
            s1 = L1[19]
            s2 = L1[23]
            s3 = L1[27]
            a123 = A3(s1,s2,s3,SenseStrand=('2s2s' in L1[0]))
            J12Inconsistent = a123.qMatchJF12[2]+a123.qMatchJF12[4]
            J12Supportive = a123.qMisMatchJF12[2]+a123.qMatchJF12[4]
            J13Inconsistent = a123.qMatchJR13[2]+a123.qMatchJR13[4]
            J13Supportive = a123.qMisMatchJR13[2]+a123.qMatchJR13[4]
            if not(J12Supportive<2 or
               J13Supportive<2 or
               J12Inconsistent>2 or
               J13Inconsistent>2 or
               J12Supportive<=J13Inconsistent+2 or
               J13Supportive<=J12Inconsistent+2 or
               a123.qMatchJF12[2]>3 or
               a123.qMatchJR13[2]>3):
                OutFile=FGood1
            if '2s2s' in L1[0]:
                PreciseEventStartLo1 = int(L1[21])+a123.s2JMin
                PreciseEventStartHi1 = int(L1[21])+a123.s2JMax
                PreciseEventEndLo1 = int(L1[25])+a123.s3JMin
                PreciseEventEndHi1 = int(L1[25])+a123.s3JMax
                DupSeq = a123.JDup
                InsSeq = a123.JIns
                JU5p = a123.s2J5p
                JU3p = a123.s2J3p
                JD5p = a123.s3J5p
                JD3p = a123.s3J3p
                if OutFile==FGood1:
                    Index1 = (L1[20],PreciseEventStartLo1,PreciseEventStartHi1,L1[24],PreciseEventEndLo1,PreciseEventEndHi1,InsSeq)
                    if not Index1 in MetaD1:
                        MetaD1[Index1] = {}
                        MetaD2[Index1] = []
                        MetaD3[Index1] = ['','','','']
                        MetaD4[Index1] = []
                        MetaD5[Index1] = []
                        MetaD6[Index1] = L1[7]
                    MetaD4[Index1].append(DupSeq)
                    MetaD5[Index1].append(a123.pivot)
                    Index2 = '_'.join(('2s',str(a123.s2JMin),SampleID1))
                    if not(Index2 in MetaD1[Index1]):
                        MetaD1[Index1][Index2] = 0
                    MetaD1[Index1][Index2] += 1
                    MetaD2[Index1].append(a123.A)
                    if len(JU5p)>len(MetaD3[Index1][0]):
                        MetaD3[Index1][0]=JU5p
                    if len(JU3p)>len(MetaD3[Index1][1]):
                        MetaD3[Index1][1]=JU3p
                    if len(JD5p)>len(MetaD3[Index1][2]):
                        MetaD3[Index1][2]=JD5p
                    if len(JD3p)>len(MetaD3[Index1][3]):
                        MetaD3[Index1][3]=JD3p
            elif '2a2a' in L1[0]:
                PreciseEventEndLo1 = -int(L1[21])-a123.s2JMax+1
                PreciseEventEndHi1 = -int(L1[21])-a123.s2JMin+1
                PreciseEventStartLo1 = -int(L1[25])-a123.s3JMax+1
                PreciseEventStartHi1 = -int(L1[25])-a123.s3JMin+1
                DupSeq = antisense(a123.JDup)
                InsSeq = antisense(a123.JIns)
                JU5p = antisense(a123.s3J3p)
                JU3p = antisense(a123.s3J5p)
                JD5p = antisense(a123.s2J3p)
                JD3p = antisense(a123.s2J5p)
                if OutFile==FGood1:
                    Index1 = (L1[24],PreciseEventStartLo1,PreciseEventStartHi1,L1[20],PreciseEventEndLo1,PreciseEventEndHi1,InsSeq)
                    if not Index1 in MetaD1:
                        MetaD1[Index1] = {}
                        MetaD2[Index1] = []
                        MetaD3[Index1] = ['','','','']
                        MetaD4[Index1] = []
                        MetaD5[Index1] = []
                        MetaD6[Index1] = L1[7]
                    MetaD4[Index1].append(DupSeq)
                    MetaD5[Index1].append(a123.pivot)
                    Index2 = '_'.join(('2a',str(a123.s2JMin),SampleID1))
                    if not(Index2 in MetaD1[Index1]):
                        MetaD1[Index1][Index2] = 0
                    MetaD1[Index1][Index2] += 1
                    MetaD2[Index1].append(alignreverse(a123.A))
                    if len(JU5p)>len(MetaD3[Index1][0]):
                        MetaD3[Index1][0]=JU5p
                    if len(JU3p)>len(MetaD3[Index1][1]):
                        MetaD3[Index1][1]=JU3p
                    if len(JD5p)>len(MetaD3[Index1][2]):
                        MetaD3[Index1][2]=JD5p
                    if len(JD3p)>len(MetaD3[Index1][3]):
                        MetaD3[Index1][3]=JD3p

            Goodies = [ PreciseEventStartLo1,
                        PreciseEventStartHi1,
                        PreciseEventEndLo1,
                        PreciseEventEndHi1,
                        InsSeq,
                        DupSeq,
                        -a123.pivot,
                        JU5p,
                        JU3p,
                        JD5p,
                        JD3p,
                        a123.A]
            L0 = L0.strip() + '\t' + '\t'.join(map(str,Goodies))+'\r'
        OutFile.write(L0)
FGood1.close()
FBad1.close()
MetaOut1 = open(MetaName1.split('.')[0]+now1+'MetaSummary.tdt',mode='w')
MetaOutHeader1 = 'Ups_Chr\tUps_Lo\tUps_Hi\tDwn_Chr\tDwn_Lo\tDwn_Hi\tInsSeq\tDupSeq\tPivot_Distance\tDup_Qual\tUps_5p\tUps_3p\tDwn_5p\tDwn_3p\tAlign_Types\tAlign_Strings\tAlign_Correlation\tAlign_Aggregate\tMatch\tMismatch\tLeftInconsistent\tRightInconsistent\tLeftSupportive\tRightSupportive\tInstances\tUniqueSettings\tRearrangment\tUCSC_Link(FirstInstance)\r'            
MetaOut1.write(MetaOutHeader1)
for Index1 in sorted(MetaD1.keys()):
    del1 = Index1[4]-Index1[1]
    if del1>0:
        Mnemonic1 = 'Del_'+str(del1)
    else:
        Mnemonic1 = 'Cir_'+str(-del1)
    aR1,aQ0 = alignreconcile(MetaD2[Index1])
    aQ1 = alignquality(aR1)
    aov1,aoq1 = consensus(MetaD4[Index1])
    if MetaD5[Index1]==[]:
        Pivot1 = 0
    else:
        Pivot1 = -sum(MetaD5[Index1])/(1.0*len(MetaD5[Index1]))
    HitsAll1 = sum(MetaD1[Index1].values())
    HitsUnique1 = len(MetaD1[Index1])
    MetaOut1.write('\t'.join(map(str,Index1))+'\t')
    MetaOut1.write(aov1+'\t'+ ("%.2f" % Pivot1) +'\t'+("%.3f" % aoq1) +'\t')
    MetaOut1.write('\t'.join(map(str,MetaD3[Index1]))+'\t')
    MetaOut1.write(str(MetaD1[Index1]).replace("'","")[1:-1]+'\t')
    if MetaD1[Index1]=={}: aardvark
    MetaOut1.write(' '.join(MetaD2[Index1])+'\t')
    MetaOut1.write("%.3f" % aQ0+'\t')
    MetaOut1.write(aR1+'\t')
    MetaOut1.write('\t'.join(map(str,aQ1))+'\t')
    MetaOut1.write(str(HitsAll1)+'\t'+str(HitsUnique1)+'\t'+Mnemonic1+'\t')
    MetaOut1.write(MetaD6[Index1]+'\r')
MetaOut1.close()

##tetritisK1 = 12
##tetritisEndFuzzy = 4
##tetritis=open('illuminatetritis.fa',mode='rU').read()
##
##tetritis=tetritis.splitlines()
##T1={}  ## keys are tetritisK mers that are in at least one illumina linker/adaptor or its complement, values (at the moment) are just zero
##for l1 in tetritis:
##    l2=l1.strip().upper()
##    a2=antisense(l2)
##    for i in range(len(l2)-tetritisK1+1):
##        T1[l2[i:i+tetritisK1]]=1
##        T1[a2[i:i+tetritisK1]]=1
##
 
##    OffTetritis1 = endL1-tetritisEndFuzzy
##    sU1 =s1[endL1-tetritisEndFuzzy:startR1+tetritisEndFuzzy]
##    tStart1 = 0
##    tEnd1 = 0
##    for t in range(len(sU1)-tetritisK1+1):
##        if sU1[t:t+tetritisK1] in T1:
##            if tStart1 == 0:
##                tStart1 = t
##            tEnd1 = t+tetritisK1
##    if tStart1>0:
##        tStart1 += OffTetritis1
##        tEnd1 += OffTetritis1
