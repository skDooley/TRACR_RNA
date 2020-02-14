import tempfile
from sys import path as spath
spath.append("scripts/")

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqIO import index as fasta_index, parse, write
from Bio.SeqRecord import SeqRecord
from CRISPRtools import * #MakeFasta, PilerCRReader, MinCEDReader
from glob import glob
from easyFunctions import BLAST_short, comma, Coordinate, dump
from InfernalResults import TRACR_RNA, ProcessInfernal
from HMMParser import *
from os import chdir, path, stat, system
from pandas import Series
from pickle import load
from random import randint
from Rho import *
from RhoTermPredict import RhoTermPredict

################################         Data          ################################
print("Loaded libraries")
chdir("/mnt/research/germs/shane/transActRNA/data")
casRelatedAssemblies = {}
gene = "Cas9"
casRelatedAssemblies = fasta_index("proteins/%s-Like-clustered.faa" % (gene),"fasta")
print("Number of %ss: %i" % (gene,len(casRelatedAssemblies)))

################################    Custom Methods     ################################
r = lambda: randint(0,255)            
def genColor():return '#%02X%02X%02X' % (r(),r(),r())
def findBiggestCluster(clusters, clusterMap): return Counter(dict((k, clusterMap[k]) for k in (clusters))).most_common(1)[0][0] 
def duplexPercIdent(tracrRNA,crRNA):
    crLen = len(crRNA)
    if len(tracrRNA) >= crLen: return
    fmatches,rmatches = 0,0
    for i in range(crLen): 
        fmatches += int(tracrRNA[i]==crRNA[i])
        rmatches += int(tracrRNA[-(i+1)]==crRNA[-(i+1)])
    revComp = RC(crRNA)
    rfmatches,rrmatches = 0,0
    for i in range(crLen): 
        rfmatches += int(tracrRNA[i]==revComp[i])
        rrmatches += int(tracrRNA[-(i+1)]==revComp[-(i+1)])
    crLen = float(crLen)
    return max(fmatches/crLen, rmatches/crLen, rfmatches/crLen, rrmatches/crLen)*100

print("Loading Cas Operons")
casOperons = load(open("pickles/%s_Operons.p" % gene,"rb"))
print("Loading Cas assemblies")
allAssemblies = dict(fasta_index("assemblies/All_%s-Like-filtered.faa" % (gene),'fasta'))
breakPoints = set(range(100,len(allAssemblies),100))
print("Total # of Cas9s: %s" % (comma(len(allAssemblies))))
noPredictedTracr = {} # {ID:Reason it didn't have a tracr}
totalSols, erpSols, breakCount, hadToGetSeq = 0, 0, 0, 0
die = False
possibleSol = open("sequences/tempAll_%s_Predicted_TracrRNAs.fasta" % (gene),"w")
with open("FindTracrsLog.txt","w") as progress: progress.write("%i\n" % (len(allAssemblies)))
try:
    for i, protID in enumerate(allAssemblies): #     protID = 'NFJO01000020_ORF977'
        if protID == "CP010309_modified_ORF105548": continue
        if i in breakPoints: 
            with open("FindTracrsLog.txt","a") as progress: progress.write("%i " % (i))
            print(i,end=' ')

        # Step 1. Get the CRISPR operon = casOperons[protID]
        operon = casOperons.operons[casOperons.seqMap[protID]]
        crispr = operon.getCRISPR(protID)

        # Step 2. Write all consensus repeats to a file
        crispr.repeatSeqs(protID,open(operon.getRepeatPath(protID),'w'))

        # Step 3. Blast the consensus repeats against the Cas-like protein-containing chromosome
        blastResults = parseBLAST(BLAST_short(operon.getRepeatPath(protID), operon.getFastaPointer(protID), "blastout/conRepeats/%s.xml" % (protID)))
        if len(blastResults) == 0: 
            noPredictedTracr[protID] = 'No BLAST results'
            continue
        
        #Step 4. Narrow the blast results down by removing crRNAs
        crispr.clusterBLASTResults(blastResults)
         
        # Step 5. Get the approriate flanking sequence for each anti-repeat candidate
        if len(crispr.antiRepeats) == 0:
            noPredictedTracr[protID] = 'No Anti-repeats'
            continue
        if operon.seq is None: hadToGetSeq += 1; operon.seq = fasta_index(operon.getFastaPointer(protID),'fasta')[protID]

        crispr.getAntiRepeatCandidates(open("tmp/possibleTracrs.fasta","w"), operon.getSeq())
        
        # Step 6. Look for termination signals
        res = system("~/bin/Arnold/erpin ~/bin/Arnold/rho-indep.epn tmp/possibleTracrs.fasta -1,4 -add 1 4 2 -pcw 3.0 -cutoff 100% >tmp/rhoInd.out")
        if res != 0:
            noPredictedTracr[protID] = 'Erpin Failed %i' % (res)
            continue
        
        # Step 7. Read the termination signals
        erpOut = ErpinOut()
        erpSols += len(erpOut.terminators)

        # Step 8. Get tracrRNA candidates with rho-ind signals
        numNewTracrs = crispr.getTracrRNA_Candidates(erpOut,possibleSol)
        if numNewTracrs == 0: noPredictedTracr[protID] = 'No terminators'
        
        # Keep track of how many solutions have been found so far and print
        totalSols += numNewTracrs
except: die = True

print("\nHad to getAssemblies:",hadToGetSeq)
possibleSol.close()
if not die:
    dump(casOperons, "pickles/%s_Operons.p" % gene)
    dump(noPredictedTracr, "pickles/%sWithNoPredictedTracr.p" % (gene))
print("\nErpin Solutions:", erpSols)
print("Found %i possible tracr solutions from %i assmeblies" % (totalSols,len(allAssemblies)-len(noPredictedTracr)))

casOperons = load(open("pickles/%s_Operons.p" % gene,"rb"))

totalTerms = 0
for protID in allAssemblies:
    if protID == "CP010309_modified_ORF105548":continue
    operon = casOperons.operons[casOperons.seqMap[protID]]
    crispr = operon.getCRISPR(protID)
    totalTerms += len(crispr.terminators)
print("DoubleCheck:",totalTerms)
