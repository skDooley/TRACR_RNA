import re
import sys
sys.path.append("/mnt/research/germs/shane/transActRNA/scripts/")

from Bio.SeqIO import write, parse, index
from CRISPRtools import *
from easyFunctions import *
from HMMParser import *
from os import chdir, listdir, system
from pickle import load #Dump has been wrapped aournd a function from easyFunctions

chdir("/mnt/research/germs/shane/transActRNA/data")
gene = "Cas9"
hmmResultsDir = "hmm/results"
# crisprFiles = load(open("pickles/CRISPRs.p","rb")) 
# casOperons = CasOperons(gene)
# casOperons.hasCas9(hmmResultsDir,crisprFiles)
casOperons = load(open("pickles/Cas9_Operons_HMM.p",'rb'))

#Get unique chrs and the proteins they are associated with
allCasAsmFile = "assemblies/All_%s_Unique_Assemblies.fasta" % (gene)
allCasAAsFile = "proteins/All_%s-Like.faa" % (gene)
casOperons.uniqueNukeSeqs(allCasAsmFile,allCasAAsFile) # Calls dump when it finishes

# Launch the domain search for the faa file created above
system("sbatch /mnt/research/germs/shane/transActRNA/scripts/hpc/DomainSearch.sb")
casAAs = dict(index(allCasAAsFile,"fasta"))
unUsed = set(casOperons.seqMap.protToAsm).difference(casAAs)
deletedOperons = {}
pres,absnt,hasSeq = 0,0,0
for protID in unUsed:
    try:
        operon = casOperons.operons[casOperons.seqMap[protID]]
        deletedOperons[protID] = operon
    except: absnt+=1

from pickle import dump
dump(deletedOperons,open("/mnt/research/germs/shane/transActRNA/data/pickles/DeletedOperons.p","wb"))
del deletedOperons

for protID in unUsed:
    try:
        operon = casOperons.operons[casOperons.seqMap[protID]]
        hasSeq+= int(operon.seq is not None)
        pres +=1
        del casOperons.operons[casOperons.seqMap[protID]]
    except: absnt+=1

print("After cleaning up we removed %i sequences from %i operons" % (hasSeq,len(unUsed)))
print("Failed:",absnt)
dump(casOperons,open('/mnt/research/germs/shane/transActRNA/data/pickles/Cas9_Operons_SeqFilter.p','wb'))
print("Done")

#casOperons = load(open('pickles/Cas9_Operons.p','rb'))
