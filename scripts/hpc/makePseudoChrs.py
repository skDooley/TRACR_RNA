from Bio.SeqIO import write
from pickle import load
from os import listdir, path, system
from sys import path as spath
spath.append("../scripts/")
from HMMParser import *

gene = "Cas9"
assemblyDir = "assemblies/assemblies_W_%s/" % (gene)
cas9Assemblies = listdir(assemblyDir)
goodDomIDS = load(open("pickles/%s_GoodDomainIDS.p"       % (gene), "rb"))
goodDomMap = load(open("pickles/%s_GoodDomMap.p"          % (gene), "rb"))
hmm_parser = load(open("pickles/%s_HMM_Parsing_Results.p" % (gene), "rb"))

print("All loaded")
#Copy unique nucleotide sequence
from Bio.SeqIO import index
nukSeqHash,protSeqHash = set(), set()
alreadyGotIt,count = 0,0
for assembly in cas9Assemblies:
    baseID = assembly[:-6]
    allAssemblySeqs = index(assemblyDir+assembly,"fasta")
    overlap = goodDomIDS.intersection(allAssemblySeqs.keys())
    for recID in overlap:
        seq = str(allAssemblySeqs[recID].seq).upper()
        if seq in nukSeqHash and len(goodDomMap[recID]) == 1:
            alreadyGotIt += 1
            continue
        nukSeqHash.add(seq)
        #There may be more than 1 protein on the pseudochromosome, save both as separate files  
        if len(goodDomMap[recID])>1: print("%i Cas9s on %s %s" %(len(goodDomMap[recID]), recID, baseID))
        for orfID in goodDomMap[recID]:
#                 protSeq = str(hmm_parser.results[baseID].proteins[orfID].seq).upper()
            with open("assemblies/pseudoChromos/%s.fasta" % (orfID),"w") as fh:
                write(allAssemblySeqs[recID],fh,"fasta")
                count +=1
        if count % 1000 == 0:print(count,end=" ")