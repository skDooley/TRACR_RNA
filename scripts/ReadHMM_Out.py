import re
from sys import path as spath
spath.append("scripts/")
import numpy as np
from Bio.SeqIO import write, parse
from collections import Counter
from CRISPRtools import *
from easyFunctions import *
from HMMParser import *
from IPython.display import display, HTML
from matplotlib import pyplot as plt
from os import chdir, listdir,system
from pandas import Series
from pickle import load #Dump has been wrapped aournd a function from easyFunctions

chdir("/mnt/research/germs/shane/transActRNA/data")
gene = "Cas9"
geneProfile = "proteins/DiverseCas9s.faa" #File containing the seed proteins for the hmmsearch
geneFile = "proteins/%s.faa" % gene
alnName  = "alignments/%s.aln" % gene
hmmName  = "hmm/%s.hmm" % gene
hmmResultsDir = "hmm/results"
baseDbsDir = "/mnt/research/germs/shane/databases/assemblies/"
refDatabases = ["NCBI/refseq/bacteria","NCBI/refseq/archaea","NCBI/genbank/bacteria","NCBI/genbank/archaea","PATRIC2/fastas"]


#1. Read the domain search report
spacePat = re.compile( r'^[\t ]*$')
hits=SamplesDict()
print("Reading Table Definitions")
for line in open("hmm/Cas9-Like_phi.faa.domtbl"):
    if line.startswith( '#' ):continue
    fields = re.split( r'(\[[^\[]*[^\S]+[^\]]*\]|[^\t ]+)', line.strip() )
    fields = DomainHit(list(filter(lambda i: not spacePat.search(i), fields)))
    hits[fields.hit]=fields
print(len(hits.samples))
dump(hits,"pickles/%s_HMM_DOMAIN_Search_Results.p" % (gene))


#2. Generate the profile for all of the proteins
hasAllDomains = hits["RuvC_1_Cas9"].intersection(hits["RuvC_2_Cas9"].intersection(hits["RuvC_3_Cas9"].intersection(hits["HNH_Cas"])))
allSamples = set(hits.samples.keys())
nSamples = len(allSamples)
noRuvC1 = allSamples.difference(hits["RuvC_1_Cas9"])
noRuvC2 = allSamples.difference(hits["RuvC_2_Cas9"])
noRuvC3 = allSamples.difference(hits["RuvC_3_Cas9"])
noHNH   = allSamples.difference(hits["HNH_Cas"])

htmlString ="""<table align='left'>
    """\
    "<tr style='background-color:#9B7E46;color:white'><td>Number of Proteins:</td><td>%i</td></tr>" % (nSamples) +\
    "<tr style='background-color:#373D20;color:white'><td>Has all domains</td><td>%i</td></tr>" % (len(hasAllDomains)) +\
    "<tr style='background-color:#BCBD8B;color:white'><td>Number of sequences with no detected RuvCI domain</td><td>%i</td></tr>" % (len(noRuvC1)) +\
    "<tr style='background-color:#717744;color:white'><td>Number of sequences with no detected RuvCII domain</td><td>%i</td></tr>" % (len(noRuvC2)) +\
    "<tr style='background-color:#766153;color:white'><td>Number of sequences with no detected RuvCIII domain</td><td>%i</td></tr>" % (len(noRuvC3)) +\
    "<tr style='background-color:#F4B266;color:white'><td>Number of sequences with no detected HnH</td><td>%i</td></tr>" % (len(noHNH)) +\
    """
    </table>
"""

print(htmlString,"\n\n\n")

ruvC1Coords = []
ruvC3Coords = []
samples = []
difs = []
for sample in hasAllDomains:
    ruc1Start = hits.samples[sample]["RuvC_1_Cas9"].start
    ruc3End = hits.samples[sample]["RuvC_3_Cas9"].dend
    ruvC1Coords.append(ruc1Start)
    ruvC3Coords.append(ruc3End)
    difs.append(ruc3End-ruc1Start)
    samples.append(sample)
    
ruvC1Dists = Series(ruvC1Coords,index=samples)
ruvC3Dists = Series(ruvC3Coords,index=samples)
distBetween = Series(difs,index=samples)

print ("\nMean distance from RuvC_1 to start of protein: %.0f std: %.0f" % (ruvC1Dists.mean(), ruvC1Dists.std()))

print()

start_outliers = set(ruvC1Dists[ruvC1Dists > ruvC1Dists.mean()+8*ruvC1Dists.std()].index) #RuvC1 outlier
end_outliers = set(ruvC3Dists[ruvC3Dists > ruvC3Dists.mean()+5*ruvC3Dists.std()].index) #RuvC3 outlier
hasGoodDomains = hasAllDomains.difference(start_outliers.union(end_outliers))
htmlString ="""<table align='left'>
    """\
    "<tr style='background-color:#373D20;color:white'><td>No outlier domains</td><td>%i</td></tr>" % (len(hasGoodDomains)) +\
    "<tr style='background-color:#BCBD8B;color:white'><td>Number of RuvC1 outliers</td><td>%i</td></tr>" % (len(start_outliers)) +\
    "<tr style='background-color:#717744;color:white'><td>Number of RuvC3 outliers</td><td>%i</td></tr>" % (len(end_outliers)) +\
    "<tr style='background-color:#766153;color:white'><td>Outlier Intersection</td><td>%i</td></tr>" % (len(start_outliers.intersection(end_outliers))) +\
    """
    </table>
"""

print(htmlString+'\n\n')

#Checked domains, removed outliers, now remove chrs that have more than 1 Cas9
corrected = {}
remove,baseMap = set(),{}
print("Before Checking chromosomes:",len(hasGoodDomains))
for orfName in hasGoodDomains:
    baseCHR = orfName[:orfName.rfind("_")]
    if baseCHR in corrected:
        try: 
            baseMap[baseCHR].add(orfName)
            baseMap[baseCHR].add(corrected[baseCHR])
        except: baseMap[baseCHR] =set([corrected[baseCHR],orfName])
        remove.add(orfName)
        remove.add(corrected[baseCHR])
    corrected[baseCHR] = orfName
hasGoodDomains = hasGoodDomains.difference(remove)
print("After Checking:",len(hasGoodDomains))

dump(hasGoodDomains,"pickles/%s_CorrectDomains.p" % gene)

#3. Done Profiling Domains now lets create fasta files
casOperons = load(open("pickles/%s_Operons.p" % gene,"rb"))
casOperons.getRepSeqs(hasGoodDomains,"proteins/All_%s-Like-filtered.faa" %(gene),"proteins/All_%s-Like.faa" % (gene)) 
casOperons.getRepSeqs(hasGoodDomains,"assemblies/All_%s-Like-filtered.fasta" %(gene),"assemblies/All_%s_Unique_Assemblies.fasta" % (gene))


#4. Now lets clean up a bit
casAAs = dict(index("proteins/All_%s-Like-filtered.faa" %(gene),"fasta"))
unUsed = set(casOperons.seqMap.protToAsm).difference(casAAs)

pres,absnt,hasSeq = 0,0,0
for protID in unUsed:
    try:
        operon = casOperons.operons[casOperons.seqMap[protID]]
        hasSeq+= int(operon.seq is not None)
        pres +=1
        del casOperons.operons[casOperons.seqMap[protID]]
    except:
        absnt+=1
print("After cleaning up we removed %i sequences from %i operons" % (hasSeq,len(unUsed)))

print("Failed:",absnt)
from pickle import dump
dump(casOperons,open('pickles/Cas9_Operons.p','wb'))
print("Done")




# Teasing out on the same chr but uniq vs duplication event (ie found orfs)
from Bio.SeqIO import index as fasta_index
seqs = fasta_index("proteins/All_%s-Like.faa" % (gene),"fasta")
prots = {}
print(len(remove))
for orfID in remove: 
    seq = str(seqs[orfID].seq) 
    try: prots[seq].add(orfID)
    except: prots[seq] = set([orfID])
print(len(prots))

















