from pickle import load
from os import system,path
import sys
sys.path.append("../scripts/")
from HMMParser import *
goodDomIDS = load(open("pickles/%s_GoodDomainIDS.p","rb"))
hmm_parser = load(open("pickles/hmm_parsing_results.p","rb"))
assembliesWCrisprs = load(open("pickles/allAssemblyCRISPRs.p","rb")) 
counter = 0
protMapper = {}
start = int(len(hmm_parser.results)/4 + 15000)
print(start)
for result in hmm_parser.results[start:0:-1]:
    fullHmmPath = result.filePath
    baseID = fullHmmPath.replace("hmm/results/","").replace(".hmmout","")
    if path.exists("assemblies/assemblies_W_Cas9s/%s.fasta" % (baseID)):
        counter+=1
        if counter % 1000 ==0: print(counter,end=" ")
        continue
    for rec in parse(assembliesWCrisprs[baseID],"fasta"):
        if rec.id in goodDomIDS:
            system("ln -s %s assemblies/assemblies_W_Cas9s/%s.fasta" % (assembliesWCrisprs[baseID],baseID))
            counter+=1
            if counter % 1000 ==0: print(counter,end=" ")
            break