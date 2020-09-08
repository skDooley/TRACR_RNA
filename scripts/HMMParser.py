from Bio.SeqIO import parse
import numpy as np
class SamplesDict:
    def __init__(self):
        self.domains = set()
        self.domainMap = {}
        self.samples = {}
    def __setitem__(self,key, value):
        self.domains.add(value)
        try:self.samples[key][value.domain_name]=value
        except: self.samples[key] = {value.domain_name:value}
        try:self.domainMap[value.domain_name].add(key)
        except:self.domainMap[value.domain_name] =set([key])
    def __getitem__(self,key): return self.domainMap[key]
    def iteritems(self):return self.domainMap.iteritems()

class DomainHit:
    def __init__(self,fields):
        self.hit = fields[0]
        self.domain_name = fields[3]
        self.tlen = int(fields[2])
        self.start = int(fields[19])
        self.end = int(fields[20])
        self.dend = self.tlen-self.end
    def __hash__(self): return hash(self.domain_name)
    
def outliers(data,stdsAwayFromMean):
    return data[abs(data - np.mean(data)) > stdsAwayFromMean * data.std()]

class HMM_Parser:
    def __init__(self,fname=""): 
        self.results = {}; self.protIDs = {}
        if fname != "":
            self.parse(fname)
            self.results = self.results[self.basePath].proteins

    def parse(self,filePath): 
        hresult = HMM_Result(filePath)
        if hresult.has_hits: 
            basePath = filePath[filePath.rfind("/")+1:].replace(".hmmout","")
            self.basePath = basePath
            self.results[basePath] = hresult
            for recID, prot in hresult: self.protIDs[recID] = basePath
    def __len__(self): return len(self.results)
    def __iter__(self): return iter(self.results.items())
    def __getitem__(self, protName): return self.results[self.protIDs[protName]].proteins[protName]
    def __contains__(self, protName): return protName in self.results

from os import path
class HMM_Result:
    def __init__(self,filePath):
        hits = []
        for line in open(filePath):
            if ">> " in line: hits.append(line.strip().replace(">> ","").split(" ")[0])
        self.has_hits, self.hits, self.proteins = (len(hits)>0), set(hits), {}
        if self.has_hits: self.__get_proteins__(filePath)
        self.filePath=filePath
    def __get_proteins__(self,filePath):
        # if not path.exists(filePath.replace(".hmmout",".orfs")): return
        orfFile = "sequences/orfs/" + filePath[filePath.rfind("/"):]
        orfFile=orfFile.replace("//","/").replace(".hmmout",".orfs")
        if not path.exists(orfFile): 
            print("Missing:",orfFile)
            return
        for rec in parse(orfFile,"fasta"): #Assumes all hmm_results end with .hmmout and their corresponding orfs are in the same place but have a .orfs extension
            if rec.id in self.hits: self.proteins[rec.id]=rec
    def __iter__(self): return iter(self.proteins.items())
