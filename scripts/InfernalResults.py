import os
from collections import Counter
from easyFunctions import Coordinate,RC
from glob import glob
from pickle import dump
from random import randint
from sys import argv

########################### Functions & Classes ################################       
class TRACR_RNA:
    def __init__(self,structName,source,start,end,strand,conf,tuncated):
        self.structName=structName
        self.sourceName = source
        self.location = Coordinate(int(start),int(end))
        self.confidence = conf
        self.strand = (strand == "+")
        self.overlapsBLAST = False
        self.truncated3Prime = ("3'" in tuncated)
        self.truncated5Prime = ("5'" in tuncated)
        self.aln = None
    def __hash__(self):
        if self.strand: return hash("%i_%i" % (self.location.start,self.location.end))
        return hash("%i_%i" % (self.location.end,self.location.start))
    def __eq__(self,other): 
        return self.location.start-other.location.start <= 10 and \
        self.location.end-other.location.end <= 10
    def __str__(self): 
        truncated = 'no'
        if self.truncated3Prime and self.truncated5Prime: truncated="5'/3'"
        elif self.truncated3Prime: truncated="3'"
        elif self.truncated5Prime: truncated="5'"
        return "%s %s %s %s" % (self.structName,str(self.location),self.confidence,truncated)
    def seq(self, fullSeq): 
        if not self.strand:return RC(str(fullSeq[self.location.start-1:self.location.end]).upper())
        return str(fullSeq[self.location.start-1:self.location.end]).upper()

class ProcessInfernalFile(dict):
    def __init__(self,filePath,recordTruncated=True):
        fh = open(filePath)
        line = ""
        counter = 0
        while "Query" not in line: 
            line = fh.readline()
        self.name = line.strip().split(" ")[-3]
        #Skip next 3 lines
        for times in range(4): line = fh.readline()
        while "inclusion threshold" not in line and "" != line and "Hit alignments" not in line:
            line = line.replace("  ", " ").replace("  ", " ").replace("  ", " ").replace("  ", " ").replace("  ", " ").strip()
            lineInfo = line.split(" ")
            newTracr = TRACR_RNA(self.name,lineInfo[5],lineInfo[6],lineInfo[7],lineInfo[8],lineInfo[2],lineInfo[-3])

            if newTracr.truncated3Prime and not recordTruncated:
                line = fh.readline().strip()
                continue
            try: self[newTracr.sourceName].append(newTracr)
            except: self[newTracr.sourceName] = [newTracr]
            line = fh.readline().strip()
        fh.close()

class Infernal_Results:
    def __init__(self,dumpResults=True):
        self.seqTracrs = {}
        self.seqMap = {}
        self.structMapping = {}
        self.dump = dumpResults
    def add(self,results):
        newSeqIDs = set(results.keys()).difference(self.seqTracrs.keys())
        for seqID in newSeqIDs: self.seqTracrs[seqID], self.seqMap[seqID] = {}, set()
        for seqID in results:
            for tracr in results[seqID]:
                try: self.seqTracrs[seqID][tracr.structName].append(tracr)
                except: self.seqTracrs[seqID][tracr.structName]=[tracr]
                self.seqMap[seqID].add(tracr.structName)
                try: self.structMapping[tracr.structName].add(seqID)
                except: self.structMapping[tracr.structName] = set([seqID])
    def __len__(self):return len(self.seqTracrs)
    def __contains__(self,item):return item in self.seqMap or item in self.structMapping
    def __iter__(self):
        for seqID in self.seqTracrs:yield seqID
    def iteritems(self):
        for seqID,tracrList in self.seqTracrs.iteritems():yield seqID,tracrList
    def __getitem__(self,item):
        if item in self.seqMap: return self.seqMap[item]
        elif item in self.structMapping: return self.structMapping[item]
    def structs(self):
        for structName, structureHit in self.structMapping.iteritems(): yield structName, structureHit
    def tracrs(self,id):
        for tracr in self.seqTracrs[id]: yield tracr


############################### MAIN METHODS ###################################
def ProcessInfernal(roundNum,gene,dumpResults=True,recordTruncated=True):
    path = "conseqs%i/" %(roundNum)
    os.chdir(path)
    # dirs = os.listdir("./")
    infernalResults = Infernal_Results()
    files = glob("*/*.out")
    for fName in files:
        if os.stat(fName).st_size == 0: 
            clName = fName[:fName.find('/')]
            print ("Empty Output File: %s %s %s" % (fName,str(os.stat(fName).st_size),clName))
            # os.system("mv %s/%s_Seqs.fasta ./%s.fasta" % (clName,clName,clName))
            # os.system("rm -rf %s" % (clName))
            continue
        results = ProcessInfernalFile(fName,recordTruncated)
        infernalResults.add(results)
    print("\tNumber of results: %i" % len(infernalResults) )
    if dumpResults:
        print ("Dumping","/mnt/research/germs/shane/transActRNA/data/pickles/%s_InfernalResults0%i.p" % (gene,roundNum))
        dump(infernalResults,open("/mnt/research/germs/shane/transActRNA/data/pickles/%s_InfernalResults0%i.p" % (gene,roundNum),"wb"))
    os.chdir("/mnt/research/germs/shane/transActRNA/data")
    return infernalResults
    
    
if __name__ == "__main__":  
    #import Infernal_Results
    ProcessInfernal(int(argv[1]),"Cas9")    
#    print "Number of assemblies covered:", len(existing.keys())
#    print "Number of assemblies w/more than one struct hit:",dupHitsInAssembly
