import os
from collections import Counter
from easyFunctions import Coordinate
from glob import glob
from pickle import dump
from random import randint


########################### Functions & Classes ################################       
class TRACR_RNA:
    def __init__(self,structName,source,start,end,strand,conf):
        self.structName=structName
        self.sourceName = source
        self.location = Coordinate(start,end)
        self.confidence = conf
        self.strand = strand == "+"
        self.overlapsBLAST = False
        self.aln = None
    def __hash__(self):
        if self.strand: return hash("%i_%i" % (self.location.start,self.location.end))
        return hash("%i_%i" % (self.location.end,self.location.start))
    def __eq__(self,other): 
        return self.location.start-other.location.start <= 10 and \
        self.location.end-other.location.end <= 10
    def __str__(self): return "%s %s %s" % (self.structName,str(self.location),self.confidence)
    def seq(self, fullSeq): return str(fullSeq[self.location.start:self.location.end]).upper()
class ProcessInfernalFile(dict):
    def __init__(self,filePath):
        fh = open(filePath)
        line = ""
        counter = 0
        while "Query" not in line: 
            counter +=1
            line = fh.readline()
            if counter == 10000:
                print(filePath)
                die
        self.name = line.strip().split(" ")[-3]
        #Skip next 3 lines
        counter +=1
        for times in range(4): line = fh.readline()
        while "inclusion threshold" not in line and "" != line:
            line = line.replace("  ", " ").replace("  ", " ").replace("  ", " ").replace("  ", " ").replace("  ", " ").strip()
            lineInfo = line.split(" ")
            newTracr = TRACR_RNA(self.name,lineInfo[5],lineInfo[6],lineInfo[7],lineInfo[8],lineInfo[2])
            try: self[newTracr.sourceName].append(newTracr)
            except: self[newTracr.sourceName] = [newTracr]
            line = fh.readline().strip()
            if counter == 10000:
                print(filePath)
                die
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
r = lambda: randint(0,255)            
def genColor(): return '#%02X%02X%02X' % (r(),r(),r())

def generateColors(infernalResults,casAssemblies,tRound,gene):
    clusterCounter,treeMap,revMap, missingIDs = {}, {}, {}, set()
    for struct,ids in infernalResults.structMapping.items(): clusterCounter[struct] = len(ids)
    # for protID in casAssemblies:
    #     # baseID = protID[:protID.rfind("_")]
    #     if protID not in infernalResults.seqTracrs: missingIDs.add(protID); continue
    #     biggestCluster =  findBiggestCluster(infernalResults.seqTracrs[protID].keys(),clusterCounter)
    #     treeMap[protID] = biggestCluster
    #     try: revMap[biggestCluster].add(protID)
    #     except: revMap[biggestCluster] = set([protID])
    clusterCounter = Counter(clusterCounter)
    usedIDs = set()
    treeColors,colors = {},{}
    for cluster, protIDs in revMap.items():
        colors[cluster]=genColor()
        for protID in protIDs: treeColors[protID] = colors[cluster]
    dump(treeColors,open("pickles/%s_StructureResultsTreeColors_Round%i.p" % (gene,tRound),'wb'))

    dump(colors,    open("pickles/%s_StructureResultsColors_Round%i.p"     % (gene,tRound),'wb'))
    print("In round %i, %i assmebly had a predicted tracr (%.1f %%) with %i tracrRNA Structures" % 
         (tRound,len(treeMap),len(treeMap)/float(len(casAssemblies))*100,len(revMap)))
    print("\t%i still remaining" % (len(missingIDs)))     
def findBiggestCluster(clusters, clusterMap): return Counter(dict((k, clusterMap[k]) for k in (clusters))).most_common(1)[0][0] 


############################### MAIN METHODS ###################################
def ProcessInfernal(roundNum,gene,dumpResults=True):
    path = "conseqs%i/" %(roundNum)
    os.chdir(path)
    # dirs = os.listdir("./")
    infernalResults = Infernal_Results()
    files = glob("*/*.out")
    for fName in files:
        if os.stat(fName).st_size == 0: 
            print ("Empty Output File: %s %s" % (fName,str(os.stat(fName).st_size)))
            continue
        results = ProcessInfernalFile(fName)
        infernalResults.add(results)

    # for dir in dirs:
    #     if ".py" in dir:continue
    #     os.chdir(dir)
    #     files = os.listdir("./")
    #     for file in files:
    #         if ".out" not in file:continue
    #         if os.stat(file).st_size == 0: 
    #             print ("Empty Output File: %s %s" % (file,str(os.stat(file).st_size)))
    #             continue
    #         results = ProcessInfernalFile(file)
    #         infernalResults.add(results)
    #     os.chdir("../")
    if dumpResults:
        print ("Dumping","/mnt/research/germs/shane/transActRNA/data/pickles/%s_InfernalResults0%i.p" % (gene,roundNum))
        print("\tNumber of results: %i" % len(infernalResults) )       
        dump(infernalResults,open("/mnt/research/germs/shane/transActRNA/data/pickles/%s_InfernalResults0%i.p" % (gene,roundNum),"wb"))
    os.chdir("/mnt/research/germs/shane/transActRNA/data")
    return infernalResults
    
    
if __name__ == "__main__":  
    #import Infernal_Results
    ProcessInfernal(1,"Cas9")    
#    print "Number of assemblies covered:", len(existing.keys())
#    print "Number of assemblies w/more than one struct hit:",dupHitsInAssembly
