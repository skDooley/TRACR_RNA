import os
import pickle
from easyFunctions import Coordinate

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
        if (self.strand == other.strand):die;return self.location.start-other.location.start <= 10
        return self.location.start-other.location.end <= 10
    def __str__(self):return "%s %s %s" % (self.structName,str(self.location),self.confidence)

class ProcessInfernalFile(dict):
    def __init__(self,filePath):
        fh = open(filePath)
        line = ""
        while "Query" not in line: line = fh.readline()
        self.name = line.strip().split(" ")[-3]
        #Skip next 3 lines
        for times in range(4): line = fh.readline()
        while "inclusion threshold" not in line and "" != line:
            line = line.replace("  ", " ").replace("  ", " ").replace("  ", " ").replace("  ", " ").replace("  ", " ").strip()
            lineInfo = line.split(" ")
            newTracr = TRACR_RNA(self.name,lineInfo[5],lineInfo[6],lineInfo[7],lineInfo[8],lineInfo[2])
            try: self[newTracr.sourceName].append(newTracr)
            except: self[newTracr.sourceName] = [newTracr]
            line = fh.readline().strip()
        fh.close()

class Infernal_Results:
    def __init__(self):
        self.seqTracrs = {}
        self.seqMap = {}
        self.structMapping = {}
    def add(self,results):
        for seqID in results:
            for tracr in results[seqID]:
                try: self.seqTracrs[seqID][tracr.structName] = tracr
                except: self.seqTracrs[seqID] = {tracr.structName:tracr}
                try: self.seqMap[seqID].add(tracr.structName)
                except: self.seqMap[seqID] = set([tracr.structName])
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
        for structName, structureHit in self.structMapping.iteritems():
            yield structName, structureHit
    def tracrs(self,id):
        for tracr in self.seqTracrs[id]: yield tracr
        

############################### MAIN METHODS ###################################
def ProcessInfernal(roundNum):
    path = "conseqs%i/" %(roundNum)
    os.chdir(path)
    dirs = os.listdir("./")
    infernalResults = Infernal_Results()
    for dir in dirs:
        if ".py" in dir:continue
        os.chdir(dir)
        files = os.listdir("./")
        for file in files:
            if ".out" not in file:continue
            if os.stat(file).st_size == 0: 
                print ("Empty Output File: %s %s" % (file,str(os.stat(file).st_size)))
                continue
            results = ProcessInfernalFile(file)
            infernalResults.add(results)
        os.chdir("../")
    
    print ("Dumping","/mnt/research/germs/shane/transActRNA/data/pickles/InfernalResults0%i.p" % (roundNum))
    print("\tNumber of results: %i" % len(infernalResults) )       
    pickle.dump(infernalResults,open("/mnt/research/germs/shane/transActRNA/data/pickles/InfernalResults0%i.p" % (roundNum),"wb"))
    os.chdir("/mnt/research/germs/shane/transActRNA/data")
    return infernalResults
    
    
if __name__ == "__main__":  
    #import Infernal_Results
    ProcessInfernal(2)    
#    print "Number of assemblies covered:", len(existing.keys())
#    print "Number of assemblies w/more than one struct hit:",dupHitsInAssembly
