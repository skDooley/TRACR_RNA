
from Bio import AlignIO, pairwise2
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SubsMat import MatrixInfo as matlist
from collections import Counter
from random import randint
import numpy as np
import os
from pickle import dump as pcklDump, load
import random
import sys
import tempfile

IUPACS = {
'[A/G]':'R','[C/T]':'Y','[A/T]':'W','[G/A]':'R','[T/C]':'Y','[T/A]':'W',
'[C/G]':'S','[A/C]':'M','[G/T]':'K','[G/C]':'S','[C/A]':'M','[T/G]':'K',
'[C/G/T]':'B','[C/T/G]':'B','[G/C/T]':'B','[G/T/C]':'B','[T/C/G]':'B','[T/G/C]':'B', 
'[A/G/T]':'D','[A/T/G]':'D','[G/A/T]':'D', '[G/T/A]':'D','[T/A/G]':'D','[T/G/A]':'D', 
'[A/C/T]':'H','[A/T/C]':'H','[C/A/T]':'H','[C/T/A]':'H','[T/A/C]':'H','[T/C/A]':'H',
'[A/G/C]':'V','[A/C/G]':'V','[G/A/C]':'V', '[G/C/A]':'V','[C/A/G]':'V','[C/G/A]':'V', 
'[A/C/T/G]':'N' }

def dump(dmpObj,fname):
    print ("\nDumping "+ fname)
    pcklDump(dmpObj,open(fname,"wb"))

class Coordinate:
    def __init__(self,start,end): self.start = int(start); self.end = int(end); self.strand = (self.end > self.start)
    def contains(self,other): return self.start <= other.start and other.end <= self.end
    def __eq__(self,other): return (self.start == other.start and self.end == other.end) or \
                                   (self.start == other.end and self.end == other.start)
        
    def overlaps(self,other):
        #Equal if they overlap or 1 contains the other
        if self.start>self.end:
            startBetween = self.end <= other.start and other.start <= self.start 
            endBetween = self.end <= other.end and other.end <= self.start
        else:
            startBetween = self.start <= other.start and other.start <= self.end 
            endBetween = self.start <= other.end and other.end <= self.end
        contains = self.contains(other) or other.contains(self)
        return startBetween or endBetween or contains

    def distance(self,other):
        if self.overlaps(other):return 0
        return min(abs(self.start-other.start),
                   abs(self.start-other.end),
                   abs(self.end-other.start),
                   abs(self.end-other.end))
    def __ne__(self,other): return not self.__eq__(other)
    def __hash__(self): return hash(self.start)
    def __str__(self): 
        dir="+"
        #if not self.strand:dir="-"
        return "[%s\t%s]%s" % (min(self.start,self.end),max(self.start,self.end),dir)
    def __sub__(self,other):return abs(self.start-other.start)+abs(self.end-other.end)
    def __len__(self):return max(self.start,self.end)-min(self.start,self.end)

class MakeFasta:
    def __init__(self,fileName): 
        if type(fileName) is str: self.fileHandle = open(fileName,"w")
        else: self.name = fileName.name; self.fileHandle = fileName
    def write(self,id,seq):self.fileHandle.write(">%s\n%s\n" %(id,seq))
    def seek(self,num):self.fileHandle.seek(num)
    def close(self): self.fileHandle.close()
    def die(self):os.system("rm "+self.fileHandle.name)

def createBLASTdb(refFileName):
    outputDBName = "blastDBs/" + refFileName[refFileName.rfind("/")+1:]
    cmd = "makeblastdb -in %s -dbtype \"nucl\" -out %s >/dev/null" %(refFileName,outputDBName)
    runStatus = os.system(cmd)
#    print cmd,runStatus
    return outputDBName,runStatus

def BLAST_short(fragmentFile,sequenceFile,outputName):
    coords = []
    sequenceFileDB,runStatus = createBLASTdb(sequenceFile)
    if runStatus == 0:
        NcbiblastnCommandline(cmd='blastn', out=outputName, outfmt=5, query=fragmentFile, strand="both", dust="no", db=sequenceFileDB, evalue=0.001, word_size=7)() #reward=1,penalty=-3,
        return NCBIXML.parse(open(outputName,'r'))

def parseSingleBLAST(results): # For single blast ID because the results are all the same subject
        coords = []
        for result in results:
            for alignment in result.alignments:
                for hsp in alignment.hsps: 
                    if hsp.strand[0] is not None: print( hsp.strand+"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW")
#                    print (hsp)
#                    if hsp.expect < 0.00005:
                    coords.append(Coordinate(hsp.sbjct_start-1,hsp.sbjct_end))
#                    coords.append([hsp.sbjct_start-1,hsp.sbjct_end])
        coords.sort(key=lambda x: x.start)
        return coords

def alignSequence(refSeq,fragment):
    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -.8
    probe = refSeq 
    ext_plus_lig = fragment
    alns = pairwise2.align.globalds(probe, ext_plus_lig, matrix, gap_open, gap_extend)
    if len(alns) == 0: return
    top_aln = alns[0]
    aln_probe, aln_arms, score, begin, end = top_aln
    return aln_probe,aln_arms

def alignSeqs(sequenceDict):
    counter = 0
    newFasta = MakeFasta(tempfile.NamedTemporaryFile(delete=False))
    for ID,seqDict in sequenceDict.iteritems(): 
        if type(seqDict) == type(1): seqDict = sequenceDict
        for seq,count in seqDict.iteritems():
            for times in range(count):
                newFasta.write("Seq%i" % counter,seq)
                counter +=1
    newFasta.close()        
    clustalw_cline = ClustalwCommandline("dooleysh/bin/clustalW-2.1/clustalw", infile=newFasta.fileHandle.name)    
    out, err = clustalw_cline()
    alignments = AlignIO.read(newFasta.fileHandle.name +".aln", "clustal") 
    os.system("rm %s %s.aln %s.dnd" % (newFasta.fileHandle.name,newFasta.fileHandle.name,newFasta.fileHandle.name))
    return alignments

def calc_consensus_seq(sequenceList,useIUPACs=False):
    consensus = ""
    if useIUPACs:
        for pos in range(len(sequenceList[0])):
            positionalBases = []
            for seq in sequenceList: positionalBases.append(seq[pos])
            baseCounter = Counter(positionalBases).most_common(6) #6 for the 4 bases, 1 IUPAC and a gap char
            bases = {}
            for base,count in baseCounter:bases[base]=count
            if "N" in bases and max(bases.values())== bases["N"]:
                print ("Does this ever happen?")
                for seq in sequenceList: print (seq)
                sys.exit(0)
            if "N" in bases: del bases["N"]
            if "-" in bases:
                percentGaps = (bases["-"]/float(len(positionalBases)))*100
                if percentGaps <= 25.0 :del bases['-'] #Assumption is there is only a small percentaget then ignore it
            if "-" not in bases and len(bases)>1:
                combinedBases = "[%s]" % ("/".join(bases))
#                if combinedBases == "[A/C/T/G]": print (bases)
                consensus += IUPACS[combinedBases]
            elif "-" not in bases: consensus += baseCounter[0][0]
            elif percentGaps >= 75.0: pass #Asssumption, if there are a greater number of - then bases skip the position
            else:#Assumption: This condition will occur if there are "-" and but less than 75% and more than 25%
                consensus = sequenceList[0].replace("-","")
                break
    else:
        for pos in range(len(sequenceList[0])):
            positionalBases = []
            try:
                for seq in sequenceList: positionalBases.append(seq[pos])
            except: return sequenceList[0]
            baseCounter = Counter(positionalBases).most_common(3) #For now designed to blow up until we see the condition of have more than 2 matching alleles with the same count
            if len(baseCounter)>3 and (baseCounter[0][1] == baseCounter[1][1]) and (baseCounter[1][1] == baseCounter[2][1]) and baseCounter[2][1] == baseCounter[3][1]:stop; consensus += baseCounter[random.randint(0,3)][0]
            elif len(baseCounter)>2 and (baseCounter[0][1] == baseCounter[1][1]) and (baseCounter[1][1] == baseCounter[2][1]):stop; consensus += baseCounter[random.randint(0,2)][0]
            elif len(baseCounter)>1 and baseCounter[0][1] == baseCounter[1][1]:consensus += baseCounter[random.randint(0,1)][0]
            elif baseCounter[0][0] != "-": consensus += baseCounter[0][0]
    if "-" in consensus: print (consensus)
    return consensus

def CalcPercIdent(seqA,seqB):
    seq1,seq2 = alignSequence(seqA,seqB)
    matches = 0
    for index in range(len(seq1)):
        if seq1[index] == seq2[index]:matches += 1
    return float(matches)/len(seqA)*100

def outliers(data,stdsAwayFromMean): return data[abs(data - np.mean(data)) > stdsAwayFromMean * np.std(data)]

r = lambda: randint(0,255)            
def color():return '#%02X%02X%02X' % (r(),r(),r())

class ClusterDict(dict):
    def __init__(self):
        self.revMap = {}
        self.curID = ""
    def newCluster(self,ID): 
        self[ID] = Cluster(ID)
        self.curID = ID
    def add(self,line):
        seqID = self[self.curID].add(line)
        self.revMap[seqID] = self.curID
    def filter(self,cutoff=100):
        for id, cluster in self.iteritems():
            for seqID, perc in cluster.members.iteritems():
                if perc < cutoff: yield seqID, cluster.dirs[seqID]
            

class Cluster():
    def __init__(self,clusterName):
        self.name = clusterName
        self.members = {} #{SeqID:perc}   The percent identity to the representative sequence
        self.dirs = {}    #{SeqID:strand} Sequence direction relative to the representative sequence
    def __len__(self): return len(self.members)
    def add(self,line):
        seqID = line[line.find(">")+1:line.find("...")] 
        percIndex = line.find("%")
        if (percIndex != -1):
            perc = float(line[line.find("/")+1:percIndex])
            direction = (line[line.find("/")-1:line.find("/")] == "+")
        else:
            perc = 0.0
            direction = True
        self.dirs[seqID] = direction
        self.members[seqID]=perc
        return seqID
    def __iter__(self): return self.members.__iter__()

def processClusterFile(tracrFile):
    clusterResults = ClusterDict()
    for line in open(tracrFile):
        if ">Cluster" in line: clusterResults.newCluster(line.strip().replace(">",""))
        else: clusterResults.add(line)
#            sequenceID = line[line.find(">")+1:line.find(".")] 
#            seqID = seqID[:seqID.rfind("_")]
#            try:allClusters[clusterID].add(seqID)
#            except:allClusters[clusterID]=set([seqID])
#            try:allClusterSeqIDs[seqID].add(clusterID)
#            except:allClusterSeqIDs[seqID]=set([clusterID])
    return clusterResults

def generateTreeColors(allClusters,cutoff=5):
    remove, colors, TreeColors, clusterDist=set(), {}, {}, []
    for cluster, ids in allClusters.iteritems():
        clusterDist.append(len(ids))
        if len(ids) < cutoff: remove.add(cluster)
        else:
            colors[cluster]=color()
            for id in ids:TreeColors[id.replace(".","_")] = colors[cluster]
    for id in remove: del allClusters[id]
    print( "\tTotal number of clusters: %i" % len(allClusters))
    print ("\tNumber of nodes covered: %i" % len(TreeColors))
    
    return TreeColors, clusterDist


#looking through gbk files
# files = os.listdir(kozak_path)
# nameCounter = {"KOZAK01":0,"KOZAK06":0}
# for filename in files:
#     if ".gb" not in filename:continue
# #     print filename
#     for rec in parse(open(kozak_path+filename),"genbank"):
#         for feature in rec.features:
#             if "KOZAK" in feature.qualifiers['label'][0]:
                
#                 fh=open(kozak_path+feature.qualifiers['label'][0]+"_"+filename.replace(".gb","")+".fa","w")
                
#                 seq = str(rec.seq)
#                 rna = seq[feature.location.start-50:feature.location.end+50].replace("T","U")
#                 fh.write(">%s_%s\n%s\n" % (feature.qualifiers['label'][0],filename,rna))
#                 fh.close()
#                 break
#         print filename,feature.qualifiers['label'][0]
