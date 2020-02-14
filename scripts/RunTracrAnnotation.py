import pickle
import sys
import tempfile
sys.path.append("scripts/")

from AnnotateTree import AnnotateTree
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqIO import index as fasta_index, parse, write
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat import MatrixInfo as matlist
from CRISPRtools import * #BLAST, MakeFasta, PilerCRReader, MinCEDReader
from glob import glob
from easyFunctions import Coordinate, BLAST_short
from os import chdir, path, stat, system
from random import randint

################################         Data          ################################
chdir("/mnt/research/germs/shane/transActRNA/data")
casRelatedAssemblies={}
for rec in parse(open("sequences/CasRelatedAssemblies_NoTracr.fa","rb"),"fasta"):casRelatedAssemblies[rec.id]=str(rec.seq)
#allCas9s = fasta_index("sequences/Cas9-Like-clustered-2.faa", "fasta")
#print "Number of Cas9s:", len(allCas9s)

################################    Custom Methods     ################################
r = lambda: randint(0,255)            
def color():return '#%02X%02X%02X' % (r(),r(),r())

def WriteSequence(tmpFasta,index,seq,rec):
    index+=1
    tmpFasta.write(">Seq_%i_%i_%i\n%s\n" % (min(rec.start,rec.end),max(rec.start,rec.end),index,seq[min(rec.start,rec.end):max(rec.start,rec.end)+500]))
    index+=1
    tmpFasta.write(">Seq_%i_%i_%i\n%s\n" % (min(rec.start,rec.end),max(rec.start,rec.end),index,seq[min(rec.start,rec.end)-500:max(rec.start,rec.end)]))
    return index


################################ Custom Data Structures ################################
class ErpinOut:
    def __init__(self, outfile="tmp/rhoInd.out", inputfile="tmp/possibleTracrs.fasta"):
        self.numRecords = 0
        self.terminators = []
        self.records={}
        with open(inputfile) as file:
            for rec in parse(file,"fasta"):
                self.numRecords += 1
                self.records[rec.id]=str(rec.seq)
        with open(outfile) as file:
            for i in range(9): file.readline()
            capture = False
            for line in file:
                if capture: 
                    line = line.strip().replace("  "," ").replace("  "," ").replace("  "," ")
                    self.terminators.append(RhoIndTerminator(seqName,line.split(" ")))
                capture = (">" in line)
                if capture: seqName = line.strip().replace(">","")
class RhoIndTerminator:
    def __init__(self,name,info):
        self.name = name
        self.strand = (info[0]=="FW")
        start,end = info[2].split("..")
        self.Rholocation = Coordinate(start,end)
        self.fwd = int(self.name[-1]) % 2 != 0
    def __str__(self): return "%s\t%s\t%s\t%s\t" % (self.name,str(self.strand),str(self.Rholocation),str(self.fwd))
    
################################ Custom Functions ################################
matrix = matlist.blosum62
gap_open = -8
gap_extend = -.8

def RC(seq): return str(Seq(seq).reverse_complement())
def scoreAlign(alignment):
    ref, frag, score, begin, end = alignment
    matches = 0
    for pos in range(len(ref)):
        if ref[pos] == frag[pos]:matches+=1
    return matches/float(len(frag.replace("-","")))         
def scoreAligns(aln1,aln2):
    score1, score2 = scoreAlign(aln1), scoreAlign(aln2)
    if score1 > score2: return aln1,score1*100
    else: return aln2, score2*100   
    
def alignSequences(refSeq,fragment):
    try: aln1 = pairwise2.align.globalds(refSeq, fragment, matrix, gap_open, gap_extend)[0]
    except: aln1 = None
    try: aln2 = pairwise2.align.globalds(refSeq, RC(fragment), matrix, gap_open, gap_extend)[0]
    except: aln2 = None
    if aln1 == None and aln2 == None: return None,0
    elif aln1 == None: top_aln = aln2
    elif aln2 == None: top_aln = aln1
    else: top_aln,alnScore = scoreAligns(aln1,aln2)
    if alnScore == None: alnScore = 0    
    if top_aln == None:print "Here"
    aln_probe, aln_arms, score, begin, end = top_aln
    return alnScore

def alignSequence(refSeq,fragment):
    aln1 = pairwise2.align.globalds(refSeq, fragment, matrix, gap_open, gap_extend)[0]
    try: aln1 = pairwise2.align.globalds(refSeq, fragment, matrix, gap_open, gap_extend)[0]
    except: aln1 = None
    try: aln2 = pairwise2.align.globalds(refSeq, RC(fragment), matrix, gap_open, gap_extend)[0]
    except: aln2 = None
    if aln1 == None and aln2 == None: return None,0
    elif aln1 == None: top_aln = aln2
    elif aln2 == None: top_aln = aln1
    else: top_aln,alnScore = scoreAligns(aln1,aln2)
    if alnScore == None: alnScore = 0    
    if top_aln == None:print "Here"
    aln_probe, aln_arms, score, begin, end = top_aln
    return '%s\t%% Matching %.2f%%\n\t%s' % (aln_probe, alnScore, aln_arms), alnScore

minCED_results  = pickle.load(open("pickles/MinCED_CRISPRS.p","rb"))
pilerCR_results = pickle.load(open("pickles/PilerCR_CRISPRS.p","rb"))

crRNALens, crisprLocations, noPredictedTracr, seqLenDist = [],[],[],[]
hasNRegionSet,crRNALens, possibleSolutions, newSolutions, index, counter = set(),[], set(), {}, 0, 0
nCounter, n_and_no_possibles = 0,0
for protID in casRelatedAssemblies:
    baseID = protID[:protID.find("_ORF")]
    if protID not in pilerCR_results and protID not in minCED_results: 
        die
        continue #Missing from both
    try: 
        if protID in pilerCR_results: assemblyLoci = pilerCR_results[protID].values()
        else:
            assemblyLoci = pilerCR_results[baseID].values() #Looking only at the most abundant crispr
        for locus in assemblyLoci:
            if locus.name == protID: break
        if locus.name != baseID: die
    except: 
        try: assemblyLoci = minCED_results[baseID].values()
        except: assemblyLoci = minCED_results[protID].values()
                 
        for locus in assemblyLoci:
            if locus.name == baseID: break
        if locus.name != baseID: die

    write([SeqRecord(id=protID,description='',seq=Seq(locus.consensusRepeat[0]))],"tmp/consFasta.fa","fasta")
    blast_results = parseSingleBLAST(BLAST_short("tmp/consFasta.fa", "assemblies/%s.fa" % (protID), "tmp/consBlast.xml"))
    
    #Step 2b
    boundaryHits,crRNALen = locus.checkArrayBoundaries(blast_results)
    crRNALens.append(crRNALen)
    locus.annotate(casRelatedAssemblies[protID], "assemblies/%s.fa" % (protID))
    if locus.hasNRegion: 
        nCounter += 1
        hasNRegionSet.add(protID)
    
    #Step 2c
    tmpFasta = open("tmp/possibleTracrs.fasta","w")
    
    ##Step 2d
    terminusSeqs, index = locus.getTerminusRepeats(casRelatedAssemblies[protID],index)
    tmpFasta.write(terminusSeqs)
    
    #Collect some stats
    startCoord = locus.repeatCoords[0]
    endCoord = locus.repeatCoords[-1]
    seqlen = len(casRelatedAssemblies[protID])
    seqLenDist.append(seqlen)
    minDist = min(seqlen-startCoord.start,seqlen-startCoord.end,seqlen-endCoord.start,seqlen-endCoord.end)
    minDist =  min(minDist,startCoord.start,startCoord.end,endCoord.start,endCoord.end)
    crisprLocations.append(minDist)

    ##Step 2e
    for rec in boundaryHits: index = WriteSequence(tmpFasta,index,casRelatedAssemblies[protID],rec)
    
    ##Step 2f
    tmpFasta.close()
    cmd = "~/bin/Arnold/erpin ~/bin/Arnold/rho-indep.epn tmp/possibleTracrs.fasta -1,4 -add 1 4 2 -pcw 1 -cutoff 100% >tmp/rhoInd.out"
    system(cmd)
    
    ##Setp 2j
    erpOut = ErpinOut()
    if len(erpOut.terminators)==0: 
        print "\nThis ref has nothing: %s\n" % protID
        noPredictedTracr.append(protID)
        if locus.hasNRegion: n_and_no_possibles +=1
    else: 
        newSolutions[protID] = erpOut
        repeatIndex = 0

pickle.dump(newSolutions,open("pickles/predictedTracrRNAs.p","wb"))
print "\n\nNumber of proteins with large N region and no terminator / with total large N region: %i / %i" % (n_and_no_possibles,nCounter)
print "%i possible solutions in %s references. Nothing found in %i references" % (index-1, len(newSolutions),len(noPredictedTracr))
