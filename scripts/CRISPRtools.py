'''
 Tools to find, annotate, and annalyze CRISPR Loci
'''
#import sys
import tempfile
#from Bio.Alphabet import generic_dna
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqIO import index as fasta_index, parse, write
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat import MatrixInfo as matlist
from copy import deepcopy
from easyFunctions import *
from GetOrfs import *
from HMMParser import HMM_Parser
from os import path, stat 
from pandas import DataFrame
from pickle import dump, load

####################### CRISPR READING TOOLS ####################### 
class FileWrapper:
    def __init__(self,filename):
        self.filename = filename
        self.index = -1
        self.lines = open(filename,"r").readlines()
    def readline(self):
        self.index += 1
        try: return self.lines[self.index]
        except:
            print ("Failed at index", self.index, "for file", self.filename)
            return "SUMMARY BY SIMILARITY--------"
    def has_more_lines(self): return len(self.lines)-1 != (self.index)

class MinCEDReader(dict):
    def __init__(self, filename): 
        self.ctype = False
        self.parse_file(FileWrapper(filename)); 
        self.calulate_consensus_seqs()
    
    def readLocus(self, line, newLocus,file):
        while "--------" not in line:
            repeat = ""
            try: pos, empty, repeat, spacer, lengths = line.strip().split("\t")
            except ValueError: pos, empty, repeat = line.strip().split("\t")
            newLocus.add_crRNA(pos, repeat, spacer)
            line = file.readline()
                
    def parse_file(self,file):
        line = ""
        newLocus = CRISPR("MinCED")
        while file.has_more_lines():
            if "Sequence " in line: newLocus.name, asmSize = self.getSeqName(line) #Line 1
            else: newLocus.name, asmSize = self.getSeqName(file.readline().strip()) #Line 1
            contID = newLocus.name
            if contID not in self: 
                self[contID] = CRISPR("MinCED")
                newLocus = self[contID]
                newLocus.name = contID
            else: newLocus = self[contID]
            file.readline() #Line 2 is always blank
            coordLine = file.readline().strip().split(" ") #Line 3
            newLocus.CRISPRarrayLocations.append(Coordinate(int(coordLine[5]), int(coordLine[7])))
            file.readline() #Line 4 is a header line: POSITION REPEAT SPACER
            file.readline() #Line 5 is a divider line w/a bunch of ---
            line = file.readline() #ResultLines
            self.readLocus(line,newLocus,file)
            file.readline() #Summary line
            file.readline() #Blank line
            line = file.readline().strip() #line with time to find repeats
            while "Time to" not in line:
                coordLine = line.split(" ")
                newLocus.CRISPRarrayLocations.append(Coordinate(int(coordLine[5]), int(coordLine[7])))
                file.readline()
                file.readline()
                line = file.readline()
                self.readLocus(line,newLocus,file)
                file.readline()#Summary line
                file.readline()#Blank line
                line = file.readline()#line with time to find repeats
            #2 blank lines
            file.readline()
            file.readline()
            self[newLocus.name] = newLocus
    
    def calulate_consensus_seqs(self): 
        for locus in self: self[locus].calc_consensus()
    
    def getSeqName(self,line):
        firstquote = line.find("'")+1
        name = line[firstquote:line.find("'",firstquote)]
        length = line[line.find("(")+1:line.find(" bp)")]
        return name,length

class PilerCRReader(dict):
    def __init__(self, filename): 
        self.ctype = True
        self.parse_file(FileWrapper(filename))
    def parse_file(self,fw):
        for i in range(11):fw.readline()
        while fw.has_more_lines():
            line = fw.readline()#Array # in file
            if "SUMMARY BY SIMILARITY" in line: break
            contID = fw.readline().strip().replace(">","") #Contig ID
            if " " in contID: contID=contID[:contID.find(" ")]
            if contID not in self: 
                self[contID] = CRISPR("PilerCR")
                newLocus = self[contID]
                newLocus.name = contID
            else: newLocus = self[contID]
            fw.readline() #Blank Line
            fw.readline() #Header
            fw.readline() #Divider line made of =======
            line = fw.readline() #1rst repeat
            while "==" not in line and "SUMMARY BY" not in line:
                crRNAInfo = line.strip().split("    ")
                # if len(newLocus.CRISPRarrayLocations) == 0: newLocus.CRISPRarrayLocations.append(Coordinate(int(crRNAInfo[0]), int(crRNAInfo[0])))
                # else: newLocus.CRISPRarrayLocations[0].end = int(crRNAInfo[0]) + len(crRNAInfo[-2])
                newLocus.add_crRNA(int(crRNAInfo[0]), crRNAInfo[-2], crRNAInfo[-1])
                line = fw.readline()
            if "SUMMARY BY" in line: break
            line = fw.readline().strip()#Summary Line
            newLocus.resolveRepeats(line)
            fw.readline() #Blank Line
            fw.readline() #Blank Line
class CRISPR:
    def __init__(self,crisprType):
        self.name = ""
        self.repeats = {}
        self.consensusRepeats = set()
        self.repeatCoords = []
        self.spacers = set()
        self.spacerCoords = []
        self.antiRepeats = {}
        self.antiRepeatCandidates = set()
        self.tracrRNACandidateSeqs = set()
        self.terminators = []
        self.numBlastRes = 0
        # self.possibleTracrs = set()
        self.crisprType = crisprType
        self.min,self.max = -1,-1
        self.hasNRegion = False
        self.CRISPRarrayLocations = []
    def add_crRNA(self,pos,repeat,spacer):
        self.spacers.add(spacer)
        if self.crisprType != 'PilerCR': 
            self.addrepeat(repeat) # pilerCR repeats don't need to be added
        self.repeatCoords.append(Coordinate(int(pos)-1,int(pos)+len(repeat)-1))
        self.spacerCoords.append(Coordinate(int(pos)+len(repeat),int(pos)+len(repeat)+len(spacer)-1))
    def addrepeat(self,repeat):
        try: self.repeats[len(self.CRISPRarrayLocations)][repeat]+=1
        except: 
            try: self.repeats[len(self.CRISPRarrayLocations)][repeat] = 1
            except: self.repeats[len(self.CRISPRarrayLocations)] = {repeat:1}
    def calc_consensus(self): #Only called for minCED arrays
        for id, seqDict in self.repeats.items(): #add 1 consensus per array id =array #
            maxCount, maxSeq = 0, ""
            for sequence, count in seqDict.items(): 
                if count>maxCount: maxSeq, maxCount = sequence, count
            self.consensusRepeats.add(maxSeq)
    def clusterBLASTResults(self, blastResults,protID):
        try:
            self.antiRepeats={}
            self.numBlastRes = len(blastResults)
            spacerLen = len(max(self.spacers, key=len)) + 10
            blastResults.sort()
            prevResult,nextResult = blastResults[0],blastResults[0]
            self.antiRepeats[prevResult] = AntiRepeatCandidate(prevResult,name=protID)
            if len(blastResults) == 1: return #If there is only 1 result stop
            self.repeats = blastResults
            for nextResult in blastResults[1:]:
                if prevResult.distance(nextResult) > spacerLen: 
                    self.antiRepeats[nextResult] = AntiRepeatCandidate(nextResult,name=protID)
                    if prevResult in self.antiRepeats: self.antiRepeats[prevResult].addDir()
                    else: self.antiRepeats[prevResult] = AntiRepeatCandidate(prevResult,'downstream',name=protID)
                prevResult=nextResult
            if nextResult in self.antiRepeats: self.antiRepeats[nextResult].addDir()
            else:  self.antiRepeats[nextResult] = AntiRepeatCandidate(nextResult,'downstream',name=protID)
        except: pass
    def getAntiRepeatCandidates(self,fh,chrSeq):
        for antiRepeat in self.antiRepeats.values(): fh.write(antiRepeat.getSeq(chrSeq)+"\n")
        fh.close()   
    def getTracrRNA_Candidates(self,erpOut,fh):
        crispr.terminators = []
        for i,terminator in enumerate(erpOut.terminators):
            seq = erpOut.records[terminator.name]
            tracrSeq = ""
            if terminator.upstream and not terminator.strand: tracrSeq = seq[terminator.Rholocation.start-1:].upper()
            elif not terminator.upstream and terminator.strand: tracrSeq = seq[:terminator.Rholocation.end].upper()
            if tracrSeq.count("N")>=4: continue
            if tracrSeq != "":
                terminator.seq = tracrSeq
                self.terminators.append(terminator)
                self.tracrRNACandidateSeqs.add(tracrSeq)
                write([SeqRecord(id="%s_%i" % (self.name,i),description='',seq=Seq(tracrSeq))],fh,'fasta')
                #fh.write(">%s_%i\n%s\n" % (,tracrSeq))
        return len(self.tracrRNACandidateSeqs)
    def repeatSeqs(self,protID,fh):
        self.name = protID
        for index, repeat in enumerate(self.consensusRepeats): fh.write(">%s_%i\n%s\n" % (self.name,index,repeat.upper()))
        fh.close()
    def resolveRepeats(self,summary_line):
        '''
        This method is to resolve repeats from pilercr which are given as matches
        being '.' and SNPs given as the actual base. @parameter summary_line looks like:
                29      46              29                GTTGTGAATAGCTTTCAAAATTGTATCTTAGTAGATGATTCACAGG
        NumRepeats RepeatLen   %id(blank)  SpacerLen    Consensus
        '''
        conRepeat = summary_line.split("    ")[-1].replace("-","")
        conRepeat = conRepeat
        self.consensusRepeats.add(conRepeat)
    def setName(self,name): 
        self.name = name
        # self.crispr
####################################################################
####################### Entire Operon Tools  ####################### 
def alignSequence(refSeq,fragment):
    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -.8
    probe = refSeq #scrubSequence(refSeq)
    ext_plus_lig = fragment
    alns = pairwise2.align.globalds(probe, ext_plus_lig, matrix, gap_open, gap_extend)
    top_aln = alns[0]
    aln_probe, aln_arms, score, begin, end = top_aln
    return [aln_probe, aln_arms]

def baseFile(fname): return fname[:fname.rfind(".")], fname[fname.rfind("."):]

def calcPercIdent(seqA,seqB):
    seq1,seq2 = alignSequence(seqA,seqB)
    matches = 0
    for index in range(len(seq1)):
        if seq1[index] == seq2[index]:matches += 1
    return float(matches)/len(seqA)*100

def duplexPercIdent(tracrRNA,crRNA):
    crLen = len(crRNA)
    if len(tracrRNA) >= crLen: return
    fmatches,rmatches = 0,0
    for i in range(crLen): 
        fmatches += int(tracrRNA[i]==crRNA[i])
        rmatches += int(tracrRNA[-(i+1)]==crRNA[-(i+1)])
    revComp = RC(crRNA)
    rfmatches,rrmatches = 0,0
    for i in range(crLen): 
        rfmatches += int(tracrRNA[i]==revComp[i])
        rrmatches += int(tracrRNA[-(i+1)]==revComp[-(i+1)])
    crLen = float(crLen)
    return max(fmatches/crLen, rmatches/crLen, rfmatches/crLen, rrmatches/crLen)*100

def findNCoords(seq):
    seq = seq.upper()
    index = -1
    coords = []
    while seq.find("N"*20,index+1) > -1:
        index = seq.find("N"*20,index+1)
        start_index = index
        while index<len(seq) and seq[index] == 'N': index+=1
        coords.append(Coordinate(start_index,index))
    return coords

def hmmHasResults(hmmFile, ext):
    fsize = stat(hmmFile).st_size
    return fsize <= 1900 or ext

def pullSeq(seqRec,coord):
    if not coord.strand:  return seqRec.seq[coord.start:coord.end].reverse_complement()
    return seqRec.seq[coord.start:coord.end]

class AntiRepeatCandidate:
    def __init__(self,coord,direction='upstream',name=""):
        self.name=name
        self.location = coord
        self.directions = set([direction])
    def addDir(self): self.directions.add('downstream')
    def getSeq(self,chrSeq):
        retSeqs = []
        buffer = 350
        seqName = ">" + self.name + "_%i_%i_%s"
        for dir in self.directions:
            if dir == 'upstream':
                start, end = max(self.location.start-buffer,0), self.location.end+1
                retSeqs.append(seqName % (start,end,'U'))
                retSeqs.append(chrSeq[start:end])
            else:
                start, end = self.location.start-1, min(self.location.end+buffer,len(chrSeq)-1)
                retSeqs.append(seqName % (start,end,'D'))
                retSeqs.append(chrSeq[start:end])
        return "\n".join(retSeqs)
    def __hash__(self): return hash(self.location)
    def __str__(self): 
        dirs = "upstream/downstream"
        if len(self.directions) == 1: dirs = self.directions.pop()
        return str(self.location) + " "+dirs
    def __lt__(self,other): return self.location < other.location
    def __gt__(self,other): return self.location > other.location

def casAnnotation(orfs,blastResultsFile,minCoord, maxCoord):
    fileHandle = FileWrapper(blastResultsFile)
    fileHandle.readline()
    
    while fileHandle.has_more_lines():
        line = fileHandle.readline().strip().split("\t") 
        bigDescription = line[-1]
        description = line[-1].lower()
        if "hypothetical" in description:continue
        words = description.split(" ")
        goodWords = set(["crispr","endonuclease","cas9","csn1","cas1","cas2","crrna","tracrrna","crispr-associated","csn2","dna2/cas4","cas4","cas1-like", "cas1/cas4","precrrna", "crrna", "tracr", "tracrrna", "tnpb", "pre-crrna"])
        goodWordFound = False
        for word in words:
            goodWordFound = (word in goodWords)
            if goodWordFound:break
        if not goodWordFound:continue
        found = False
        ID = line[0][:line[0].rfind("_")]
        ORF_ID = line[0]
        for feature in orfs.records[ID].features:
            if len(feature.qualifiers['label']) == 0:continue
            if feature.qualifiers['label'][0] == ORF_ID and not found:
                feature.qualifiers['label'] = bigDescription[bigDescription.find(" ")+1:bigDescription.find("[")]
                feature.qualifiers['note'] = bigDescription + " " + ORF_ID
                feature.qualifiers['blast_hit'] = "True"
                minCoord = min(feature.location.start,feature.location.end,minCoord)
                maxCoord = max(feature.location.start,feature.location.end,maxCoord)
                found = True
                break
    for id in orfs.records:
        remove = []
        for i in range(len(orfs.records[id].features)):
            feature = orfs.records[id].features[i]
            if 'blast_hit' not in feature.qualifiers:remove.append(i)
        for index in remove[::-1]: del orfs.records[id].features[index]
    return orfs, minCoord, maxCoord 


class CasOperon:
    def __init__(self,asmName,crisprFile):
        self.assembly = asmName
        self.chrWoperon = {}
        self.crispr = {}
        self.crisprFiles = set([crisprFile])
        self.prots = {}
        self.seq = None
        self.protID = ""
    def addCRISPR(self,crisprFile): self.crisprFiles.add(crisprFile)
    # def _setCrispr(self):
    #     if type(self.crispr) != type({}): return self.crispr
    #     try: self.crispr = self.crispr[self.protID]['pilerCR']
    #     except:self.crispr = self.crispr[self.protID]['minCED']
    def setName(self,protID):
        self.protID = protID
        if '%s' in self.chrAsmName: 
            self.chrAsmName  = self.chrAsmName % (protID)
            self.repeatFPath = self.repeatFPath % (protID)
    def getSeq(self): return str(self.seq.seq)
    def getCRISPR(self,protID):
        try:    return self.crispr[protID][True]
        except: return self.crispr[protID][False] #'minCED'
    def getFastaPointer(self,protID): return "assemblies/pseudoChromos/%s.fasta" % (protID)
    def getRepeatPath(self,protID): return "sequences/conRepeats/%s.fasta" % (protID)
    def annotate(self,protID,infernalResults):
        minCoord, maxCoord = len(self.seq), 0 # Set oppisite so it can shrink/grow respectively

        #BLAST results for possible protein matches in NR

        blastResultsFile = "/mnt/research/germs/shane/transActRNA/data/annotations/blastout/%s.blastout" % (protID)

        if not path.exists(blastResultsFile):
           # print ("!File exists", blastResultsFile)
           # die
            return

        orfs = GetORFs(self.getFastaPointer(protID),cutoff=120,writeOutPut=False)         #GetOrfs("proteins/orfs/%s.orfs" % (protID),120) # GetORFs(self.getFastaPointer(protID),120) # annotations/blastout/$f.blastout
        orfs, minCoord, maxCoord = casAnnotation(orfs,blastResultsFile, minCoord, maxCoord)
        
        crispr = self.getCRISPR(protID)
        #ADD the CRISPR ARRAY HITs/BLAST
        for repeatIndex, repeat in enumerate(crispr.repeats):
            fid = "Repeat_%i" % (repeatIndex)
            minCoord =  min(repeat.start,repeat.end,minCoord)
            maxCoord =  max(repeat.start,repeat.end,maxCoord)
            repeatFeature = SeqFeature(FeatureLocation(min(repeat.start,repeat.end), max(repeat.start,repeat.end)), type="repeat", strand=0)
            repeatFeature.qualifiers['label'] = [fid]
            repeatFeature.qualifiers['note'] = "color: #339966"
            orfs.records[protID].features.append(repeatFeature)
        
        for repeatIndex,repeat in enumerate(crispr.antiRepeats):
            fid = "Theoretical antirepeat_%i" % (repeatIndex)
            minCoord =  min(repeat.start,repeat.end,minCoord)
            maxCoord =  max(repeat.start,repeat.end,maxCoord)
            repeatFeature = SeqFeature(FeatureLocation(min(repeat.start,repeat.end), max(repeat.start,repeat.end)), type="CDS", strand=0)
            repeatFeature.qualifiers['label'] = [fid]
            repeatFeature.qualifiers['note'] = "color: #ff9900"
            orfs.records[protID].features.append(repeatFeature)

        #Get N coordinates
        NCoords = findNCoords(self.getSeq())
        for coord in NCoords:
            minCoord =  min(coord.start,coord.end,minCoord)
            maxCoord =  max(coord.start,coord.end,maxCoord)
            ambigSeq = SeqFeature(FeatureLocation(coord.start,coord.end), type="ambigious_seq", strand=0)
            ambigSeq.qualifiers['label'] = ["AmbigiousSequence"]
            orfs.records[protID].features.append(ambigSeq)
            self.hasNRegion = True

        self.min,self.max = minCoord, maxCoord
        
        ##Look through terminators
        crispr.terminators.sort(key=lambda x: x.genomeLocation.start)
        for repeatIndex, terminator in enumerate(crispr.terminators):
            strand = -1
            if terminator.strand: strand = 1
            repeatFeature = SeqFeature(FeatureLocation(terminator.genomeLocation.start, terminator.genomeLocation.end), type="terminator", strand=strand)
            fid = "RhoTerm_%i" % (repeatIndex)
            repeatFeature.qualifiers['label'] = [fid]
            orfs.records[protID].features.append(repeatFeature)

        for clusterID, tracr in infernalResults.items():
            strand = -1
            if tracr.strand: strand = 1
            tracrFeature = SeqFeature(FeatureLocation(tracr.location.start, tracr.location.end), type="misc_feature", strand=strand)
            tracrFeature.qualifiers['label'] = [clusterID]
            orfs.records[protID].features.append(tracrFeature)

    #     #shift annoations by buffer
    #     orfs.records[protID].seq = orfs.records[protID].seq[max(minCoord-500,0): min(maxCoord+500,len(operon.seq))]
    #     for rec in orfs.records[protID].features:
    #         rec.location = FeatureLocation(rec.location.start - max(self.min-500,0), rec.location.end - max(self.min-500,0))
        
        # Step 2k Create setup for GB file for Vector NTI
        for id in orfs.records:
            print ("annotations/genbank/%s.gb" % (protID))
            fh = open("annotations/genbank/%s.gb" % (protID),"w")
            biopythonID = id[:id.find(":")-1] #Biopython limits genbank file ids to 16 characters
            orfs.records[id].id = biopythonID[:16]
            orfs.records[id].name = biopythonID[:16]
            write([orfs.records[id]],fh,'genbank')
            fh.close()
        return protID

class CRISPRs(dict):
    def __init__(self): self.numCrisprFiles = 0
    def hasCrispr(self,crisprFiles,toolType,crisprPath,assemblyPath):
        #Variables to check crispr files
        nCrisprs = len(crisprFiles)
        validExts = set([".pcrout",".mnout"])
        percCutoff = int(nCrisprs*.15)
        counter = 0
        self.numCrisprFiles += nCrisprs
        if toolType: minSize = 200 #PilerCR File
        else: minSize = 0 #MinCED File

        #What needs tp be done?
        print("Working on checking %i CRISPRs from %s" % (nCrisprs,crisprPath))
        # needToCheck = set()
        # for fileName in crisprFiles:
        #     baseAsmName, ext = baseFile(fileName)
        #     if (baseAsmName not in self) and (ext in validExts): needToCheck.add(fileName)
        # crisprFiles = needToCheck
        
        #Checking CRISPRs
        for fileName in crisprFiles:
            counter += 1
            if counter % percCutoff == 0: print("\t%i%% of the way through with %i CRISPRs found" % (int((counter/float(nCrisprs))*100),len(self)))
            fsize = stat(crisprPath+"/"+fileName).st_size
            baseAsmName, ext = baseFile(fileName)
            if fsize <= minSize or ext not in validExts:continue
            try: self[baseAsmName].addCRISPR(crisprPath + fileName)
            except: self[baseAsmName] = CasOperon(assemblyPath[baseAsmName],crisprPath+fileName)

def readCRISPR(crisprFile):
    if ".pcrout" in crisprFile: return PilerCRReader(crisprFile)
    return MinCEDReader(crisprFile)

class CasOperons:
    def __init__(self,gene):
        self.operons = {}
        self.gene = gene
        self.multipleCas9s = {}
        self.casOnMultipleChrs = set()
        self.seqMap = seqDict()
    def items(self): return self._iter()
    def keys(self): return self.operons.keys()
    def __len__(self): return len(self.operons)
    def __setitem__(self, key, locus): self.operons[key] = locus
    def __iter__(self): return iter(self.operons.items())
    def __getitem__(self, key):
        return self.operons[self.seqMap[key]]

        # if key in self.operons: 
        #     return self.operons[self.seqMap[key]]
            #self.operons[key].setName(key)
            # return self.operons[key] #super(CasOperons, self).__getitem__(key)
        #self.operons[self.revMap[key]].setName(key)
        #return self.operons[self.revMap[key]] #super(CasOperons, self).__getitem__(self.revMap[key])       
    def hasCas9(self,hmmResultsDir,crisprs):
        hmmFiles = os.listdir(hmmResultsDir)
        nHmmFiles = len(hmmFiles)
        print("Processing %i hmm files" % (nHmmFiles))
        
        # Variables for checking hmm files
        validExts = set([".hmmout"])
        fivePercent = int(nHmmFiles*.05)
        counter = 0

        #Go through each file
        for fileName in hmmFiles:
            counter += 1
            if counter % fivePercent == 0: print("\t%i%% of the way through" % (int((counter/float(nHmmFiles))*100)))
            
            #Check that there are hmm hits
            baseAsmName, ext = baseFile(fileName)
            if hmmHasResults(hmmResultsDir+"/"+fileName ,ext not in validExts): continue
            
            ##Parse the hmm results
            proteins = HMM_Parser(hmmResultsDir+fileName).results
            
            ##Get the CHRs with proteins
            chrsWithCas = {}
            for protID,rec in proteins.items():
                baseProtID = protID[:protID.rfind("_")]
                try: chrsWithCas[baseProtID].add(protID) #remove the _ORF## extension of the protID to get the pseudo CHR
                except: chrsWithCas[baseProtID] = set([protID])
            
            ##If there are hmm hits then check that the hit has a CRISPR on the same CHR
            ##Check for crispr overlap and store results
            operon = crisprs[baseAsmName]
            for crisprFile in operon.crisprFiles:
                #Read the CRISPR file
                crisprResults = readCRISPR(crisprFile)

                #Get the chrs that have both a protein and a crispr
                chrOverlap = set(chrsWithCas.keys()).intersection(crisprResults.keys())
                if len(chrOverlap) == 0: continue

                #For each chromosome and each protein on the chromosomes
                for chromosome in chrOverlap:
                    for protID in chrsWithCas[chromosome]:

                        #add mapping information to all operons object 
                        self.seqMap[baseAsmName] = [chromosome,protID]

                        #Store the protein
                        operon.prots[protID] = proteins[protID] 

                        #Store the Crispr results
                        try: operon.crispr[protID][crisprResults.ctype] = crisprResults[chromosome]
                        except: operon.crispr[protID] = {crisprResults.ctype:crisprResults[chromosome]}

                        #Keep track of chr-prot relationship
                        try: operon.chrWoperon[chromosome].add(protID)
                        except: operon.chrWoperon[chromosome] = set([protID])

            # If there are proteins that also have a crispr, store them
            if len(operon.prots) > 0: self[baseAsmName] = operon
        self.save()
    
    def uniqueNukeSeqs(self,allCasAssembliesFile,allCasAminoAcidsFile):
        allCasSeqs = open(allCasAssembliesFile,"w")
        allAASeqs = open(allCasAminoAcidsFile,"w")
        uniqNukSeqMap, uniqNukSeqs, addedIDs, protSeqs = {}, {}, set(), set()
        counter = 0
        percCutoff = int(len(self)*.10)
        casOnMultipleChrs = set()
        metadata = open("tables/Cas9_Chr_Metadata.tsv","w")
        # for each assembly
        for asmName, operon in self:
            counter += 1
            if counter % percCutoff == 0: print("Made it through %i of the operons" % (counter))

            #If the assembly has multiple Cas9s on different chrs
            if len(operon.chrWoperon) > 1: 
                casOnMultipleChrs.add(asmName)
                continue

            #Get the sequence of the assembly
            seqs = fasta_index(operon.assembly,"fasta")

            #For every chr that has a hit
            for chrName, protIDs in operon.chrWoperon.items():

                #Hash the nucleotide seq
                chrSeqHash = hash(str(seqs[chrName].seq).strip().upper())

                #Map the seqID to the hash
                uniqNukSeqMap[chrName] = chrSeqHash

                #If the nucleotide sequence is unique
                if chrSeqHash not in uniqNukSeqs:
                    operon.seq = seqs[chrName]
                    operon.seq.seq = operon.seq.seq.upper()
                    if len(protIDs) > 1: self.multipleCas9s[chrName] = protIDs
                    else: operon.seq.id = list(protIDs)[0]
                    for protID in protIDs:
                        #Prep the nucleotide and protein records
                        nuclRec = operon.seq
                        nuclRec.id = protID
                        protIndex = 0
                        while nuclRec.id in addedIDs:
                            nuclRec.id = protID + "_DuplicateID_%i" % (protIndex)
                            protIndex+=1
                        protRec = operon.prots[protID]
                        protRec.id = nuclRec.id
                        metadata.write("\t".join([protID,nuclRec.description])+'\n')
                        nuclRec.description = ''
                        protRec = ''

                        #Write the sequences to files
                        write(protRec,allAASeqs,"fasta")
                        write(nuclRec,allCasSeqs,"fasta")
                        if not path.exists("assemblies/pseudoChromos/%s.fasta" % (nuclRec.id)):
                            with open("assemblies/pseudoChromos/%s.fasta" % (nuclRec.id),'w') as fh: write(nuclRec,fh,"fasta")
                        
                        #Store the nucleotide sequence in the set of known nucleotide seqs and AA seqs
                        try: uniqNukSeqs[chrSeqHash].add(chrName)
                        except: uniqNukSeqs[chrSeqHash] = set([chrName])
                        protSeqs.add(str(hash(protRec.seq)))
                else: uniqNukSeqs[chrSeqHash].add(chrName)

        self.casOnMultipleChrs = casOnMultipleChrs
        allCasSeqs.close()
        allAASeqs.close()
        metadata.close()
        print("There were %i unique nucleotide sequences with %i unique proteins" % (len(uniqNukSeqs),len(protSeqs)))
        print("There are %i assemblies with a Cas9 in multiple pseudochromosomes" % (len(casOnMultipleChrs)))
        print("There are %i pseudochromosomes the have more than 1 Cas9" % (len(self.multipleCas9s)))
        dump(uniqNukSeqs,  open("pickles/%s_uniqSeqMap.p"    % self.gene,"wb"))
        dump(uniqNukSeqMap,open("pickles/%s_uniqSeqRevMap.p" % self.gene,"wb"))
        self.save()
    
    def getRepSeqs(self,hasGoodDomains,repSeqsFile,allCasRepsFile):
        repSeqs = open(repSeqsFile,'w')
        for rec in parse(allCasRepsFile,"fasta"):
            if rec.id in hasGoodDomains: write(rec,repSeqs,"fasta")
        repSeqs.close()
        print("Saved %i sequences to %s" % (len(hasGoodDomains),repSeqsFile))

    def save(self):
        nCasOps = len(self)
        if nCasOps>0:
            print("Saving progress for %i Cas Operons" % nCasOps)
            dump(self,open("pickles/%s_Operons.p" % self.gene, "wb"))
        
    def loadProgress(self):
        self.casOperons = load(open("pickles/casOperons.p","rb"))
        self.revMap   = load(open("pickles/protAssemblyMap.p","rb"))
                      
    # def CRISPR_assemblies(self): return iter(self.crisprs.items())


#######################      Depricated      #######################
# class TRACR_Locus(dict):
#     def __init__(self, name): 
#         self.name = name
#         self.structuralHit = {} # TRACRStructureName:[Coordinates]
#         self.overlaps = False
#         self.topHit = None
#     def _setTopHit(self):
#         #Only 1 hit
#         if len(self.structuralHit)==1:
#             [(self.topHit, v)] = self.structuralHit.items()
#             return
#         #Check overlap
#         df = DataFrame(index=self.structuralHit,columns=self.structuralHit)
#         for hitName1,hit1 in self.structuralHit.items():
#             for hitName2,hit2 in self.structuralHit.items():
#                 if hitName1 == hitName2: continue
#                 if hit1.location.overlaps(hit2.location):
#                     dif = hit1.location - hit2.location
#                     if df[hitName2][hitName1] != dif: df[hitName1][hitName2] = dif
#         maxID1,maxID2 = df.max(axis=1).idxmax(),df.max(axis=0).idxmax()
#         #If no overlap
#         if maxID1 != maxID1:
#             print (df)
#             return #This case is interesting because overlap of blast hit with tracr but no overlap of tracrs
#         #else if overlap
#         if len(self.structuralHit[maxID1].location)>len(self.structuralHit[maxID2].location): self.topHit = maxID1
#         else: self.topHit = maxID2    

#     def structureHit(self,structureHit):
#         for structList in structureHit:
#             for structHit in structList:self.structuralHit[structHit.structName] = structHit
#     def addBLASTHit(self,blastHit):
#         for hitID,coord in blastHit: #hit = [StructName,coordinate]
#             if hitID in self.structuralHit:
#                 if self.structuralHit[hitID].location.overlaps(coord): 
#                     self.structuralHit[hitID].overlapsBLAST = True
#                     self.overlaps = True
#         self._setTopHit()
#     def addHit(self,hitType,structureName,coordinate):
#         try: self[hitType][structureName] = coordinate
#         except: self[hitType] = {structureName:coordinate}
#     def __hash__(self):return hash(self.name)

# class repeatStruct:
#     def __init__(self,id,seq): self.id = id; self.seq = seq
#     def __hash__(self): return hash(str(self.id)+self.seq)
#     def __eq__(self,other): return self.id == other.id and self.seq == other.seq

# class CRISPRLocus:
#     def __init__(self,type):
#         self.name = ""
#         self.assembly_size = 0
#         self.CRISPRarrayLocations = []
#         self.repeats = {}
#         self.consensusRepeats = set()
#         self.repeatCoords = []
#         self.spacers = set()
#         self.spacerCoords = []
#         self.antiRepeats = set()
#         self.possibleTracrs = set()
#         self.type = type
#         self.orfs = None
#         self.min,self.max = -1,-1
#         self.hasNRegion = False
    
#     def annotate(self,seq,assemblyFile):           
#         minCoord, maxCoord = len(seq), 0
        
#         #BLAST for possible protein matches in NR
#         orfs = GetORFs(assemblyFile,90)
#         blastResults = "/mnt/research/germs/shane/transActRNA/data/annotations/blast/%s_orfs.blastout" % (self.name)
#         from os.path import exists
#         if not exists(blastResults):
#             print ("!File exists", blastResults)
#             die
#             return
        
#         fileHandle = FileWrapper(blastResults)
#         fileHandle.readline()
        
#         while fileHandle.has_more_lines():
#             line = fileHandle.readline().strip().split("\t") 
#             bigDescription = line[-1]
#             description = line[-1].lower()
#             if "hypothetical" in description:continue
#             words = description.split(" ")
#             goodWords = set(["crispr","endonuclease","cas9","csn1","cas1","cas2","crrna","tracrrna","crispr-associated","csn2","dna2/cas4","cas4","cas1-like", "cas1/cas4","precrrna", "crrna", "tracr", "tracrrna", "tnpb", "pre-crrna"])
#             goodWordFound = False
#             for word in words:
#                 goodWordFound = (word in goodWords)
#                 if goodWordFound:break
#             if not goodWordFound:continue
#             found = False
#             ID = line[0][:line[0].rfind("_")]
#             ORF_ID = line[0]
#             for feature in orfs.records[ID].features:
#                 if len(feature.qualifiers['label']) == 0:continue
#                 if feature.qualifiers['label'][0] == ORF_ID and not found:
#                     feature.qualifiers['label'] = bigDescription[bigDescription.find(" ")+1:bigDescription.find("[")]
#                     feature.qualifiers['note'] = bigDescription + " " + ORF_ID
#                     feature.qualifiers['blast_hit'] = "True"
#                     minCoord = min(feature.location.start,feature.location.end,minCoord)
#                     maxCoord = max(feature.location.start,feature.location.end,maxCoord)
#                     found = True
#                     break
#         for id in orfs.records:
#             remove = []
#             for i in range(len(orfs.records[id].features)):
#                 feature = orfs.records[id].features[i]
#                 if 'blast_hit' not in feature.qualifiers:remove.append(i)
#             for index in remove[::-1]: del orfs.records[id].features[index]
        
#         #ADD the CRISPR ARRAY HITs/BLAST
#         repeatIndex = 0
#         for repeat in self.repeatCoords:
#             fid = "Repeat_%i" % (repeatIndex)
#             minCoord =  min(repeat.start,repeat.end,minCoord)
#             maxCoord =  max(repeat.start,repeat.end,maxCoord)
#             repeatFeature = SeqFeature(FeatureLocation(min(repeat.start,repeat.end), max(repeat.start,repeat.end)), type="repeat", strand=0)
#             repeatFeature.qualifiers['label'] = [fid]
#             repeatFeature.qualifiers['note'] = "color: #339966"
#             orfs.records[self.name].features.append(repeatFeature)
#             repeatIndex+=1
#         repeatIndex = 0
#         for repeat in self.possibleTracrs:
#             fid = "Theoretical antirepeat_%i" % (repeatIndex)
#             minCoord =  min(repeat.start,repeat.end,minCoord)
#             maxCoord =  max(repeat.start,repeat.end,maxCoord)
#             repeatFeature = SeqFeature(FeatureLocation(min(repeat.start,repeat.end), max(repeat.start,repeat.end)), type="terminator", strand=0)
#             repeatFeature.qualifiers['label'] = [fid]
#             repeatFeature.qualifiers['note'] = "color: #ff9900"
#             orfs.records[self.name].features.append(repeatFeature)
#             repeatIndex+=1
        
#         #Get N coordinates
#         NCoords = findNCoords(seq)
#         for coord in NCoords:
#             if ((minCoord-500) <= coord.start and coord.start <= (maxCoord+500)) or (max(minCoord-500,0) <= coord.end and coord.end <= (maxCoord+500)):
#                 minCoord =  min(coord.start,coord.end,minCoord)
#                 maxCoord =  max(coord.start,coord.end,maxCoord)
#                 ambigSeq = SeqFeature(FeatureLocation(coord.start,coord.end), type="ambigious_seq", strand=0)
#                 ambigSeq.qualifiers['label'] = ["AmbigiousSequence"]
#                 orfs.records[self.name].features.append(ambigSeq)
#                 self.hasNRegion = True
        
#         self.min,self.max, self.orfs = minCoord, maxCoord,orfs
        
#     def confirmedAnti(self,anti_repeat):
#         remove = set()
#         foundAnti = False
#         for rep in self.antiRepeats:
#             foundAnti = anti_repeat.overlaps(rep)
#             if not foundAnti:remove.add(rep)
#             if foundAnti: break
#         if foundAnti: return foundAnti
#         # for coord in remove: self.antiRepeats.remove(coord)
#         return False

#     def combineResults(self,other):
#         if self.type != "Piler": self.consensusRepeats = other.consensusRepeats[0]
#         self.type = "Combined"
#         #Compare CRISPR Array Locations
#         locA_index = 0
#         for locationA in other.CRISPRarrayLocations:
#             foundLocation = False
#             locB_index = 0
#             for locationB in self.CRISPRarrayLocations:
#                 foundLocation = (locationA == locationB)
#                 if foundLocation:
#                     self.CRISPRarrayLocations[locA_index].start = min(self.CRISPRarrayLocations[locA_index].start,self.CRISPRarrayLocations[locB_index].start)
#                     self.CRISPRarrayLocations[locA_index].end = max(self.CRISPRarrayLocations[locA_index].end,self.CRISPRarrayLocations[locB_index].end)
#                     break
#                 locB_index += 1
#             locA_index += 1
#             if not foundLocation: self.CRISPRarrayLocations.append(locationA)
        
#         #Combine anti-repeats
#         for anti_repeatA in other.antiRepeats:
#             #If the anti-repeat is in the crispr array of the other found array then continue on to the next repeat
#             inCRISPRArray = False
#             for location in self.CRISPRarrayLocations:
#                 inCRISPRArray = (location.contains(anti_repeatA) or anti_repeatA.overlaps(location))
#                 if inCRISPRArray:break
#             if inCRISPRArray:continue
            
#             #Check to see if the anti-repeat is a duplicate of an existing anti    
#             foundAntiRepeat = False
#             for anti_repeatB in self.antiRepeats:
#                 foundAntiRepeat = (anti_repeatA == anti_repeatB)
#                 if foundAntiRepeat: break
#             if not foundAntiRepeat: self.antiRepeats.add(anti_repeatA)
        
#         #Combine Repeats
#         lenRepCoords = len(self.repeatCoords)
#         if lenRepCoords != len(set(self.repeatCoords)):
#             print (lenRepCoords,len(set(self.repeatCoords)))
#             raise Exception("How can a repeat coord have a start position that is the same?")
#         theRepeats = set(self.repeatCoords)
#         numFound = 0
#         for coordA in other.repeatCoords:
#             if coordA in theRepeats:pass
#             else:
#                 foundCoord = False
#                 for coordB in theRepeats:
#                     foundCoord = (coordA == coordB)
#                     if foundCoord:break
#                 if not foundCoord: self.repeats["NewRepeat_%i" %(numFound)] = coordA
    
#     def checkExtendsArray(self,coord,crRNALen):
#         print ("\t\t" +str(coord),end=" ")
#         if self.CRISPRarrayLocations[0].distance(coord) <= max(crRNALen,100):
#             self.extendArray(coord)
#             print ("Extends") #"\tExtends = TRUE for ",coord
#         else:
#             print ("tracrRNA Candidate")
#             return coord
    
#     def checkArrayBoundaries(self, blast_results):
#         print (self.name + " CRISPR Array:",self.CRISPRarrayLocations[0])
#         before,after = [],[]
#         for hit in blast_results:
#             if hit.start < self.CRISPRarrayLocations[0].start: before.append(hit)
#             else: after.append(hit) 

#         retHits = []
#         avgSpacerLen = 0
#         for spacer in self.spacers: avgSpacerLen += len(spacer)
#         avgSpacerLen = int(avgSpacerLen/float(len(self.spacers))) 
#         crRNALen = avgSpacerLen + len(list(self.consensusRepeats)[0])
#         print ("\tAverage spacer len: %i + Repeat len: %i = %i basepairs" % (avgSpacerLen,len(max(self.consensusRepeats,key=len)),crRNALen))
#         #For the each of the blast hits
#         ##Ignore any hits that are completely contained within the crispr array
#         ##If there is overlap but not complete containment, update the array
#         ##If the coordinate doesn't overlap, check to see if it extends
#         ###If it extends update the array and the repeats
#         ###Else add to list of hits to return
#         print ("\tBefore:")
#         for coord in before[::-1]: #For the each of the spacially ORDERED blast hits
#             if self.CRISPRarrayLocations[0].contains(coord): print ("\t\t",coord,"Within"); continue #Ignore any hits that are completely contained within the crispr array
#             elif self.CRISPRarrayLocations[0].overlaps(coord): 
#                 print ("\t\t",coord,"Overlap Extension")
#                 coord = self.extendArray(coord) #Check for a blast hit outside of the crispr array 
#             else: coord = self.checkExtendsArray(coord,crRNALen) #Check that the coordinate is completly contained within the CRISPR
#             if coord is not None: retHits.append(coord)
        
#         print ("\tAfter:")
#         for coord in after: 
#             if self.CRISPRarrayLocations[0].contains(coord): print ("\t\t",coord,"Within"); continue ##Ignore any hits that are completely contained within the crispr array
#             elif self.CRISPRarrayLocations[0].overlaps(coord): 
#                 print ("\t\t",coord,"Overlap Extension")
#                 coord = self.extendArray(coord) #Check for a blast hit outside of the crispr array 
#             else: coord = self.checkExtendsArray(coord,crRNALen) #Check that the coordinate is completly contained within the CRISPR
#             if coord is not None: retHits.append(coord)
       
#         print ("\tUpdated CRISPR Array:",self.CRISPRarrayLocations[0])
#         print ("\tPossible TRACRS:")
#         print ("\t\t",self.repeatCoords[0])
#         print ("\t\t",self.repeatCoords[-1])

#         for hit in retHits: print ("\t\t",hit)
#         print()
#         self.possibleTracrs.add(self.repeatCoords[0])
#         self.possibleTracrs.add(self.repeatCoords[-1])
#         for coord in retHits:self.possibleTracrs.add(coord)
#         return retHits, crRNALen    
        
#     def getTerminusRepeats(self,seq,index):
#         index+=1
#         retStr = ""
#         fastaSeq = ">Seq_%s_%s_%i_%i_%i\n%s\n" 
#         startCoord = self.repeatCoords[0]
#         endCoord = self.repeatCoords[-1]
#         retStr += fastaSeq % ("E",self.consensusRepeats[0],min(startCoord.start,startCoord.end),max(startCoord.end,startCoord.start),index,seq[min(startCoord.start,startCoord.end)-500:max(startCoord.end,startCoord.start)])
#         index+=1
#         retStr += fastaSeq % ("S",self.consensusRepeats[0],min(startCoord.start,startCoord.end),max(startCoord.end,startCoord.start),index,seq[min(endCoord.start,endCoord.end):max(endCoord.end,endCoord.start)+500])
#         return retStr,index
            
#     def extendArray(self,newCoord):
#         if newCoord.start < self.CRISPRarrayLocations[0].start or newCoord.end < self.CRISPRarrayLocations[0].start:
#             self.CRISPRarrayLocations[0].start = min(newCoord.start,newCoord.end)
#             self.repeatCoords = [newCoord] + self.repeatCoords
#         else: 
#             self.CRISPRarrayLocations[0].end = max(newCoord.start,newCoord.end)
#             self.repeatCoords.append(newCoord)

#     def addrepeat(self,repeat):
#         try: self.repeats[len(self.CRISPRarrayLocations)][repeat]+=1
#         except: 
#             try: self.repeats[len(self.CRISPRarrayLocations)][repeat] = 1
#             except: self.repeats[len(self.CRISPRarrayLocations)] = {repeat:1}
        
#     def add_crRNA(self,pos,repeat,spacer):
#         self.spacers.add(spacer)
#         self.addrepeat(repeat)
#         self.repeatCoords.append(Coordinate(int(pos)-1,int(pos) + len(repeat)-1))
#         self.spacerCoords.append(Coordinate(int(pos) + len(repeat),int(pos) + len(repeat) + len(spacer)-1))
    
#     def combineRepeats(self,repeatA,repeatB):
#         repeatAID,repeatASeq,repeatBID,repeatBSeq,countA,countB = "","","","",0,0
#         for ID,seqDict in self.repeats.items():
#             for repeat in seqDict:
#                 if repeat == repeatA: 
#                     repeatAID = ID
#                     repeatASeq = repeat
#                     countA = self.repeats[ID][repeat]
#                     break
#                 elif repeat == repeatB: 
#                     repeatBID = ID
#                     repeatBSeq = repeat
#                     countB = self.repeats[ID][repeat]
#                     break
#         if repeatASeq != "" and repeatBSeq != "":
#             if countA > countB:
#                 self.repeats[repeatAID][repeatASeq] += countB
#                 del self.repeats[repeatBID][repeatBSeq]
#                 if len(self.repeats[repeatBID]) == 0: del self.repeats[repeatBID]
#             else:
#                 self.repeats[repeatBID][repeatBSeq] += countA
#                 del self.repeats[repeatAID][repeatASeq]
#                 if len(self.repeats[repeatAID]) == 0: del self.repeats[repeatAID]
    
#     def removeDuplicates(self):
#         remove = set()
#         for IDA,seqDictA in self.repeats.items():
#             for IDB,seqDictB in self.repeats.items():
#                 if IDA != IDB:
#                     for seqA in seqDictA:
#                         for seqB in seqDictB:
#                             if seqA == seqB: 
#                                 if self.repeats[IDA][seqA] > self.repeats[IDB][seqB]: remove.add(repeatStruct(IDB,seqB))
#                                 else: remove.add(repeatStruct(IDA,seqA))
#         if len(remove) != 0:
#             for coord in remove:
#                 del self.repeats[coord.id][coord.seq]
#                 if len(self.repeats[coord.id]) == 0: del self.repeats[coord.id]
        
#     def calc_consensus(self): 
#         maxCount, maxSeq = 0, ""
#         # print(self.repeats)
#         for id, seqDict in self.repeats.items():
#             for sequence, count in seqDict.items(): 
#                 if count>maxCount: 
#                     maxSeq = sequence
#                     maxCount = count
#                # print sequence, count
#                # allRepeats.add(sequence)
#         self.consensusRepeats=[maxSeq]
#        # self.consensusRepeats=list(allRepeats)
#        # #Strategy 1 (Depricated): 
#        # '''
#        # If there is is only one crispr array location, then compute a consensus
#        # by aligning all of the repeats and taking the IUPAC for the allelels, ignore
#        # gaps if  gap_percentage >= .75 or gap_percentage <= .25 (typical cases meaning really rare or almost completly gapped
#        # '''
#        # if self.type == "Piler": raise Exception("This should not happen")
#        # if len(self.CRISPRarrayLocations) == 1:
#        #     alignments = alignSeqs(self.repeats)
#        #     seqs=[]
#        #     for alnmnt in alignments: seqs.append(str(alnmnt.seq))
#        #     self.consensusRepeats.append(calc_consensus_seq(seqs,True))
#        # else:
#        #     for id,seqs in self.repeats.items():
#        #         alignments = alignSeqs(seqs)
#        #         seqs=[]
#        #         for alnmnt in alignments: seqs.append(str(alnmnt.seq))
#        #         self.consensusRepeats.append(calc_consensus_seq(seqs,True))
               
#        #     AllPossibleRepeats = set()
#        #     for id, repeatDict in self.repeats.items():
#        #         for repeat in repeatDict: AllPossibleRepeats.add(repeat)
#        #     ignore = set()
#        #     for repeatA in AllPossibleRepeats:
#        #         for repeatB in AllPossibleRepeats:
#        #             if repeatA != repeatB and ("%s_%s" % (repeatA,repeatB) not in ignore):
#        #                 percIdent = CalcPercIdent(repeatA,repeatB)
#        #                 if percIdent > 90.0: self.combineRepeats(repeatA,repeatB)
#        #                 ignore.add("%s_%s" % (repeatA,repeatB))
#        #                 ignore.add("%s_%s" % (repeatB,repeatA))
#        #     self.removeDuplicates()
#        #     if len(self.CRISPRarrayLocations) == 1 and len(self.repeats[1]) == 1: #expect to blow up when dupes in 1-2 and 1 is deleted cause more in 2
#        #         for id,seqDict in self.repeats.items():
#        #             self.consensusRepeats = seqDict.keys()
#        #             if len(self.consensusRepeats) != 1: raise Exception("Problem here " +str(self.repeats))
#        #             return
#        #     alignments = alignSeqs(self.repeats)
#        #     for alnmnt in alignments: seqs.append(str(alnmnt.seq))
#        #     for id,seqDict in self.repeats.items():
#        #         for seq in seqDict:
#        #             self.consensusRepeats.append(seq)
#        #     self.consensusRepeats = calc_consensus_seq(seqs)
            
#     def resolveRepeats(self,summary_line):
#         '''
#         This method is to resolve repeats from pilercr which are given as matches
#         being '.' and SNPs given as the actual base. @parameter summary_line looks like:
#                 29      46              29                GTTGTGAATAGCTTTCAAAATTGTATCTTAGTAGATGATTCACAGG
#         NumRepeats RepeatLen   %id(blank)  SpacerLen    Consensus
#         '''
#         conRepeat = summary_line.split("    ")[-1].replace("-","")
#         conRepeat = conRepeat
#         self.consensusRepeats.add(conRepeat)
#         #TODO: Could resolve repeats into sequence but for now will leave them as ....A....T......

#     def findAntiRepeats(self,record):
#         if len(self.consensusRepeats) == 0 and len(self.repeatCoords)!= 0: self.calc_consensus()
#         if len(self.consensusRepeats) == 0: raise Exception("No repeats for", self.name, len(self.repeatCoords))
#         write([record],"tmp/assemblyFasta.fa","fasta")
#         repeatFasta = MakeFasta(tempfile.NamedTemporaryFile(delete=False))
#         counter = 0
#         for repeat in self.consensusRepeats:
#             repeatFasta.write("%s_%i" %(self.name,counter), repeat)
#             counter += 1
#         repeatFasta.close()
#         results = BLAST(repeatFasta.name,"tmp/assemblyFasta.fa","tmp/blast_results.out")
#         repeatFasta.die()
#         remainingResults = self.compareResults(results) # This assumes that a hit is the same if it is contained or within the CRISPR
#         for result in remainingResults: self.antiRepeats.add(result)
#             #if Orientation(result): self.antiRepeats.append(result)
#             #else: self.repeatCoords.append(result)

#     def compareResults(self,results):
#         #Assumption: The anti repeat is not in the crispr array and doesn't overlap it
#         retList = []
#         for result in results:
#             foundInArray = False
#             overlapsLocation = False
#             for location in self.CRISPRarrayLocations:
#                 foundInArray = location.contains(result)
#                 if foundInArray: break
#                 overlapsLocation = location.overlaps(result)
#                 if overlapsLocation: break
#             minSeparationFromArray = min(abs(location.start-result.start),abs(location.start-result.end),abs(location.end-result.end),abs(location.end-result.start))
#             if not foundInArray and minSeparationFromArray > 1 and not overlapsLocation: retList.append(result)
#         return retList

#     def consensus(self):return self.calc_consensus()
    
#     def __str__(self): 
#        # crRNAs = ""
#        # for index in range(0,len(self.repeatCoords)): crRNAs += "%i\t%i\t%i\t%i\n" %(self.repeatCoords[index].start, self.repeatCoords[index].end, self.spacerCoords[index].start,self.spacerCoords[index].end)
#        #  AssemblyName \t CRISPR_Start \t CRISPR_End \t Consensus Repeat \n
#        # print "%s\t%i\t%i\t%s\n" % (self.name,self.CRISPRarrayLocations.start,self.CRISPRarrayLocations.end,self.consensusRepeats)
#         retstr = ""
#         for crArray in self.CRISPRarrayLocations: retstr += "%s\t%i\t%i\t%s\tcrRNAs:%i\n" % (self.name,crArray.start,crArray.end,str(self.consensusRepeats),len(self.spacers))         
#         return retstr
'''
        print "Working on",self.name
        seqFile = open("tmp/%s_genomicSequence.fasta" % (self.name),"w")
        seqFile.write(">%s\n%s\n" % (self.name,seq))
        seqFile.close()
        newFasta = open("tmp/%s_orfs.fasta" %(self.name),"w")
        orfs = getOrfs("tmp/%s_genomicSequence.fasta" % (self.name),90)
        for seq in orfs.Seqlist: writeFasta(seq, newFasta, "fasta")
        newFasta.close()
        allBlast = open("AnnotateCMDs.sh","a")
        cmd = "bsub -P gis-jbrowse-blastall -J CasBLAST -q prod -R \"rusage[mem=1048,scr=1000]\" " \
        "-R \"span[hosts=1]\" -n 12 -e tmp/%s_blast.err -o tmp/%s_blast.out \"bash tmp/%s.sh; rm -f tmp/%s.sh tmp/%s_genomicSequence.fasta tmp/%s_orfs.fasta\""% (self.name,self.name,self.name,self.name,self.name,self.name)
        blastCMDs = "blastp -db blast/db/nr "\
        "-query tmp/%s_orfs.fasta -outfmt \"7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames stitle\" " \
        "-num_threads 12 -out data/Annotations/%s_orfs.blastout\n" % (self.name,self.name)
        cmdFile = open("tmp/%s.sh"%(self.name),"w")
        cmdFile.write(blastCMDs)
        cmdFile.close()
        allBlast.write(cmd+"\n")
        allBlast.close()
'''

'''
write(sequences, handle, format)
    Write complete set of sequences to a file.

    Arguments:
     - sequences - A list (or iterator) of SeqRecord objects, or (if using
       Biopython 1.54 or later) a single SeqRecord.
     - handle    - File handle object to write to, or filename as string
       (note older versions of Biopython only took a handle).
     - format    - lower case string describing the file format to write.

    Note if providing a file handle, your code should close the handle
    after calling this function (to ensure the data gets flushed to disk).

    Returns the number of records written (as an integer).
'''