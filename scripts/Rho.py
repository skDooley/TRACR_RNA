from Bio.SeqIO import parse
from easyFunctions import Coordinate
# from xlrd import open_workbook as read_excel
# from xlrd import empty_cell
from pyexcel import get_array as read_excel


class ErpinOut:
    def __init__(self, outfile="tmp/rhoInd.out", inputfile="tmp/possibleTracrs.fasta"):
        self.terminators = []
        self.records={}
        for rec in parse(inputfile,"fasta"): self.records[rec.id]=str(rec.seq)
        self.numRecords = len(self.records)
        with open(outfile) as file:
            for i in range(10): file.readline()
            skip=False
            for line in file:
                if line.strip() == '':break
                if skip:
                    skip = not skip
                    continue
                if (">" in line): seqName = line.strip().replace(">","")
                else: 
                    skip = not skip
                    line = line.strip().replace("  "," ").replace("  "," ").replace("  "," ")
                    newRho = RhoIndTerminator(seqName,line.split(" "))
                    if (newRho.strand and newRho.upstream) or (not newRho.strand and not newRho.upstream): continue
                    try:
                        lastRho = self.terminators[-1]
                        if lastRho.name == newRho.name and lastRho.strand == newRho.strand: #same location and strand
                            if not newRho.strand: self.terminators[-1]=newRho #Replace because this on is shorter
                        else: self.terminators.append(newRho)
                    except: self.terminators.append(newRho)
                    # try:
                    #     lastRho = self.terminators[-1]
                    #     if lastRho.upstream == newRho.upstream and lastRho.name == newRho.name and lastRho.strand == newRho.strand:
                    #         lastShorter = lastRho.Rholocation.end < newRho.Rholocation.end
                    #         if newRho.upstream and lastShorter: self.terminators[-1] = newRho
                    #         elif not newRho.upstream and not lastShorter: self.terminators[-1] = newRho
                    #     else: self.terminators.append(newRho)
                    # except: self.terminators.append(newRho)

class RhoTermPredictOut:
    def __init__(self,seqFile="tmp/possibleTracrs.fasta",rhoPredictFile="tmp/RhoPredictions.xlsx"):
        self.terminators = []
        self.records={}
        for rec in parse(seqFile,"fasta"): self.records[rec.id]=str(rec.seq)
        self.numRecords = len(self.records)
        predictions = read_excel(file_name=rhoPredictFile)
        for rec in predictions[1:]: #Start in 2nd element because first is header
            self.terminators.append(RhoIndTerminator(rec[0],rec,'rhoterm'))
                   
class RhoIndTerminator:
    def __init__(self,name,info,rtype = "erpin"):
        self.name = name
        self.upstream = self.name[-1] == 'U'
        self.seq = None
        locRec = name.split("_")
        self.genomeLocation = Coordinate(locRec[-2],locRec[-3]) #Location w/in the genome

        if rtype == "erpin":
            self.strand = (info[0]=="FW")
            start, end = info[2].split("..")
            self.Rholocation = Coordinate(start,end) #Location within the subsequence
        else:
            self.strand = (info[3]=="plus")
            start, end = info[1], info[2]
            self.Rholocation = Coordinate(start,end)
    def __str__(self): return "%s\t%s\t%s\t%s\t" % (self.name,str(self.strand),str(self.Rholocation),str(self.upstream))