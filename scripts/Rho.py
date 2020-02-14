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
            for i in range(9): file.readline()
            capture = False
            for line in file:
                if capture: 
                    line = line.strip().replace("  "," ").replace("  "," ").replace("  "," ")
                    self.terminators.append(RhoIndTerminator(seqName,line.split(" ")))
                capture = (">" in line)
                if capture: seqName = line.strip().replace(">","")

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