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
        start, end = info[2].split("..")
        self.Rholocation = Coordinate(start,end)
        self.fwd = int(self.name[-1]) % 2 != 0
    def __str__(self): return "%s\t%s\t%s\t%s\t" % (self.name,str(self.strand),str(self.Rholocation),str(self.fwd))