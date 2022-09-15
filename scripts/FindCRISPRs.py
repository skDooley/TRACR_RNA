from datetime import date
import glob
import sys
from os import access, chdir, listdir, mkdir, path, system, W_OK
import pickle

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--fastaDict",  dest="FASTA_DICT", help="The pickled (.p) file containing the dictionary with AssemblyID:path/to/fasta_file.")
parser.add_option("-c", "--CRISPRDir",  dest="CRISPRDir", help="The output location of the CRISPRs")  
parser.add_option("-e", "--envPath",  dest="env", help="The name of the env that has pilercr and minced.E.g.\n\t~/anaconda3/bin/activate\n\t\"~/anaconda3/bin/activate cas\"")  
parser.add_option("-q", "--queue",  dest="q", help="The name slurm queue")  

(opts, args) = parser.parse_args()

try:
    if not path.exists(opts.CRISPRDir): mkdir(opts.CRISPRDir)
    if not access(opts.CRISPRDir,W_OK): parser.error("The directory \'%s\' in not writeable." % (opts.SCRIPTS_DIR)); sys.exit(-2)
except: parser.error("Unable to create a directory at \'%s\'" % (opts.CRISPRDir))

#Paths
dba_base = opts.CRISPRDir
refDatabases = pickle.load(open(opts.FASTA_DICT,'rb'))
crisprSCRIPT_Path = path.realpath("./")

#Get Todays date for the scripts
theDate = str(date.today())

#HPC Header
header = ""
for line in open(crisprSCRIPT_Path + "/scripts/hpc/SlurmHeader.qs"): header+=line
header.replace("bioinformatics.q", opts.q)

firstFileName = crisprSCRIPT_Path+"/scripts/hpc/CRISPR_run_%i_%s.sb" % (0,theDate)
bashLog = open(firstFileName, 'w')
bashLog.write(header)
bashLog.write("source " + opts.env + '\n')
logNumber,count =0,0
validExts = set([".fna",".fasta",".fa"])
asmCount = 0

for asmID, assembly in refDatabases.items():
    asmID = asmID[:asmID.rfind(".")]
    ext = assembly[assembly.rfind("."):]

    if ext not in validExts: continue
    
    #Add PilerCR command
    pcrPath = "%s/%s.pcrout" % (dba_base, asmID)
    if not path.exists(pcrPath):
        asmCount += 1
        pcmd = "pilercr -minid 0.85 -mincons 0.8 -noinfo -in %s -out %s 2>/dev/null" % (assembly, pcrPath)
        count+=1
        bashLog.write(pcmd+"\n")

    #Add minCED command
    mincedPath = pcrPath.replace("pcrout","mnout")
    if not path.exists(mincedPath): 
        asmCount += 1
        mcmd = "minced -minRL 16 -maxRL 64 -maxSL 64 -searchWL 6 %s %s" % (assembly, mincedPath)
        count+=1
        bashLog.write(mcmd+"\n")
print("Prepared to run minced and pilercr on a total of %i searches" % (asmCount))
print("To run, navigate to the TRACR_RNA folder and submit the job to the SLURM manager:\n\tqsub "+ firstFileName) 
print("Note, if you don't want to run the search on a slurm hpc, just run the script directly from the TRACR_RNA folder with\n\tbash "+ firstFileName[firstFileName.find("TRACR_RNA")+10:])
bashLog.close() 