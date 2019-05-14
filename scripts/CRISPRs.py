import datetime
import glob
import sys
#sys.path.append("scripts/")
#from CRISPRtools import PilerCRReader, MinCEDReader
from os import chdir, listdir, path, system

count = 0 
theDate = str(datetime.date.today())
# errorLogName = "CRISPR_ErrorLog_%s.log" % theDate
# while (path.exists(errorLogName)): errorLogName = "CRISPR_ErrorLog_%s_%s.log" % (theDate,count); count += 1
# error_log = open(errorLogName,"w")

outdir = "/mnt/research/germs/shane/databases/crisprs"
bashLog = open(outdir+"/scripts/CRISPR_run_%i_%s.sh" % (2,theDate), 'w')

assembly_dir = "/mnt/research/germs/shane/databases/assemblies/NCBI/refseq/bacteria"
chdir(assembly_dir)
assemblies = glob.glob("*.fasta")





logNumber = 3
print ("Number of assemblie to run CRISPR array detection on:",len(assemblies))
for assembly in assemblies:

    #Run PilerCR
    if path.exists("%s/%s.pcrout" % (outdir, assembly.replace(".fasta",""))):retCode1 = 0
    else: 
        pcmd = "pilercr -minid 0.85 -mincons 0.8 -noinfo -in %s -out %s/%s.pcrout 2>/dev/null" %(assembly, outdir, assembly.replace(".fasta",""))
        count+=1
        # print cmd
        bashLog.write(pcmd+"\n")
    if (count % 500 == 0):
        bashLog.close()
        bashLog = open(outdir+"/scripts/CRISPR_run_%i_%s.sh" % (logNumber,theDate), 'w')
        logNumber+=1
bashLog.close() 
        #retCode1 = system(pcmd)
        
    # #Run minced
    # if path.exists("%s/%s.mnout" % (outdir, assembly)):retCode2 = 0
    # else: 
    #     mcmd = "minced -minRL 16 -maxRL 64 -minSL 8 -maxSL 64 -searchWL 6 %s %s/%s.mnout" % (assembly_dir+assembly, outdir, assembly)
    #     # print cmd
    #     retCode2 = system(mcmd)
    
    # if retCode1 != 0: 
    #     error_log.write(str(retCode1) + "  " + pcmd+"\n")
    #     break
    # if retCode2 != 0: 
    #     print "%i\nminced -maxSL 75 --maxRL 75 -minRL 16 -minSL 20 -searchWL 6 %s %s/%s.mnout" % (retCode2,assembly_dir+assembly, outdir, assembly)
    #     break
