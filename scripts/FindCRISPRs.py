from datetime import date
import glob
import sys
from os import access, chdir, listdir, mkdir, path, system, W_OK
import pickle

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--fastaDict",  dest="FASTA_DICT", help="The pickled (.p) file containing the dictionary with AssemblyID:path/to/fasta_file.") 
parser.add_option("-s", "--scriptsDir",  dest="SCRIPTS_DIR", help="The directory containing where the bash script(s) will be built to find the CRISPR arrays. If the directory does not exist, it will automatically be created.") 
parser.add_option("-r", "--removeOld",  dest="REMOVE_OLD", action="store_true", default=False, help="Removes all '*.sh' files from the '--scriptsDir'. Default = False") 
(opts, args) = parser.parse_args()

try:
    if not opts.SCRIPTS_DIR: parser.error("The directory where the scripts will be written is required."); sys.exit(-1)
    if not path.exists(opts.SCRIPTS_DIR): mkdir(opts.SCRIPTS_DIR)
    if not access(opts.SCRIPTS_DIR,W_OK): parser.error("The directory \'%s\' in not writeable." % (opts.SCRIPTS_DIR)); sys.exit(-2)
except: parser.error("Unable to create a directory at \'%s\'" % (opts.SCRIPTS_DIR))

crisprSCRIPT_Path = path.realpath()

print(workingDirectory)
die

#Remove old hpc scrips
if opts.REMOVE_OLD: system("rm -f " +path.join(opts.SCRIPTS_DIR,"*.sh"))

#Paths
basepath = "/mnt/research/germs/shane/databases/crisprs"
dba_base = "/mnt/research/germs/shane/databases/assemblies/"
refDatabases = {"Corteva":".fasta","NCBI/refseq/archaea":".fasta","NCBI/refseq/bacteria":".fasta","NCBI/genbank/bacteria":'.fasta',"NCBI/genbank/archaea":'.fasta',"PATRIC2/fastas":".fna"}



#Get Todays date for the scripts
theDate = str(date.today())

#HPC Header
header = ""
for line in open("CRISPR_HPC_Header.sb"): header+=line

firstFileName = base+"/scripts/hpc_scripts/CRISPR_run_%i_%s.sb" % (0,theDate)
bashLog = open(firstFileName, 'w')
bashLog.write(header)
logNumber,count =0,0
validExts = set([".fna",".fasta"])

pCompletedAssemblies = pickle.load(open("PilerCR_Files.p","rb"))
mCompletedAssemblies = pickle.load(open("MinCED_Files.p","rb"))

for assembly_dir, ext in refDatabases.items():
    bashLog.write("cd "+ dba_base + assembly_dir+"\n")
    # print("Cleaning up")
    # system("rm -f %s/core* %s/*.tmp %s/*.gz" % (dba_base+assembly_dir,dba_base+assembly_dir,dba_base+assembly_dir))
    print("Getting assemblies from " + assembly_dir)
    assemblies = listdir(dba_base+assembly_dir)
    print("\t%i assemblies in %s" % (len(assemblies),assembly_dir))

    pilerOut = base + "/pilerCR/pat2"; mincedOut = base + "/minCED/pat2"
    if "refseq" in assembly_dir:  pilerOut = base + "/pilerCR/refseq"; mincedOut = base + "/minCED/refseq"
    if "genbank" in assembly_dir: pilerOut = base + "/pilerCR/genbank"; mincedOut = base + "/minCED/genbank"
    asmCount = 0
    for assembly in assemblies:
        if assembly[assembly.rfind("."):] not in validExts: continue
        #Add PilerCR command
        pcrPath = "%s/%s.pcrout" % (pilerOut, assembly.replace("%s"%(ext),""))
        if not pcrPath in pCompletedAssemblies and not path.exists(pcrPath):
            asmCount += 1
            pcmd = "pilercr -minid 0.85 -mincons 0.8 -noinfo -in %s -out %s/%s.pcrout 2>/dev/null" %(assembly, pilerOut, assembly.replace("%s" %(ext),""))
            count+=1
            bashLog.write(pcmd+"\n")
            if (count % 100 == 0):
                logNumber+=1
                bashLog.write("rm -f core* *.tmp *.gz\n")
                bashLog.close()
                bashLog = open(base+"/scripts/hpc_scripts/CRISPR_run_%i_%s.sb" % (logNumber,theDate), 'w')
                bashLog.write(header)
                bashLog.write("cd "+ dba_base + assembly_dir + "\n")
        else: pCompletedAssemblies.add(pcrPath)

        #Add minCED command
        miCDPath = "%s/%s.mnout" % (mincedOut, assembly.replace("%s" %(ext),""))
        if not miCDPath in mCompletedAssemblies and not path.exists(miCDPath): 
            asmCount += 1
            mcmd = "minced -minRL 16 -maxRL 64 -maxSL 64 -searchWL 6 %s %s/%s.mnout" %(assembly, mincedOut, assembly.replace("%s" %(ext),""))
            count+=1
            bashLog.write(mcmd+"\n")
            if (count % 100 == 0):
                logNumber+=1
                bashLog.write("rm -f core* *.tmp *.gz\n")
                bashLog.close()
                bashLog = open(base+"/scripts/hpc_scripts/CRISPR_run_%i_%s.sb" % (logNumber,theDate), 'w')
                bashLog.write(header)
                bashLog.write("cd "+ dba_base + assembly_dir + "\n")
        else: mCompletedAssemblies.add(miCDPath)
    print("\t\t%i assemblies from %s" % (asmCount,assembly_dir))

bashLog.close() 

pickle.dump(pCompletedAssemblies,open("PilerCR_Files.p","wb"))
pickle.dump(mCompletedAssemblies,open("MinCED_Files.p","wb"))

print("\nThere are %i logs for %i assemblies" % (logNumber,count))


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
