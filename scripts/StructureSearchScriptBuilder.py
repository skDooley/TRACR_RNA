from os import chdir,system
from sys import argv
from glob import glob

def Main(runNum):
    print ("Running conseqs%i" % (runNum))
    BASE_DIR = "/mnt/research/germs/shane/transActRNA"
    SCRIPTS_DIR = "%s/scripts" % (BASE_DIR)
    DATA_DIR = "%s/data" % (BASE_DIR)
    LOGS = "/mnt/research/germs/shane/transActRNA/data/logs/structureSearches/"

    ####################################  Commands #################################### 
    MAP_CMD = "python ~/bin/mapCommand.py %s"
    MAFFT_CMD = "mafft --thread 20 --clustalout --maxiterate 1000 --globalpair %s.fasta > %s.aln 2>/dev/null"
    RNAailifold_CMD = "RNAalifold -r -d 2 < %s.aln > %s.bt"
    CONVERT_BT = "python %s/convertBT.py %s.bt"
    CREATE_STO = "perl %s/Bracket2Stockholm.pl ./"
    BUILD_CM = "cmbuild -F %s.cm %s.sto >/dev/null"
    CM_CALIBRATE = "cmcalibrate --cpu 20 %s.cm  >/dev/null"
    #CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/assemblies/Cas9_Representative_Assemblies.fasta > %s.out"
    if runNum == 0:   CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/assemblies/All_Cas9-Like-filtered.fasta &> %s.out"
    elif runNum == 1: CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/sequences/TracrRNA_PredictionOverlap.fasta &> %s.out"
    elif runNum == 3: CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/assemblies/UnknownTracrSystems.fasta &> %s.out"
    elif runNum == 4: CM_SEARCH = "cmsearch --cpu 1 %s.cm %s/assemblies/PioneerAssemblies.fasta &> %s.out"
    else:             CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/sequences/SgRNA_StructureResults.fasta &> %s.out"
    # CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/MasterTracrSeqs.fasta > %s.out"
    REV_CMD = "python ~/bin/mapCommand.py %s.map.p %s.out -r"
    FOLD_CMD = "perl %s/varna/HTP_RNAfold.pl %s.fasta > %s.dbn"
    GEN_IMAGE_CMD = "perl %s/varna/HTP_VARNA_images.pl %s.dbn" 
    chdir("/mnt/research/germs/shane/transActRNA/data/conseqs%i" % (runNum))
    cmOnly = False
    if len(argv)==3:
        cmOnly=True
        MAP_CMD = "#python ~/bin/mapCommand.py %s"
        MAFFT_CMD = "#mafft --thread 8 --clustalout --maxiterate 1000 --globalpair %s.fasta > %s.aln 2>/dev/null"
        RNAailifold_CMD = "#RNAalifold -r -d 2 < %s.aln > %s.bt"
        CONVERT_BT = "#python %s/convertBT.py %s.bt"
        CREATE_STO = "#perl %s/Bracket2Stockholm.pl ./"
        BUILD_CM = "#cmbuild -F %s.cm %s.sto >/dev/null"
        CM_CALIBRATE = "#cmcalibrate --cpu 8 %s.cm  >/dev/null"
        files = ["Cluster_%s" %(argv[2])]
    else:
        files = glob("*.fasta")
    files.sort()
    print ("There are %i files to process" % len(files))

    ####################################  Build Command Scripts ####################################
    HEADER = ""
    for line in open("/mnt/research/germs/shane/transActRNA/scripts/hpc/StructHeader.sb"): HEADER+=line
    if runNum == 0 and cmOnly: HEADER=HEADER.replace("time=1:30:00","time=2:00:00")
    if runNum == 1: HEADER=HEADER.replace("cpus-per-task=20","cpus-per-task=20").replace("time=1:30:00","time=0:30:00")
    if runNum == 4: HEADER=HEADER.replace("cpus-per-task=20","cpus-per-task=1").replace("time=1:30:00","time=0:30:00")

    
    
    sbScripts = []
    for i,fName in enumerate(files):
        id = fName.replace(".fasta","")
        scriptID = "/mnt/research/germs/shane/transActRNA/scripts/hpc/structSearches/%s.sb" % (id)
        with open(scriptID,'w') as fh:
            fh.write(HEADER)
            fh.write("cd %s/conseqs%i \n" % (DATA_DIR, runNum))
            # if not cmOnly:
            #     fh.write("mkdir %s \n" % (id))
            #     fh.write("mv %s %s/ \n" % (fName,id))
            fh.write("cd %s \n" % (id))

            name = fName.replace(".fasta","")
            cmds = [MAP_CMD %(fName),MAFFT_CMD % (name,name), RNAailifold_CMD % (name,name), CONVERT_BT % (SCRIPTS_DIR,name), CREATE_STO %(SCRIPTS_DIR), BUILD_CM % (name,name), CM_CALIBRATE % (name), CM_SEARCH % (name,DATA_DIR,name), REV_CMD % (name,name)] # ,FOLD_CMD %(SCRIPTS_DIR,name,name),GEN_IMAGE_CMD % (SCRIPTS_DIR,name)]
            for cmd in cmds: fh.write(cmd+'\n')
            fh.write("echo -en \"$(date) \tCompleted Run \"\n")
        sbScripts.append(scriptID)
    print("Done writing scripts.")
    if cmOnly:
        logsPath = "/mnt/research/germs/shane/transActRNA/data/logs/structureSearches"
        for sb in sbScripts:
            sName = sb.replace("/mnt/research/germs/shane/transActRNA/scripts/hpc/structSearches/","").replace(".sb","")
            system("sbatch --job-name=%s --error=%s/%s.out --output=%s/%s.out %s" % (sName,logsPath,sName,logsPath,sName,sb))
            print(sb)
          

if __name__ == "__main__":
    Main(int(argv[1]))

