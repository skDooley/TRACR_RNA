import os
import time



def Main(runNum,times=0):
    print ("Running conseqs%i" % (runNum))
    BASE_DIR = "/mnt/research/germs/shane/transActRNA"
    SCRIPTS_DIR = "%s/scripts" % (BASE_DIR)
    DATA_DIR = "%s/data" % (BASE_DIR)

    ####################################  Commands #################################### 
    MAP_CMD = "python ~/bin/mapCommand.py %s"
    MAFFT_CMD = "mafft --thread 20 --clustalout --maxiterate 1000 --globalpair %s.fasta > %s.aln 2>/dev/null"
    RNAailifold_CMD = "RNAalifold -r -d 2 < %s.aln > %s.bt"
    CONVERT_BT = "python %s/convertBT.py %s.bt"
    CREATE_STO = "perl %s/Bracket2Stockholm.pl ./"
    BUILD_CM = "cmbuild -F %s.cm %s.sto >/dev/null"
    CM_CALIBRATE = "cmcalibrate --cpu 20 %s.cm  >/dev/null"
    #CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/assemblies/Cas9_Representative_Assemblies.fasta > %s.out"
    CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/assemblies/All_Cas9-Like-filtered.fasta > %s.out"
    #CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/sequences/tracrRNAStructureResultsConfirmed.fasta > %s.out"
    # CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/sequences/All_Predicted_TracrRNAs2.fasta > %s.out"
    # CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/sequences/All_predictedSgRNAs_3-17-2020.fasta > %s.out"
    # CM_SEARCH = "cmsearch --cpu 20 %s.cm %s/MasterTracrSeqs.fasta > %s.out"

    

    REV_CMD = "python ~/bin/mapCommand.py %s.map.p %s.out -r"
    FOLD_CMD = "perl %s/varna/HTP_RNAfold.pl %s.fasta > %s.dbn"
    GEN_IMAGE_CMD = "perl %s/varna/HTP_VARNA_images.pl %s.dbn" 

    ####################################  Process Commands ####################################
    os.chdir("%s/conseqs%i" % (DATA_DIR, runNum))
    from glob import glob
    files =glob("*.fasta")
    files.sort()
    print ("There are %i files to process" % len(files))
    from random import randint as r
    for i,file in enumerate(files):
        time.sleep(r(1,5))
        if os.path.exists(file):
            id = file.replace(".fasta","")
            os.system("mkdir %s" % (id))
            os.system("mv %s %s/%s" % (file,id,file))
            os.chdir("%s" % (id))
            name = file.replace(".fasta","")
            cmds = [MAP_CMD %(file),MAFFT_CMD % (name,name), RNAailifold_CMD % (name,name), CONVERT_BT % (SCRIPTS_DIR,name), CREATE_STO %(SCRIPTS_DIR), BUILD_CM % (name,name), CM_CALIBRATE % (name), CM_SEARCH % (name,DATA_DIR,name), REV_CMD % (name,name)] #,FOLD_CMD %(SCRIPTS_DIR,name,name),GEN_IMAGE_CMD % (SCRIPTS_DIR,name)]
            for cmd in cmds:
                print ("\t",cmd)
                res = os.system(cmd)
                if res != 0: 
                    print ("\tFailed:",cmd)
                    time.sleep(1)
                    break
            os.chdir("%s/conseqs%i" % (DATA_DIR, runNum))
            if times!=0 and i+1 == times: break
    os.chdir("%s/conseqs%i" % (DATA_DIR, runNum))
print("Done")            

if __name__ == "__main__":
    from sys import argv
    if len(argv) >= 3:
        print("Running %i times" % (int(argv[2])))
        Main(int(argv[1]),int(argv[2]))
    else: Main(int(argv[1]))

