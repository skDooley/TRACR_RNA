import os
import time



def Main(runNum):
    print "Running conseqs%i" % (runNum)
    MAP_CMD = "python ~/bin/mapCommand.py %s"
    MAFFT_CMD = "mafft --clustalout --maxiterate 1000 --globalpair %s.fasta > %s.aln"
    RNAailifold_CMD = "RNAalifold -r -d 2 < %s.aln > %s.bt"
    CONVERT_BT = "python /mnt/research/germs/shane/transActRNA/scripts/convertBT.py %s.bt"
    CREATE_STO = "perl /mnt/research/germs/shane/transActRNA/scripts/Bracket2Stockholm.pl ./"
    BUILD_CM = "cmbuild %s.cm %s.sto"
    CM_CALIBRATE = "cmcalibrate %s.cm"
    CM_SEARCH = "cmsearch %s.cm /mnt/research/germs/shane/transActRNA/data/sequences/CasRelatedAssemblies.fa > %s.out"
    REV_CMD = "python ~/bin/mapCommand.py %s.map.p %s.out -r"
    
    os.chdir("/mnt/research/germs/shane/transActRNA/data/conseqs%i" %(runNum))
    files = os.listdir("./")
    print files
    for file in files:
        if ".fasta" in file and os.path.exists(file):
            id = file.replace(".fasta","")
            os.system("mkdir %s" % (id))
            os.system("mv %s %s/%s" % (file,id,file))
            os.chdir("%s" % (id))
            name = file.replace(".fasta","")
            cmds = [MAP_CMD %(file),MAFFT_CMD % (name,name), RNAailifold_CMD % (name,name), CONVERT_BT % (name), CREATE_STO, BUILD_CM % (name,name), CM_CALIBRATE % (name), CM_SEARCH % (name,name), REV_CMD % (name,name)]
            for cmd in cmds:
                print "\t",cmd
                res = os.system(cmd)
                if res != 0: print "\tFailed:",cmd
                time.sleep(3)
            os.chdir("/mnt/research/germs/shane/transActRNA/data/conseqs%i" %(runNum))

if __name__ == "__main__":
    Main(0)
