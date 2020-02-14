from os import system,path
from pickle import load
count = 0
assembliesWCrisprs = load(open("data/pickles/allAssemblyW_CRISPRs.p","rb"))
asmLinkPath = "/mnt/research/germs/shane/transActRNA/data/assemblies/assemblies_W_CRISPRs"
for asmName, asmPath in assembliesWCrisprs.items():
    cmd = "ln -s %s %s/%s.fasta" % (asmPath,asmLinkPath,asmName)
    if not path.exists("%s/%s.fasta" % (asmLinkPath,asmName)): 
        count+=1
        if count%1000 == 0:print(count,end=" ")
        system(cmd)    
print(count)