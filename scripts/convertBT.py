import sys

fileName = sys.argv[1]
lines = []
for line in open(fileName):
    lines.append(line.strip())
    
ID = fileName.split("/")[-1].split(".")[0]
newFile = open(fileName,"w")
newFile.write(">%s\n" % ID)
newFile.write("%s\n" % lines[0])
newFile.write("%s\n" % (lines[1].split(" ")[0]))