basedir = '/mnt/research/germs/shane/transActRNA/data'
logsDir = basedir+"/logs/blastLogs"
from glob import glob
logs = glob(logsDir+"/*.o")
from os import system, path
ids = {}
allParts = {}
complete = {}
blastFiles = {}
blastFilePath = "/mnt/research/germs/shane/transActRNA/data/annotations/blastout/%s.blastout"
from os import stat as fstat
for log in logs:
	if fstat(log).st_size == 0: continue

	splitID = log.replace(logsDir+"/","")
	splitID = splitID[:splitID.find(".NB.SPLIT")]
	protID = splitID[:splitID.rfind(".split")]
	parts = [0,0]
	for i,line in enumerate(open(log)): parts[i] = int(line.strip()[line.rfind(" "):])
	try: complete[protID][splitID] = parts[0]
	except:complete[protID] = {splitID:parts[0]}
	try: ids[protID] -= 1
	except: ids[protID] = parts[1] - 1
	allParts[protID]=parts[1]
	try: blastFiles[protID].add(blastFilePath % (splitID))
	except: blastFiles[protID] = set([blastFilePath % (splitID)])


for protID,numParts in allParts.items():
	if ids[protID] != 0:
		print(protID,"\t",ids[protID],numParts)
		continue
	if path.exists(blastFilePath % (protID)):continue
	with open(blastFilePath % (protID),'w') as fh:
		for splitFile in blastFiles[protID]:
			for line in open(splitFile): fh.write(line)
	res = system("rm %s/%s" % (basedir+'/assemblies/ids',protID))

	# Stuff not working and why
	ls -alh *.e | grep -v "0 Ju" | awk '{print $9}'| xargs cat