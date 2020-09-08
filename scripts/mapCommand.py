from optparse import OptionParser
from Bio.SeqIO import parse, write
from pickle import dump, load
import os

parser = OptionParser()

parser.add_option("-r", "--reverse",
                  action="store_true", dest="reverse", default=False,
                  help="Perform a reverse mapping of short ids to long IDs")

(opts, args) = parser.parse_args()

if __name__ == "__main__":
	#check to see if this is a mapping or a reverse mapping
	if opts.reverse:
		print "Reversing"
		# seqMap = load(open(args[0]+".map.p","rb"))

	else:
		print "Mapping"
		index = 0
		tmpfile = open("tmpfile.fasta","w")
		map = {}
		for seq in parse(args[0],"fasta"):
			newName = "Seq_%i" % (index)
			tmpfile.write(">%s\n%s\n" % (newName,str(seq.seq)))
			map[newName] = seq.id
			index +=1
		tmpfile.close()
		os.system("cp %s %s" % (args[0],args[0].replace(".fasta","_Original.fasta")))
		os.system("mv tmpfile.fasta %s" % (args[0]))
		# dump(map,open(args[0]+".map.p","wb"))
