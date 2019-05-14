from optparse import OptionParser
from os import system
from Bio.SeqIO import parse
from Bio import AlignIO

parser = OptionParser()

parser.add_option("-i", "--in", dest="IN_NAME", help="The stockholm formatted consensus structure file")
parser.add_option("-o", "--out", dest="OUT_NAME", help="The jpg file output name")

(opts, args) = parser.parse_args()

VARNA_CMD = "java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd  \
-algorithm radiate -baseOutline \"#FFFFFF\" \ -baseInner \"#FFFFFF\" \
-sequenceDBN \"%s\" -structureDBN \"%s\" -o \"%s.jpg\""

print "Hello"
align = AlignIO.read(opts.IN_NAME, "stockholm")
print "Good bye"
for rec in align:
	print rec
