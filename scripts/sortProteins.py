import re
import sys
sys.path.append("scripts/")

from Bio.SeqIO import write, parse
from collections import Counter
from CRISPRtools import CasOperon
from easyFunctions import *
from HMMParser import *
from IPython.display import display, HTML
from matplotlib import pyplot as plt
from os import chdir, listdir
from pandas import Series
from pickle import load #Dump has been wrapped aournd a function from easyFunctions

chdir("data")
gene = "Cas9"
geneProfile = "proteins/DiverseCas9s.faa" #File containing the seed proteins for the hmmsearch
geneFile = "proteins/%s.faa" % gene
alnName  = "alignments/%s.aln" % gene
hmmName  = "hmm/%s.hmm" % gene
hmmResultsDir = "hmm/results"
baseDbsDir = "/mnt/research/germs/shane/databases/assemblies/"
refDatabases = ["NCBI/refseq/bacteria","NCBI/refseq/archaea","NCBI/genbank/bacteria","NCBI/genbank/archaea","PATRIC2/fastas"]


casOperons = load(open("pickles/casOperonDataStructure.p","rb")) #all this has is the crisprs
casOperons.hasCas9(hmmResultsDir+"/")
casOperons.saveProgress()
dump(casOperons, open("pickles/casOperonDataStructureW%s.p" % gene ,"wb")) 