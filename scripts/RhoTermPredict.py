import numpy as np
import matplotlib.pyplot as plt
import re
from Bio.SeqIO import parse
from Bio.Seq import Seq
import openpyxl
np.warnings.filterwarnings('ignore')
with np.warnings.catch_warnings():
    np.warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')

# from:https://github.com/MarcoDiSalvo90/RhoTermPredict
def palindrome_control(sequenza,GC_genoma,strand):
	palindr = 0
	sequenza = Seq(sequenza)
	sequenza_r = sequenza.reverse_complement()
	for k in range (4,8):
		for i in range (150):
			x1 = i
			y1 = i + k
			if y1 > len(sequenza):
				break
			else:
				window = sequenza[x1:y1]
				window = str(window)
				gc_content = 100*(window.count("C")+window.count("G"))/len(window)
				sequenza_r = str(sequenza_r)
				sub_sequenza_r= sequenza_r
				par = 0
				while True:
					if re.search(window, sub_sequenza_r):
						ricerca = re.search(window, sub_sequenza_r)
						positions = ricerca.span()
						x2 = positions[0]
						y2 = positions[1]
						sub_sequenza_r=sub_sequenza_r[y2:]
						x2 += par
						y2 += par
						par += y2
						if 4 <= len(sequenza)-y2-y1+1 <= 8:
							palindr +=1
							break
					else: break
	return palindr

def palindrome_finder(sequenza,file,GC_genoma,strand,score):
	lista_scores = []
	sequenza = Seq(sequenza)
	sequenza_r = sequenza.reverse_complement()
	pattern_pause_site1 = "GG\D{8}[C,T]G"
	pattern_pause_site2 = "C[G,A]\D{8}CC" 
	for k in range (4,8):
		for i in range (150):
			x1 = i
			y1 = i + k
			if y1 > len(sequenza):
				break
			else:
				window = sequenza[x1:y1]
				window = str(window)
				gc_content = 100*(window.count("C")+window.count("G"))/len(window)
				sequenza_r = str(sequenza_r)
				score_p = score
				sub_sequenza_r= sequenza_r
				par = 0
				while True:
					if re.search(window, sub_sequenza_r):
						ricerca = re.search(window, sub_sequenza_r)
						positions = ricerca.span()
						x2 = positions[0]
						y2 = positions[1]
						sub_sequenza_r=sub_sequenza_r[y2:]
						x2 += par
						y2 += par
						par += y2
						if 4 <= len(sequenza)-y2-y1+1 <= 8:
							loop = len(sequenza)-y2 - y1
							file.write("\nPalindromic sequences found at coordinates ")
							file.write(str(x1))
							file.write("-")
							file.write(str(y1-1))
							file.write(" e ")
							file.write(str(len(sequenza)-y2))
							file.write("-")
							file.write(str(len(sequenza)-x2-1))
							file.write("(Sequence:   ")
							file.write(str(window))
							file.write(")")
							score_p += 3
							if gc_content > GC_genoma+20:
								score_p += 2
							elif gc_content > GC_genoma+10:
								score_p += 1
							if len(window) > 4:
								score_p += 1
							if loop < 6:
								score_p += 1
							if strand == 1:
								a = x1-5
								b = len(sequenza)-x2+5
								seq_pause = sequenza[a:b]
								seq_pause = str(seq_pause)
								if re.search(pattern_pause_site1,seq_pause):
									score_p += 3
									file.write(" (PAUSE-CONSENSUS PRESENT)")
							if strand == -1:
								a = x1-5
								b = len(sequenza)-x2+5
								seq_pause = sequenza[a:b]
								seq_pause = str(seq_pause)
								if re.search(pattern_pause_site2,seq_pause):
									score_p += 3
									file.write(" (PAUSE-CONSENSUS PRESENT)")
							file.write(" (SCORE: ")
							file.write(str(score_p))
							file.write(")")
							lista_scores.append(score_p)
					else:
						break
	if len(lista_scores) > 0:
		return np.max(lista_scores)
	else:
		return 0



# print("RhoTermPredict is a genome-wide predictor of transcription Rho-dependent terminators in bacterial genomes. It analyzes both the strands.\n\n")
# print("Input: Genome sequences file")
# print("Output: a xlsx file containing Rho-dependent terminators coordinates and a txt file containing informations about them")
# print("Genome file must be in fasta format")
# File_genome = input("Enter the input genome file name: ")

def GC_content(genome):
	gcCount, genomeLen = 0,0
	for rec in SeqIO.parse(genome, "fasta"):
		gcCount += rec.seq.count("G") + rec.seq.count("C")
		genomeLen += len(rec.seq)
	return gcCount/genomeLen*100

class RhoTermPredict:
	def __init__(self,subseqs,outfile,genome):
		Sequences = {}
		for rec in parse(subseqs,"fasta"): Sequences[rec.id] = str(rec.seq).upper()
		gcCount, genomeLen = 0,0
		for rec in parse(genome, "fasta"):
			gcCount += rec.seq.count("G") + rec.seq.count("C")
			genomeLen += len(rec.seq)
		GC_whole_genome = gcCount/genomeLen*100
		# try:
		# 	for seq_record in SeqIO.parse(File_genome, "fasta"): Sequences.append(str(seq_record.seq))
		# except IOError: print ("File %s inexistent in the current directory!" %(File_genome))
		# print ("RhoTermPredict is working, please wait...")
		# genome = Sequences[0] 
		# genome = genome.upper()
		# GC_whole_genome = 100*(genome.count("G")+genome.count("C"))/len(genome)
		pattern1 = "C\D{11,13}C\D{11,13}C\D{11,13}C\D{11,13}C\D{11,13}C"
		pattern2 = "G\D{11,13}G\D{11,13}G\D{11,13}G\D{11,13}G\D{11,13}G"
		pattern_pause_site1 = "GG\D{8}[C,T]G"
		pattern_pause_site2 = "C[G,A]\D{8}CC" 
		predictions = 0
		cg_value = []
		scores = []
		num = 1
		cod = 1
		
		nuovo_file = openpyxl.Workbook()
		sheet = nuovo_file.active
		sheet.title = "predictions"
		sheet["A1"] = "Region"
		sheet["B1"] = "Start RUT"
		sheet["C1"] = "End RUT"
		sheet["D1"] = "Strand"
		
		p = open("tmp/InformationsAboutRhoPredictions.txt","w") #TODO: Remove this superflious log file
		p.write("Sequences of predicted Rho-dependent terminators")
		for seqID, genome in Sequences.items():
			#positive strand
			scale = 0
			for index in range(len(genome)):
				x1 = scale + index
				x2 = scale + index + 78
				if x2 > len(genome): break
				else:
					region = genome[x1:x2]
					numG = region.count("G")
					if numG > 0: gcRatio = region.count("C")/numG
					else: numG = 1; gcRatio = (region.count("C")+1)/numG
					if gcRatio < 1: continue
					if re.search(pattern1,region):
						dati = []
						for g in range(50):
							if x2+g <= len(genome):
								region = genome[x1+g:x2+g]
								if region.count("G") > 0:
									numG = region.count("G")
									gcRatio = region.count("C")/numG
								else:
									numG = 1
									gcRatio = (region.count("C")+1)/numG
								if gcRatio > 1:
									if re.search(pattern1,region): dati.append(gcRatio)
									else: dati.append(0)
								else: dati.append(0)
							else: break
						maxP = np.argmax(dati)
						gcRatio = np.max(dati)
						x1 = x1 + maxP
						x2 = x2 + maxP
						scale = x2 - index - 1
						s = genome[x2:x2+150]
						region = genome[x1:x2]
						score = 3
						ctrl = palindrome_control(s,GC_whole_genome,1)
						if ctrl > 0 or re.search(pattern_pause_site1,s):
							sheet["A%d"%(num+1)] = seqID #"T%d"%(cod)
							sheet["B%d"%(num+1)] = x1
							sheet["C%d"%(num+1)] = x2
							sheet["D%d"%(num+1)] = "plus"
							num += 1
							predictions += 1
							scale = x2 + 150 - index - 1
							cg_value.append(gcRatio)
							if gcRatio > 2: score += 3
							elif gcRatio > 1.50: score += 2
							elif gcRatio > 1.25: score += 1
							p.write("\n\n\nPREDICTED REGION NUMBER ")
							p.write(str(cod))
							p.write(" (STRAND POSITIVE) ")
							p.write("\n\nGenomic sequence of putative RUT site (Coordinates:  ")
							p.write(str(x1))
							p.write("-")
							p.write(str(x2))
							p.write(", c/g = ")
							p.write(str(gcRatio))
							p.write(" ):	  ")
							p.write(str(region))
							p.write("\n\nThe 150 nt long genomic region immediately downstream is ")
							p.write(str(s))
							p.write("\n")
							cod += 1
							final_score = 0
							final_score = palindrome_finder(s,p,GC_whole_genome,1,score)
							par = 0
							while True:
								if re.search(pattern_pause_site1,s):
									ricerca = re.search(pattern_pause_site1,s)
									positions = ricerca.span()
									x2 = positions[0]
									y2 = positions[1]
									s = s[y2:]
									p.write("\n\nPAUSE-CONSENSUS present at the coordinates ")
									p.write(str(x2))
									p.write("-")
									p.write(str(y2))
									x2 += par
									y2 += par
								else: break
							if final_score == 0: final_score = score + 3
							scores.append(final_score)	
			#negative strand
			scale = 0
			for index in range(len(genome)):
				x1 = scale + index
				x2 = scale + index + 78
				if x2 > len(genome): break
				else:
					region = genome[x1:x2]
					if region.count("C") > 0:
						numC = region.count("C")
						gcRatio = region.count("G")/numC
					else:
						numC = 1
						gcRatio = (region.count("G")+1)/numC
					if gcRatio > 1:
						if re.search(pattern2,region):
							dati = []
							for g in range(50):
								if x2+g <= len(genome):
									region = genome[x1+g:x2+g]
									if region.count("C") > 0:
										numC = region.count("C")
										gcRatio = region.count("G")/numC
									else:
										numC = 1
										gcRatio = (region.count("G")+1)/numC
									if gcRatio > 1:
										if re.search(pattern2,region):
											dati.append(gcRatio)
										else: dati.append(0)
									else: dati.append(0)
								else: break
							maxP = np.argmax(dati)
							gcRatio = np.max(dati)
							x1 = x1 + maxP
							x2 = x2 + maxP
							scale = x2 - index - 1
							s = genome[x1-150:x1]
							region = genome[x1:x2]
							score = 3
							ctrl = palindrome_control(s,GC_whole_genome,-1)
							if ctrl > 0 or re.search(pattern_pause_site2,s):
								sheet["A%d"%(num+1)] = seqID #"T%d"%(cod)
								sheet["B%d"%(num+1)] = x1
								sheet["C%d"%(num+1)] = x2
								sheet["D%d"%(num+1)] = "minus"
								num += 1
								predictions += 1
								scale = x2 + 150 - index - 1
								cg_value.append(gcRatio)
								if gcRatio > 2: score += 3
								elif gcRatio > 1.50: score += 2
								elif gcRatio > 1.25: score += 1
								p.write("\n\n\nPREDICTED REGION NUMBER ")
								p.write(str(cod))
								p.write(" (STRAND NEGATIVE)")
								p.write("\n\nGenomic sequence of putative RUT site (Coordinates:  ")
								p.write(str(x1))
								p.write("-")
								p.write(str(x2))
								p.write(", c/g = ")
								p.write(str(gcRatio))
								p.write(" ):	  ")
								p.write(str(region))
								p.write("\n\nThe 150 nt long genomic region immediately downstream is ")
								p.write(str(s))
								p.write("\n")
								cod += 1
								final_score = 0
								final_score = palindrome_finder(s,p,GC_whole_genome,-1,score)
								par = 0
								while True:
									if re.search(pattern_pause_site2,s):
										ricerca = re.search(pattern_pause_site2,s)
										positions = ricerca.span()
										x2 = positions[0]
										y2 = positions[1]
										s = s[y2:]
										p.write("\n\nPAUSE-CONSENSUS present at the coordinates ")
										p.write(str(x2))
										p.write("-")
										p.write(str(y2))
										x2 += par
										y2 += par
									else: break
								if final_score == 0: final_score = score + 3
								scores.append(final_score)
	   
		p.write("\n\n\n\n\n\nTotal number of predicted Rho-dependent terminators: ")
		self.hasHits = (predictions > 0)
		self.predictions = predictions
		p.write(str(predictions))
		p.write("\n\n\nMean C/G content of predicted terminators: ")
		p.write(str(np.mean(cg_value)))
		p.write("\nStandard deviation: ")
		p.write(str(np.std(cg_value)))
		p.close()
		nuovo_file.save('%s.xlsx' % (outfile))
		# print("Work finished, see output files in the current directory")
			

	
if __name__ == "__main__":
	from sys import argv
	t = RhoTermPredict(argv[1],argv[2])
	