#Cock et al 2009. Biopython: freely available Python tools for computational
#molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
import sys
import re
try:
	from Bio.Alphabet import generic_dna
	from Bio.Seq import Seq, reverse_complement, translate
	from Bio.SeqRecord import SeqRecord
	from Bio.SeqFeature import SeqFeature, FeatureLocation
	from Bio import SeqIO
	from Bio.Data import CodonTable
except ImportError:
	sys.exit("Missing Biopython library")

class GetORFs:
	table_obj = CodonTable.ambiguous_generic_by_id[1]
	seq_format = "fasta"
	starts = sorted(table_obj.start_codons)
	re_starts = re.compile("|".join(starts))
	stops = sorted(table_obj.stop_codons)
	re_stops = re.compile("|".join(stops))

	def __init__(self, inputfile, outputfile, cutoff=50):
		self.inputFile = inputfile
		self.records = {}
		self.cutoff = cutoff
		self.outputFile = outputfile
		self.Seqlist = []
		self.writeOutPut()
		#self.data = []

	def start_chop_and_trans(s, strict=True):
		"""Returns offset, trimmed nuc, protein."""
		if strict:
			assert s[-3:] in stops, s
		assert len(s) % 3 == 0
		for match in re_starts.finditer(s):
			# Must check the start is in frame
			start = match.start()
			if start % 3 == 0:
				n = s[start:]
				assert len(n) % 3 == 0, "%s is len %i" % (n, len(n))
				if strict: t = translate(n, 1, cds=True)
				else: t = "M" + translate(n[3:], 1, to_stop=True) # Use when missing stop codon,
				return start, n, t
		return None, None, None
	
	def break_up_frame(self, s):
		"""Returns offset, nuc, protein."""
		start = 0
		for match in GetORFs.re_stops.finditer(s):
			index = match.start() + 3
			if index % 3 != 0: continue
			n = s[start:index]
		   
			offset = 0
			t = translate(n, 1, to_stop=True)
			if n and len(t) >= 8:
				yield start + offset, n, t
			start = index
		if "open" == "open":
			# No stop codon, Biopython's strict CDS translate will fail
			n = s[start:]
			# Ensure we have whole codons
			# TODO - Try appending N instead?
			# TODO - Do the next four lines more elegantly
			if len(n) % 3: n = n[:-1]
			if len(n) % 3: n = n[:-1]
			offset = 0
			t = translate(n, 1, to_stop=True)
			if n and len(t) >= 8:
				yield start + offset, n, t
	
	
	def get_all_peptides(self, nuc_seq):
		"""Returns start, end, strand, nucleotides, protein.
	
		Co-ordinates are Python style zero-based.
		"""
		# TODO - Refactor to use a generator function (in start order)
		# rather than making a list and sorting?
		answer = []
		full_len = len(nuc_seq)
		for frame in range(0, 3):
			for offset, n, t in self.break_up_frame(nuc_seq[frame:]):
				start = frame + offset  # zero based
				answer.append((start, start + len(n), +1, n, t))
		rc = reverse_complement(nuc_seq)
		for frame in range(0, 3):
			for offset, n, t in self.break_up_frame(rc[frame:]):
				start = full_len - frame - offset  # zero based
				answer.append((start - len(n), start, -1, n, t))
		answer.sort()
		return answer

	def writeOutPut(self):	
		#out_file = open(self.outputFile, "w")		
		in_count = 0
		out_count = 0		
		for record in SeqIO.parse(self.inputFile, GetORFs.seq_format, generic_dna):
			self.records[record.id] = record
			for i, (f_start, f_end, f_strand, n, t) in enumerate(self.get_all_peptides(str(record.seq).upper())):
				if len(t) < self.cutoff:continue
				out_count += 1
				if f_strand == +1: loc = "%i..%i" % (f_start + 1, f_end)
				else: loc = "complement(%i..%i)" % (f_start + 1, f_end)
				descr = "length %i aa, %i bp, from %s of %s" % (len(t), len(n), loc, record.description)
				fid = record.id + "_%s%i" % ("ORF", i + 1)
				t = SeqRecord(Seq(t), id=fid, name="", description=descr)
				ORF_Feature = SeqFeature(FeatureLocation(f_start + 1, f_end), type="CDS", strand=f_strand)
				ORF_Feature.qualifiers['label'] = [fid]
				self.records[record.id].features.append(ORF_Feature)
				self.Seqlist.append(t)
			in_count += 1
		self.write()   #nice_strand = '+' if f_strand == +1 else '-'
		# print "Found %i %ss in %i sequences" % (out_count, "ORF", in_count)

	def write(self):
		with open(self.outputFile,"w") as fh:
			SeqIO.write(self.Seqlist,fh,"fasta")
		
if __name__ == "__main__":
	fileName = sys.argv[1]
	outFile = sys.argv[2]
	cutoff = int(sys.argv[3])
	GetORFs(fileName, outFile, cutoff)
	#print ("Hello world")		
