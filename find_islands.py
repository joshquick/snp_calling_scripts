from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import re
import os

#Script to extract the gaps from a mauve xmfa file, currently only works on a two sample alignment

#Read the alignments into a list of seq records
alignments = []
content = open(sys.argv[1], 'r').read()
regex = re.compile(r'>([^\n]+)[\n]([ACGT\-\n]+)[\>]([^\n]+)[\n]([ACGT\-\n]+)', re.S)
matches = [m.groups() for m in regex.finditer(content)]
for m in matches:
	alignment = []
	seq = SeqRecord(Seq(m[1].replace(' ','').replace('\n', '')), id=m[0])
	alignment.append(seq)
	seq = SeqRecord(Seq(m[3].replace(' ','').replace('\n', '')), id=m[2])
	alignment.append(seq)
	alignments.append(alignment)


#Find gaps using a regex and write out the corresponding sequence
outh = sys.stdout #open('gaps.fasta', 'w')
regex = re.compile(r'([\-]{20,10000})')
for each in alignments:
	for m in regex.finditer(str(each[0].seq)):
		print 'gap found in sample 1', m.start(), m.end()
		outseq = SeqRecord(Seq(str(each[0].seq)[m.start():m.end()]), id=each[0].id, description='')
		SeqIO.write(outseq, outh, 'fasta')
		outseq = SeqRecord(Seq(str(each[1].seq)[m.start():m.end()]), id=each[1].id, description='')
		SeqIO.write(outseq, outh, 'fasta')
		print >>outh, '='
	for m in regex.finditer(str(each[1].seq)):
		print 'gap found in sample 2', m.start(), m.end()
		outseq = SeqRecord(Seq(str(each[1].seq)[m.start():m.end()]), id=each[1].id, description='')
		SeqIO.write(outseq, outh, 'fasta')
		outseq = SeqRecord(Seq(str(each[0].seq)[m.start():m.end()]), id=each[0].id, description='')
		SeqIO.write(outseq, outh, 'fasta')
		print >>outh, '='	
