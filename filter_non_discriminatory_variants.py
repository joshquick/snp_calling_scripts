#!/usr/bin/python
# Josh Quick 2012

import sys
import vcf
import numpy
import argparse
import time
#from cProfile import run
#from pstats import Stats

def filter(args):
	vcf_in = open(args.file_in, 'r')
	vcf_in.seek(-2, 2)			# Jump to the second last byte.
	while vcf_in.read(1) != '\n':		# Until EOL is found...
		vcf_in.seek(-2, 1)		# ...jump back the read byte plus one more
	last_entry = int(vcf_in.readline().split('\t')[1])
	filtered_pos = numpy.zeros(last_entry+1, dtype=bool)	
	
	#Discriminatory position filter
	vcf_in.seek(0)
	vcf_reader = vcf.Reader(vcf_in)
	discrim_records = []
	start = time.time()
	for i, record in enumerate(vcf_reader):
		samples = record.samples
		ref, alt = [0] * 2
		for sample in samples:
			if sample['GT'] == '0/0':	ref += 1
			if sample['GT'] == '1/1':	alt += 1
		if ref > 0 and alt > 0:
			discrim_records.append(record)
			filtered_pos[record.POS] = True

	print >> sys.stderr, '%s records in VCF' %i+1
	print >> sys.stderr, '%s records retained by discriminatory filter in %.3f' %(numpy.count_nonzero(filtered_pos), time.time()-start)

	#Allele frequency filter
	if args.use_het == True:
		start = time.time()
		max_het = args.het_freq
		min_het = 100 - max_het
		for record in discrim_records:
			freqs = [float(sample['FREQ'].split('%')[0]) for sample in record.samples if not sample['FREQ'] == None]
			if [freq for freq in freqs if (freq >= min_het and freq <= max_het)]:
				filtered_pos[record.POS] = False

		print >> sys.stderr, '%s records retained by allele frequency filter in %.3f' %(numpy.count_nonzero(filtered_pos), time.time()-start)
	
	#Sliding window density filter
	if args.use_density == True:
		start = time.time()
		window = args.window
		max_density = args.density
		for n in xrange(last_entry):
			count = numpy.count_nonzero(filtered_pos[n:n + window])
			if count > max_density:
				for pos in numpy.nonzero(filtered_pos[n:n + window])[0].tolist():
					filtered_pos[pos + n] = False
	
		print >> sys.stderr, '%s records retained by density filter in %.3f' %(numpy.count_nonzero(filtered_pos), time.time()-start)
	
	#Write out remaining records
	start = time.time()
	vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
	for record in discrim_records:
		if filtered_pos[record.POS] == True:
			vcf_writer.write_record(record)

	print >> sys.stderr, '%s records written out in %.3f seconds' %(numpy.count_nonzero(filtered_pos), time.time()-start)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='VCF filter to remove non-discriminatory positions, can also optionally remove positions with heterozygous calls or high SNP density ')
	parser.add_argument('-i', dest='file_in', help='File in')
	parser.add_argument('--use-het-filter', dest='use_het', action='store_true', help='Use het filter')
	parser.add_argument('--max-freq-for-het', dest='het_freq', default='90', type=float, help='Maximum allele frequency for het')
	parser.add_argument('--use-density-filter', dest='use_density', action='store_true', help='Use density filter') 
	parser.add_argument('--window', dest='window', default='1000', type=int, help='Window size')
	parser.add_argument('--density', dest='density', default='3', type=int, help='Max density')
	args = parser.parse_args()
	#run('filter(args)', 'filter_profile.stats')
	#stats = Stats('filter_profile.stats')
	filter(args)
