#!/usr/bin/python
# Josh Quick 2012

import sys
import vcf
import numpy
import argparse

def filter(args):
	vcf_reader = vcf.Reader(sys.stdin)
	orig_records = [record for record in vcf_reader]
	last_entry = orig_records[-1].POS
	print >> sys.stderr, len(orig_records), 'records in VCF'
	filtered_pos = numpy.zeros(last_entry+1, dtype=bool)

	#Discriminatory position filter
	for record in orig_records:
		samples = record.samples
		ref, alt, het = [0] * 3
		for each in samples:
			if each['GT'] == '0/0':	ref += 1
			if each['GT'] == '1/1':	alt += 1
			if each['GT'] == '0/1':	het += 1
		if ref > 0 and alt > 0 and het >= 0:
			filtered_pos[record.POS] = True

	print >> sys.stderr, numpy.count_nonzero(filtered_pos), 'records retained by discriminatory filter'
	
	#Allele frequency filter
	if args.use_het == True:
		max_het = args.het_freq
		min_het = 100 - max_het
	
		for n, record in enumerate(orig_records):
			samples = record.samples
			freqs = [float(sample['FREQ'].split('%')[0]) if not sample['FREQ'] == None else '' for sample in samples]
			if [freq for freq in freqs if freq is not None and (freq >= min_het and freq <= max_het)]:
				filtered_pos[record.POS] = False

		print >> sys.stderr, numpy.count_nonzero(filtered_pos), 'records retained by allele frequency filter'
	
	#Sliding window density filter
	if args.use_density == True:
		window = args.window
		max_density = args.density
		for n in xrange(last_entry):
			count = numpy.count_nonzero(filtered_pos[n:n + window])
			if count > max_density:
				for pos in numpy.nonzero(filtered_pos[n:n + window])[0].tolist():
					filtered_pos[pos + n] = False
	
		print >> sys.stderr, numpy.count_nonzero(filtered_pos), 'records retained by density filter'
	
	#Write out remaining records
	vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
	for record in orig_records:
		if filtered_pos[record.POS] == True:
			vcf_writer.write_record(record)	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='VCF filter to remove non-discriminatory positions, can also optionally remove positions with heterozygous calls or high SNP density ')
	parser.add_argument('--use-het-filter', dest='use_het', action='store_true', help='Use het filter')
	parser.add_argument('--max-freq-for-het', dest='het_freq', default='90', type=float, help='Maximum allele frequency for het')
	parser.add_argument('--use-density-filter', dest='use_density', action='store_true', help='Use density filter') 
	parser.add_argument('--window', dest='window', default='1000', type=int, help='Window size')
	parser.add_argument('--density', dest='density', default='3', type=int, help='Max density')
	args = parser.parse_args()
	filter(args)
