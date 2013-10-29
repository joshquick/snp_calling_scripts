#!/usr/bin/python
# Josh Quick 2012

import sys
import vcf
import numpy
import argparse

##samtoolsVersion=0.1.18 (r982:295)
##source=VarScan2

def filter(args):
	vcf_reader = vcf.Reader(sys.stdin)
	orig_records = [record for record in vcf_reader]
	last_entry = orig_records[-1].POS
	print >> sys.stderr, len(orig_records), 'records in VCF'
	filtered_pos = numpy.zeros(last_entry+1, dtype=bool)
	max_hets = args.hets
	max_ncs = args.ncs
	for record in orig_records:
		samples = record.samples
		ref, alt, het, nc = [0] * 4
		for each in samples:
			if each['GT'] == '0/0':	ref += 1
			if each['GT'] == '1/1':	alt += 1
			if each['GT'] == '0/1':	het += 1
			if each['GT'] == 'None': nc += 1

		if ref > 0 and alt > 0 and het < max_hets + 1:
			filtered_pos[record.POS] = True

	print >> sys.stderr, numpy.count_nonzero(filtered_pos), 'records retained by discrim filter'
	
	if args.use_density == True:
		window = args.window
		max_density = args.density
		for n in xrange(last_entry):
			count = numpy.count_nonzero(filtered_pos[n:n + window])
			if count > max_density:
				for pos in numpy.nonzero(filtered_pos[n:n + window])[0].tolist():
					filtered_pos[pos + n] = False
	
		print >> sys.stderr, numpy.count_nonzero(filtered_pos), 'records retained by density filter'
	
	vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
	for record in orig_records:
		if filtered_pos[record.POS] == True:
			vcf_writer.write_record(record)	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='VCF filter')
	parser.add_argument('--hets', dest='hets', default='3', type=int, help='Max hets')
	parser.add_argument('--nocalls', dest='ncs', default='3', type=int, help='Max nocalls')
	parser.add_argument('--use-density', dest='use_density', action='store_true', help='Use density filter') 
	parser.add_argument('--window', dest='window', default='1000', type=int, help='Window size')
	parser.add_argument('--max_snps', dest='density', default='3', type=int, help='Max density')
	args = parser.parse_args()
	filter(args)
