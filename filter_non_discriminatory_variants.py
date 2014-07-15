#!/usr/bin/python
# Josh Quick 2012

import sys
import vcf
from pysam import VCF
import numpy
import argparse
import time
#from cProfile import run
#from pstats import Stats

def filter(args):
	# get the chromosomes and last entry
	vcf_in = open(args.file_in, 'r')
	chroms = []
	current_cols = cached_cols = ''
	first_record = True
	for line in vcf_in:
		if line.startswith('#'):
			continue
		cached_cols = current_cols
		current_cols = line.split('\t')
		if not first_record:
			if current_cols[0] != cached_cols[0]:
				chroms.append({'name': cached_cols[0], 
					      'last': int(cached_cols[1]) + 1})
		first_record = False
	chroms.append({'name': current_cols[0], 'last': int(current_cols[1]) + 1})
	
	# discriminatory position filter
	filtered_pos = {}
	if args.exclude:
		exclude_list = [line.strip() for line in open(args.exclude, 'r')]
	else:
		exclude_list = []
	for chrom in chroms:
		filtered_pos[chrom['name']] = numpy.zeros(chrom['last']+1, dtype=bool)
	vcf_in.seek(0)
	vcf_reader = vcf.Reader(vcf_in)
	#records = list(vcf_reader)
	start = time.time()
	max_het = args.het_freq
	min_het = 100 - max_het
	for i, record in enumerate(vcf_reader):
		record.samples = [sample for sample in record.samples if sample.sample not in exclude_list]
		ref = alt = 0
		for sample in record.samples:
			if sample['GT'] == '0/0': ref += 1
			if sample['GT'] == '1/1': alt += 1
			freqs = [float(sample['FREQ'].split('%')[0]) for sample in \
                                 record.samples if not sample['FREQ'] == None]
		if ref > 0 and alt > 0:
			filtered_pos[record.CHROM][record.POS] = True
			if args.use_het:
				if [freq for freq in freqs if (freq >= min_het and freq <= max_het)]:
					filtered_pos[record.CHROM][record.POS] = False
		if i % 1000 == 0:
                        print >> sys.stderr, '%i lines processed' %i

	print >> sys.stderr, '%i records in VCF' %(i+1)
	print >> sys.stderr, '%i samples removed prior to applying filters' %(len(exclude_list))
	print >> sys.stderr, '%i records retained by discriminatory filter in %.2f seconds' \
			      %(sum([numpy.count_nonzero(filtered_pos[chrom['name']]) for chrom in chroms]), time.time()-start)

	# sliding window density filter
	if args.use_density == True:
		start = time.time()
		win = args.window
		max_density = args.density
		for chrom in chroms:
			for n in xrange(chrom['last']):
				count = numpy.count_nonzero(filtered_pos[chrom['name']][n:n + win])
				if count > max_density:
					for pos in numpy.nonzero(filtered_pos[chrom['name']][n:n + win])[0].tolist():
						filtered_pos[chrom['name']][pos + n] = False
		print >> sys.stderr, '%i records retained by density filter in %.2f seconds' \
				      %(sum([numpy.count_nonzero(filtered_pos[chrom['name']]) for chrom in chroms]), time.time()-start)
	
	# write out remaining records
	start = time.time()
	vcf_reader.samples = [sample for sample in vcf_reader.samples if sample not in exclude_list]
	vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
	vcf_in.seek(0)
	vcf_reader = vcf.Reader(vcf_in)
	i = 0
	for record in vcf_reader:
		record.samples = [sample for sample in record.samples if sample.sample not in exclude_list]
		if filtered_pos[record.CHROM][record.POS] == True:
			vcf_writer.write_record(record)
			i += 1
		else:
			continue

	print >> sys.stderr, '%i records written out in %.2f seconds' \
			      %(i, time.time()-start)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='VCF filter to remove non- \
					 discriminatory positions, can also optionally \
					 remove positions with heterozygous calls or high \
					 SNP density')
	parser.add_argument('file_in', help='File in')
	parser.add_argument('--use-het', dest='use_het', action='store_true', 
			    help='Use het filter')
	parser.add_argument('--max-freq', dest='het_freq', default='90', type=float, 
			    help='Maximum allele frequency for het')
	parser.add_argument('--use-density', dest='use_density', action='store_true', 
			    help='Use density filter') 
	parser.add_argument('--window', dest='window', default='1000', type=int, 
			    help='Window size')
	parser.add_argument('--density', dest='density', default='3', type=int, 
			    help='Max density')
	parser.add_argument('--exclude', dest='exclude', type=str,
                            help='File containing samples to exclude')
	args = parser.parse_args()
	#run('filter(args)', 'filter_profile.stats')
	#stats = Stats('filter_profile.stats')
	filter(args)
