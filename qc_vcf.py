#!/usr/bin/env python
# Josh Quick 2014

#Script to perform multiple QC steps on your samples/variants
#Output information on SNPs and INDELs, some basic information about mapping/variant quality
#Count the nocalls like vcf2phyloviz
#Calculate the mean coverage of variant positions for each sample
#Use the callable script to plot variants vs coverage (bedtools?)
#Display the genotypes present in all samples 
#Plot the number of Nocalls vs homs for each sample
#Plot SNP density to visualise recombination/mis-mapped regions
#Option to analyse a subset of records to save time
#Bacterial specific

from collections import defaultdict
import argparse
import time
import random
import vcf
import numpy
import gzip
import zlib
import sys
import os
from subprocess import call
from cStringIO import StringIO
from operator import itemgetter
import csv

def get_chroms(args):
#	magic_bytes = '\x1f\x8b\x08'
#	vcf_in = open(args.file_in, 'rb')
#	file_start = vcf_in.read(3)
#	if file_start == magic_bytes:
#		print 'Detected gzip compression...'
#		vcf_in.close()
#		print 'Finding number of records, can take a while with large files...'
#		CHUNK_SIZE = 1024*1024
#		MAX_LINE = 50000	
#		decompress = zlib.decompressobj(-zlib.MAX_WBITS)
#		in_file = gzip.open(args.file_in, 'r')
#		in_file._read_gzip_header()
#		chunk = prior_chunk = ''
#		while 1:
#			buf = in_file.fileobj.read(CHUNK_SIZE)
#			if not buf:
#				break
#			d_buf = decompress.decompress(buf)
#			if d_buf:
#				prior_chunk = chunk
#				chunk = d_buf
#
#		if len(chunk) < MAX_LINE:
#			chunk = prior_chunk + chunk
#		
#		line = chunk[-MAX_LINE:].splitlines(True)[-1]
#		last_record = int(line.split('\t')[1])
#       else:
#               vcf_in.seek(-2, 2)
#               while vcf_in.read(1) != '\n':
#                       vcf_in.seek(-2, 1)
#               last_record = int(vcf_in.readline().split('\t')[1])
#               vcf_in.close()
#       return last_record
	
	chroms = []
	vcf_in = open(args.file_in, 'rb')
	current_cols = cached_cols = ''
	first_record = True
	for line in vcf_in:
		if line.startswith('#'):
			continue
		cached_cols = current_cols
		current_cols = line.split('\t')
		if not first_record:
			if current_cols[0] != cached_cols[0]:
				chroms.append({'name': cached_cols[0], 'last': int(cached_cols[1]) + 1})
		first_record = False
	chroms.append({'name': current_cols[0], 'last': int(current_cols[1]) + 1})
	return chroms
			
def iterate_all(args, chroms):
	# set up a dict of arrays with the right lengths
	filtered_pos = {}
	for chrom in chroms:
		filtered_pos[chrom['name']] = numpy.zeros(chrom['last']+1, dtype=bool)
	vcf_reader = vcf.Reader(open(args.file_in, 'r'))
	print '%i samples' %len(vcf_reader.samples)
	no_calls = defaultdict(int)
	het_calls = defaultdict(int)
	[no_calls[sample] for sample in vcf_reader.samples]
	[het_calls[sample] for sample in vcf_reader.samples]
	start = time.time()
	for num_records, record in enumerate(vcf_reader):
		samples = record.samples
		ref = alt = 0
		for sample in samples:
			if sample['GT'] == '0/0': ref += 1
			if sample['GT'] == '1/1': alt += 1
			if sample['GT'] == '0/1': het_calls[sample.sample] += 1
			if sample.gt_bases is None: no_calls[sample.sample] += 1
		if ref > 0 and alt > 0:
			filtered_pos[record.CHROM][record.POS] = True

	print '%s records in VCF' %(num_records+1)
	#print '%s records retained by discriminatory filter in %.2f seconds' %(numpy.count_nonzero(filtered_pos), time.time()-start)
	print 'sample,no_calls'
	for sample, num_nocalls in sorted(no_calls.iteritems(), key=itemgetter(1)):
		perc = float(num_nocalls) / num_records * 100
		print '%s,%.2f' %(sample, perc)
	print 'sample,het_calls'
	for sample, num_hetcalls in sorted(het_calls.iteritems(), key=itemgetter(1)):
		perc = float(num_hetcalls) / num_records * 100
		print '%s,%.2f' %(sample, perc)
	return filtered_pos

def sliding_window(args, chrom, filtered_pos):
	counts = []
	window = 1000
	for n in xrange(0, chrom['last'], window):
		counts.append((n, numpy.count_nonzero(filtered_pos[chrom['name']][n:n+window])))
	return counts

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='QC tools for bacterial variant calling')
	parser.add_argument('file_in', help='File in')
	args = parser.parse_args()
	chroms = get_chroms(args)
	filtered_pos = iterate_all(args, chroms)
	for chrom in chroms:
		print chrom
		counts = sliding_window(args, chrom, filtered_pos)
		with open('%s_density.csv' %(chrom['name'].split('|')[3][:-2]), "wb") as fh_out:
			#print 'chrom,position,counts'
			print >>fh_out, 'position,counts'
			for each in counts:
				#print '%s %5d %s' %(chrom['name'], each[0], '-'*(each[1]/4))
				print >>fh_out, '%i,%i' %(each[0], each[1])
		
