#!/usr/bin/python
# Josh Quick 2013

import vcf
import numpy
import argparse
import sys
import csv

def count(args):
	count_list = []
	win_out = []
	window = args.win
	files = [line.split('\n')[0] for line in open(args.vcfs, 'r')]
	for i, each in enumerate(files):
		records = [record for record in vcf.Reader(open(each, 'r'))]
		last_entry = records[-1].POS
		snp_positions = numpy.zeros(last_entry+1, dtype=bool)
		for record in records:
			snp_positions[record.POS] = True
		temp_counts = []
		for n in xrange(0, last_entry, window):
			count = numpy.count_nonzero(snp_positions[n:n+window])
			temp_counts.extend([count])
			if i == 0:
				win_out.extend([n])
		if i == 0:
			count_list.append(win_out)
		count_list.append(temp_counts)

	with open("test_count.csv", "wb") as fh_out:
		writer = csv.writer(fh_out)
		files.insert(0, 'position')
		writer.writerow(files)
		writer.writerows(zip(*count_list))	
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Count SNP density')
	parser.add_argument('-vcfs', dest='vcfs', help='List of vcf files to count')
        parser.add_argument('--win', dest='win', default='1000', type=int, help='Window size')
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit()
	args = parser.parse_args()
	count(args)
