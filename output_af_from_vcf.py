#!/usr/bin/python
# Josh Quick 2013

import vcf
import sys
import csv

def main():
	allele_freq = []
	vcf_reader = vcf.Reader(sys.stdin)
	records = [record for record in vcf_reader]
	positions = [record.POS for record in records]
	for i, record in enumerate(records):
		samples = record.samples
		if i == 0:
			allele_freq.append([sample.sample for sample in samples])
		freqs = [sample['FREQ'].split('%')[0] if not sample['FREQ'] == None else '' for sample in samples]
		allele_freq.append(freqs)

	with open("allele_freqs.csv", "wb") as fh_out:
		writer = csv.writer(fh_out)
		positions.insert(0, 'sample')
		writer.writerow(positions)
		writer.writerows(zip(*allele_freq))

if __name__ == '__main__':
	main()
