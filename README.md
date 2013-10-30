snp_calling_scripts
===================

Collection of scripts for bacterial SNP calling

Usage: cat merged.vcf | python filter_non_discriminatory_variants.py --use-density --window 1000 --max_snps 3 > merged_nondis.vcf

------

Usage: python count_snp_density.py -vcfs vcfs.txt

Where vcfs.txt contains a list of .vcf files.

------

Usage: cat merged.vcf | python output_af_from_vcf.py

This will output a table of allele frequencies, samples v positions extracted from a vcf file containing allele frequencies i.e. VarScan.

sample,1,2,3,4,5
sample_1,0,0,0,0,0
sample_2,0,0,0,0,0
