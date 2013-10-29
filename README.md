snp_calling_scripts
===================

Collection of scripts for bacterial SNP calling

Usage: cat merged.vcf | python filter_non_discriminatory_variants.py --use-density --window 1000 --max_snps 3 > merged_nondis.vcf

Usage: python count_snp_density.py -vcfs vcfs.txt

Where vcfs.txt contains a list of .vcf files.
