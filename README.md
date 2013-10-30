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

sample,14270,17152,18532,19568,19572,32507,36448,60760,64630,64706,70594,72619,76953,110085,
sample_900,0,3.85,0,100,100,0,0,76.19,97.87,97.83,0,100,0,0,
sample_901,0,0,0,100,100,0,0,72,100,100,0,100,0,0,
sample_902,0,0,0,100,100,0,0,60,100,95.24,0,100,0,0,
sample_903,0,0,0,100,100,0,0,73.33,94.55,94.64,0,100,0,0,
sample_910,0,0,0,100,100,0,0,70.97,100,97.14,0,100,0,0,
sample_916,0,0,0,100,100,0,0,53.57,88.57,93.88,0,100,0,0,
sample_917,0,0,0,100,100,0,0,68.42,96.08,88.64,0,100,0,0,
sample_930,0,0,0,100,100,0,0,69.44,94.87,95.83,0,100,0,0,
sample_931,0,0,0,100,100,0,0,76.19,93.75,79.31,0,100,0,0,
