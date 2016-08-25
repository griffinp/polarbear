# polarbear
This repository contains supporting material for the EMBL-ABR bioinformatics tutorial 'Using genomic markers to study population divergence'. 

This tutorial uses data generated as part of the study 'Liu et al. (2014) Population genomics reveal recent speciation and rapid evolutionary adaptation in polar bears. Cell 157(4):785-794'. doi: 10.1016/j.cell.2014.03.054
and published in the online repository GigaDB at http://gigadb.org/dataset/view/id/100008/

The document 'convert_bear_snp_data_to_vcf.R' in this repository performs the following tasks using the files 'polar_bear.pooled.snp.txt' and 'brown_bear.pooled.snp.txt':
- extracts only SNPs from scaffold37 of the polar bear genome assembly
- retains only the first three individuals of each species
- outputs SNP information in vcf format for each species
- performs a simple merge (intersecting SNP locations only) on the two files and returns the result as vcf

