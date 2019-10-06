#!/bin/bash
#BSUB -W 168:00
#BSUB -n 1
#BSUB -J bed.ger.1to4

module load plink/1.90.beta2k

for i in 1 2 3 4
do
   plink --gen /SNPtest/GER_bed/GerGWAS_chr${i}.imp.hdb.gen 
		--sample /SNPtest/GER_bed/GerGWAS_chr${i}.imp.hdb.sample 
		--oxford-single-chr ${i} --oxford-pheno-name plink_pheno 
		--maf 0.01 
		--geno 0.05 
		--hard-call-threshold 0.1 
		--keep /SNPtest/GER_bed/Ger.cases.GWAS.txt 
		--make-bed --out /SNPtest/GER_bed/GerGWAS_chr${i}.imp.hdb
done
