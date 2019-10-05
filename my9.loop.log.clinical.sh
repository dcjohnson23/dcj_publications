#!/bin/bash
#BSUB -W 168:00
#BSUB -n 1
#BSUB -J cd.my9.lg

module load plink/1.07

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
   plink --bfile /scratch/haem/johnson/SNPtest/UK_bed_my9/UKGWAS_chr${i}.imputed.my9 --maf 0.01 --logistic --geno 0.05 --nonfounders --ci 0.95 --pheno /scratch/haem/johnson/Clin.demo/my9.UK/My9.clin.pheno.txt --all-pheno --out /scratch/haem/johnson/Clin.demo/my9.UK/chr${i}.my9 --noweb
done
