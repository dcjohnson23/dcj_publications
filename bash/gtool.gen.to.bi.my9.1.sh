 #!/bin/bash
#BSUB -J "gentobigen1"
#BSUB -W 23:00

module load gtool/0.7.5
gtool -S --g /SNPtest/UK_sub_my9/UKGWAS_chr1.imp.my9.gen 
--s /SNPtest/UK_sub_my9/UKGWAS_chr1.imp.my9.sample 
--og /SNPtest/UK_ped_my9/UKGWAS_chr1.impbi.my9.gen 
--os SNPtest/UK_ped_my9/UKGWAS_chr1.impbi.my9.sample 
--inclusion SNPtest/UK_gen/bisnps_final_chr1.txt

