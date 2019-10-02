import sqlite3

path="/scratch/HLA/CASE/GER/GERGWAS_MM_HLA_imputed"
study="GGWAS_MM"
conn = sqlite3.connect(":memory:")
c=conn.cursor()

c.execute("create table fam_dat (sample_id text, pheno integer)")
with open("{}/{}_HLA_imputed.fam".format(path, study)) as fam_file:
    for line in fam_file:
        bits=line.rstrip().split()
        samp_id=bits[1] # rsid
        pheno=bits[5]
        c.execute("insert into fam_dat values (?,?)",(samp_id, pheno))
c.execute('create index fam_dat_idx on fam_dat (sample_id)')
conn.commit()

c.execute("create table bim_dat (rsid text, bp text)")
with open("{}/{}_HLA_imputed.bim".format(path, study)) as bim_file:
    for line in bim_file:
        bits=line.rstrip().split()
        samp_id=bits[1] # rsid
        pheno=bits[3]
        c.execute("insert into bim_dat values (?,?)",(samp_id, pheno))
c.execute('create index bim_dat_idx on bim_dat (rsid)')
conn.commit()

with open("{}/{}_HLA_imputed.bgl.gprobs".format(path, study)) as ifile, open("{}/{}_HLA_imputed.geno".format(path, study),"w") as ofile:
    for line in ifile:
        bits=line.rstrip().split()
        this_snp=bits[0]
        #print(this_snp)
        bp=c.execute("select bp from bim_dat where rsid=?", (this_snp,)).fetchone()[0]
        snp_dat=["---", this_snp, bp]+bits[1:]
        ofile.write(" ".join(snp_dat)+"\n")

allrows=c.execute("select * from fam_dat").fetchall()
with open("{}/{}_HLA_imputed.sample".format(path, study),"w") as ofile:
    ofile.write("ID_1 ID_2 missing phenotype\n0 0 0 B\n")
    for row in allrows:
        sample, phen = row
        samp_dat=[sample, sample, '0', str(phen-1)]
        ofile.write(" ".join(samp_dat)+"\n")
