### eQTL mapping
### Compute PC
python3 bin/compute_genotype_pcs.py GWAS_eQTL.vcf.gz --keep --output_dir PC

### Prepare input file
python3 bin/eqtl_phenotype_preparation.py TPM.txt count.txt Gallus_gallus.GRCg6a.106.gtf sample_participants.list chicken.chrs.list pheno
# expression.tmp_new.txt: TPM matrix
# expression.counts_new.txt : raw counts of genes
# gallus.gtf: annotation file in gtf format
# sample_participants.list: a file containing sample list
# chicken.chrs.list: a file containing chromosomes to be analyzed
# pheno: prefix of output

### Compute peer factors
Rscript bin/run_PEER.R TMM.txt peer 60 --output_dir peer
# TMM.txt: output file exclueded the chromosome, start and end columns from last step

# Parse parameters
cis_output=${outpath}/ChickenGTEx.cis_qtl_mapping_step2.txt
## Creat output directory
mkdir -p ${outpath}
prefix=${outpath}/ChickenGTEx
# 1) cis-QTL mapping: compute cis nominal associations for all variant-gene pairs
python3 -m tensorqtl \
            ${genotype} ${phenotype} ${prefix} \
            --mode cis_nominal \
            --covariates ${covariates} 
# 2) cis-QTL mapping: compute cis permutations and define eGenes
python3 -m tensorqtl \
            ${genotype} ${phenotype} ${prefix} \
            --mode cis \
            --covariates ${covariates}
# 3) define eGenes						
Rscript bin/run_tensorqtl_step2.R output/ChickenGTEx.cis_qtl.txt.gz output/ChickenGTEx.cis_qtl_mapping_step2.txt
# 4) cis-QTL mapping: conditionally independent QTLs
python3 -m tensorqtl \
            ${genotype} ${phenotype} ${prefix} \
            --mode cis_independent \
            --covariates ${covariates} \
            --cis_output ${cis_output}
# 5) trans-QTL mapping
python3 -m tensorqtl \
           ${genotype} ${phenotype} ${prefix} \
           --covariates ${covariates} \
           --mode trans
echo "Everything is done!"			
python3 bin/run_tensorqtl_step4_all.py





### SMR
cd ${workPath}
## Creat output directory
mkdir -p ${myqtl_path} ${qtl_bed}

#### Prepare eqtl file 
zcat ${qtl_path}/ChickenGTEx.Liver.cis_qtl_pairs.csv.gz > ${myqtl_path}/Liver.cis_qtl_pairs.txt
cat ${myqtl_path}/Liver.cis_qtl_pairs.txt | awk -F ' '  '$7>0{print $0}' > ${myqtl_path}/Liver.cis_qtl_pairs_2.txt
cat ${myqtl_path}/Liver.cis_qtl_pairs_2.txt | awk 'NR >= 2 {print $1,$2,$3,$7,$8}' > ${myqtl_path}/ChickenGTEx.Liver.cis_qtl_pairs.txt

#### Make a BESD file from fastqtl output 
${smr_path}/smr_Linux --eqtl-summary ${myqtl_path}/ChickenGTEx.Liver.cis_qtl_pairs.txt --fastqtl-nominal-format --make-besd --out ${myqtl_path}/ChickenGTEx.Liver.cis_qtl_pairs

#### Make update file 
join -1 1 -2 2 <(cat ${myqtl_path}/Liver.cis_qtl_pairs_2.txt | awk -v OFS="\t" '{if (NR >= 2) print $2, $1, $4, $8}' | sort -k1,1) <(cat ${genotype_plink}.bim | sort -k2,2) | sed 's/ /\t/g' | awk -v OFS="\t" '{print $1,$8,$9,$3}' | sort | uniq > ${myqtl_path}/Liver.cis.freq
join -1 1 -2 2 <(cat ${myqtl_path}/Liver.cis_qtl_pairs_2.txt | awk -v OFS="\t" '{if (NR >= 2) print $2, $1, $4, $8}' | sort -k1,1) <(cat ${genotype_plink}.bim | sort -k2,2) | sed 's/ /\t/g' | awk -v OFS="\t" '{print $5,$1,$6,$7,$8,$9}' | sort | uniq > ${myqtl_path}/Liver.cis.esi
join -1 1 -2 2 <(cat ${myqtl_path}/Liver.cis_qtl_pairs_2.txt | awk -v OFS="\t" '{if (NR >= 2) print $1}' | sort -k1,1) <(cat ${all_gene} | awk -v OFS="\t" '{if (NR >= 2) print $1,$2,$3,$4,$5,$6}' | sort -k2,2) | sed 's/ /\t/g' | awk -v OFS="\t" '{print $2, $1, $3, $4, $5, $6}' | sort | uniq  > ${myqtl_path}/Liver.cis.epi

####Update coordinates of SNPs and genes, frequency of effect allele 更新esi，epi，freq
${smr_path}/smr_Linux --beqtl-summary ${myqtl_path}/ChickenGTEx.Liver.cis_qtl_pairs --update-esi ${myqtl_path}/Liver.cis.esi
${smr_path}/smr_Linux --beqtl-summary ${myqtl_path}/ChickenGTEx.Liver.cis_qtl_pairs --update-epi ${myqtl_path}/Liver.cis.epi
${smr_path}/smr_Linux --beqtl-summary ${myqtl_path}/ChickenGTEx.Liver.cis_qtl_pairs --update-freq ${myqtl_path}/Liver.cis.freq

###cis
yhrun -N 1 -n 1 -c 20 ${smr_path}/smr_Linux --bfile ${genotype_plink} \
--gwas-summary xyz_AFP.txt \
--beqtl-summary ${myqtl_path}/ChickenGTEx.Liver.cis_qtl_pairs \
--diff-freq-prop 0.5 \
--heidi-mtd 1 \
--peqtl-smr 1e-6 \
--out ${outPath}/xyz_AFP_cis
###trans
yhrun -N 1 -n 1 -c 20 ${smr_path}/smr_Linux --bfile ${genotype_plink} \
--gwas-summary xyz_AFP.txt \
--beqtl-summary ${myqtl_path}/ChickenGTEx.Liver.trans_qtl_pairs \
--diff-freq-prop 0.5 \
--heidi-mtd 1 \
--peqtl-smr 1e-6 \
--trans \
--peqtl-trans 3.3156e-7 \
--out ${outPath}/xyz_AFP_trans