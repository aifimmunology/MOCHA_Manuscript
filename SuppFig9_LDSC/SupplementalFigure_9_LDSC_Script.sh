#!/bin/bash

source activate './ldsc/ldsc'


##### Pull down reference files
wget  https://zenodo.org/records/8292725/files/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
wget https://zenodo.org/records/8292725/files/1000G_Phase3_EAS_weights_hm3_no_MHC.tgz
tar -xvzf 1000G_Phase3_EAS_weights_hm3_no_MHC.tgz

wget https://zenodo.org/records/8292725/files/1000G_Phase3_frq.tgz
tar  -xvzf 1000G_Phase3_frq.tgz

wget https://zenodo.org/records/8292725/files/1000G_Phase3_plinkfiles.tgz
tar -xvzf 1000G_Phase3_plinkfiles.tgz

wget https://zenodo.org/records/8292725/files/hm3_no_MHC.list.txt



################ Download GWAS summary stats and transform it into a sumstats format
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90029001-GCST90030000/GCST90029015/GCST90029015_buildGRCh37.tsv

../../ldsc/munge_sumstats.py \
--sumstats GCST90029015_buildGRCh37.tsv \
--merge-alleles w_hm3.snplist \
--snp variant_id \
--chunksize 500000 \
--a1 effect_allele \
--a2 other_allele \
--signed-sumstats beta,0 \
--out sumstats/Autoimmune 

wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90043001-GCST90044000/GCST90043674/GCST90043674_buildGRCh37.tsv.gz

gunzip GCST90043674_buildGRCh37.tsv.gz

head GCST90043674_buildGRCh37.tsv

../../ldsc/munge_sumstats.py \
--sumstats GCST90043674_buildGRCh37.tsv \
--merge-alleles w_hm3.snplist \
--chunksize 500000 \
--snp variant_id \
--a1 effect_allele \
--a2 other_allele \
--signed-sumstats beta,0 \
--out sumstats/Abnormal_Immune 


wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90038001-GCST90039000/GCST90038603/GCST90038603_buildGRCh37.tsv

head GCST90038603_buildGRCh37.tsv


awk '$1 ~ "rs" || $1 ~ "variant_id" {print $0}' GCST90038603_buildGRCh37.tsv > GCST90038603.tsv

tr ' ' \\t < GCST90038603_hg19.tsv > GCST90038603_v2.tsv

head GCST90038603_hg19.tsv

cat GCST90038603_hg19.tsv | wc -l

../../ldsc/munge_sumstats.py \
--sumstats GCST90038603.tsv \
--N 9367966 \
--merge-alleles w_hm3.snplist \
--chunksize 500000 \
--snp variant_id \
--a1 effect_allele \
--a2 other_allele \
--signed-sumstats beta,0 \
--out Immune_Aging



############# Now process bed files from the Figure 2 Source Data

### Using the Soure Data for Figure 2
#### Split out the bed files by method. 
for XX in $(ls | grep '.bed'); do
    awk '$4 ~ "Common" || $4 ~ "MOCHA" {print $0}' $XX > "MOCHA_${XX}";
    awk '$4 ~ "Common" || $4 ~ "HOMER" {print $0}' $XX > "HOMER_${XX}";
    awk '$4 ~ "Common" || $4 ~ "MACS2" {print $0}' $XX > "MACS2_${XX}";
done

## Create a new sub folder for bed files by Chromosome for generating annotation files
## Then create bed files for each sample by chromosome
for bedFile in $(ls | grep '.bed');
do
    for chr in {1..22};
    do 
         awk -v name="chr$chr" '$1 == name {print $0}' $bedFile > "BedSubsets/${bedFile//.bed/.}${chr}.bed";         
  done;
done


## Running a for loop in order to create LDSC annot files for all bed files
for bedFile in $(ls | grep '.bed' | grep 'all_peaks');
do
    echo $bedFile
    for chr in {1..22};
    do 
         python ../../../ldsc/make_annot.py \
    		--bed-file "BedSubsets/${bedFile//.bed/.}${chr}.bed" \
    		--bimfile ../1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
    		--annot-file ../LDScores/${bedFile//.bed/.}${chr}.annot.gz &
    done;
done

######################### Now run LDSC regression

cd bedFiles

for bedFile in $(ls | grep '.bed');
do
    for chr in {1..22};
    do 
         python ../../../ldsc/ldsc.py \
    		--l2 \
    		--bfile ../1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    		--ld-wind-cm 1 \
    		--thin-annot \
    		--annot ../LDScores/${bedFile//.bed/.}${chr}.annot.gz \
    		--out ../LDScores/${bedFile//.bed/.}${chr} \
    		--print-snps ../hm3_no_MHC.list.txt &
    done;
done

process_id=$!
echo "Waiting for LDSC" &
wait $process_id

echo "Now running heritability analysis" &

### Now run heritability

for bedFile in $(ls | grep '.bed');
    do
    cd ..;
    cd sumstats;
    for sumstat in $(ls | grep 'sumstats.gz');
        mkdir results/${sumstats}
        python ../../../ldsc/ldsc.py \
                --h2 $sumstat \
                --ref-ld-chr ../LDScores/${bedFile//.bed/.},../1000G_Phase3_baselineLD/baselineLD. \
                --out $sumstat/results \
                --overlap-annot  \
                --frqfile-chr ../1000G_Phase3_frq/1000G.EUR.QC. \
                --w-ld-chr ../1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC.  \
                --print-coefficients &
    done;
done

exit 0
