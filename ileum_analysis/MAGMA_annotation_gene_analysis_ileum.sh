#!/usr/bin/bash

#script to run MAGMA annotation and Gene analysis

#set variables 

folder="data"
file="update_EUR.CD.gwas_info03_filtered.assoc"
snp_col="SNP"
p_col="P"
ncol="SampleSize"

#annotation
magma --annotate widow=35,10 \
--snp-loc ${folder}/snploc_${file} \
--gene-loc aux/NCBI37.3.gene.loc \
--out ${folder}/${file}.step1

#gene analysis
magma --bfile aux/g1000_eur --pval ${folder}/${file} use=${snp_col},${p_col} N=${ncol} --gene-annot ${folder}/${file}.step1.genes.annot --out ${folder}/${file}.step2