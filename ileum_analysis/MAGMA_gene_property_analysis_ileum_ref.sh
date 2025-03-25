#!/usr/bin/bash

#script to run MAGMA gene property analysis - ileum reference

#set variables 
folder="data"
file="update_EUR.CD.gwas_info03_filtered.assoc"
outfile_conti="helm_batch1_13_ileum_ref_conti_spe_2025"
covar_file_conti="data/ileum_ref_conti_specificity_matrix.txt"

magma --gene-results ${folder}/${file}.step2.genes.raw --model direction=greater --gene-covar ${covar_file_conti} --out ${folder}/${outfile_conti} 