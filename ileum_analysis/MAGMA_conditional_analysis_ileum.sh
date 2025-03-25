#!/usr/bin/bash

#script to run MAGMA conditional analysis

#set variables 
folder="data"
file="update_EUR.CD.gwas_info03_filtered.assoc"
covar_file_spe="helm_batch1_13_ileum_conti_matrix_sig"
out_file_spe="helm_batch1_13_ileum_joint_spe_sig"

magma --gene-results ${folder}/${file}.step2.genes.raw --gene-covar ${folder}/${covar_file_spe}.txt --model joint-pairs direction=greater --out ${folder}/${out_file_spe}