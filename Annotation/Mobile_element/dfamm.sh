#!/bin/bash
path_dfam=${1}
vcf_fasta_format=${2}
hmm_model=${3}
otp_file=${4}


perl $path_dfam  --fastafile $vcf_fasta_format --hmmfile $hmm_model --dfam_outfile $otp_file