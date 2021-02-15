 #!/bin/bash 

inpt=${1}
otp=${2}+"_illumina.vcf"
otp_two=${2}+"_other_tech.vcf"
grep 'Illumina' $inpt > $otp
 grep -v 'Illumina' $inpt > $otp_two