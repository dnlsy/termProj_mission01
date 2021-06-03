#!/bin/bash
SRAIN=$1
ENDPERIOD=$(($2-1))
OUTDIR=$3

for i in `eval echo {0..$ENDPERIOD}`; do
  srai=${1:0:3}`expr "${1:3:10}" "+" "$i"`
  fastq-dump -I --split-files ${srai} -O ${3}/
done

#ENDPERIOD=`expr "$2" "-" "1"`