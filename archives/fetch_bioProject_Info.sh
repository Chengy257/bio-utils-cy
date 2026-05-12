#!/bin/bash
#########################################################################
# File Name: /home/chengyu/soft/tools/get_BioProjectInfo.sh
# Author: ChengYu
# Description: 
# Created Time: Wed 10 May 2023 10:25:48 PM CST
#########################################################################
cat ${1} |while read id ;
do

#url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${id}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_md5,fastq_ftp,fastq_aspera,submitted_ftp,sra_md5,sra_ftp,sample_title&format=tsv&download=true&limit=0";
url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${id}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,library_name,nominal_length,library_selection,first_public,last_updated,experiment_title,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_ftp,sra_ftp,sample_alias&format=tsv&download=true&limit=0"

axel -n 10 $url -o SraRunInfo.$id ;

done
