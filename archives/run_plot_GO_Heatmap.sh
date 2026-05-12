#!/usr/bin/bash
GOcluster_IDs=$1
EXP=$2
Sample_Info=$3
#GMT="/share/data/reference/osa/GO/gProfiler/gprofiler_full_osativa.ENSG.gmt"
GMT="/home/chengyu/references/osa/GO/gProfiler/gprofiler_full_osativa.ENSG.gmt"
#Gene_Name_anno=$4
# input a GO ID list 
cd `pwd`
rm -rf ./tempExp.*
cat ${GOcluster_IDs}|while read GO
do
    # 1. get each GO pathways' gene ID
    GO_name=`cat ${GMT} |fgrep -w ${GO}|cut -f2|uniq|sed 's/ /_/g'`
    cat ${GMT} |fgrep -w ${GO}|cut -f3-|tr "\t" "\n"|sort -k1 > tempGENE.${GO_name}
    #join ${Gene_Name_anno} tempGENE.${GO_name}|cut -f2 > tempExp.${GO_name}.Gene_Name_anno
    # 2. subset gene exp
    head -1 ${EXP} > tempExp.${GO_name}
    cat ${EXP}|fgrep -w -f tempGENE.${GO_name} >> tempExp.${GO_name}
done
# 3. get all filepaths
ls `pwd`/tempExp.*|fgrep -v "Gene_Name_anno"|fgrep -v "tempExp.allfilenames" > tempExp.allfilenames
# 4. sample group 
cat ${EXP} |head -1|tr "\t" "\n" >  ${EXP}_headers
cat ${Sample_Info}|sed '1d'|fgrep -w -f ${EXP}_headers |cut -d, -f2 > temp.${GOcluster_IDs}.sampleHeader
# 5. plot heatmap
## plotHeatmap.R [filename] [sampleHeader] [] [output name]
 /opt/anaconda3/envs/R/bin/Rscript /home/chengyu/myscripts/plotHeatmap.R tempExp.allfilenames temp.${GOcluster_IDs}.sampleHeader "Heatmap_${GOcluster_IDs}"
# clean up
# rm -rf ./temp*
