#!/usr/bin/python
#########################################################################
# File Name: /home/chengyu/myscripts/batch_join.py
# Author: ChengYu
# Description: 
# Created Time: Fri 01 Mar 2024 11:09:04 AM CST
#########################################################################
import os
import pandas as pd
import csv 
import sys
filename_list = sys.argv[1]  
selected_columns= sys.argv[2]   ##  e.g.:   1,3  # for selected the No. 1 and 3 columns 
outfile = sys.argv[3]
# outfile = "merged_output.xls" 
dfs = []
selected_columns = [int(x)-1 for x in selected_columns.split(",")]
# print(selected_columns)
with open(filename_list, 'r') as file:
    filenames = file.readlines()
    for filename in filenames:
        df = pd.read_csv(filename.strip("\n"),sep="\t",quoting=csv.QUOTE_NONE, quotechar=None)
        name = os.path.basename(filename.strip("\n"))
        df = df.iloc[:,selected_columns]
        df.columns = ["ID",name]
        dfs.append(df)        
    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on= 'ID', how='outer')
        #print(merged_df)
    #merged_df.fillna(method="ffill")
file.close()
merged_df.to_csv(outfile,index=False,sep="\t",na_rep='NA')
