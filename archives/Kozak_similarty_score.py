#!/usr/bin/python
import numpy as np 
import sys
import argparse

# From https://github.com/Agleason1/TIS-Predictor/blob/main/Koazk_Similarity_Score_Algorithm.ipynb

# Kozak Consensus Scoring System

#0=A, 1=T, 2=G, 3=C, 4=N (Missing)
weights = np.array([
       [0.04210526, 0.        , 0.03157895, 0.05263158, 0.        ],
       [0.04210526, 0.05263158, 0.10526316, 0.0625    , 0.        ],
       [0.03157895, 0.04210526, 0.05263158, 0.07368421, 0.        ],
       [0.03157895, 0.01052632, 0.04210526, 0.05263158, 0.        ],
       [0.08421053, 0.07368421, 0.18947368, 0.10526316, 0.        ],
       [0.04210526, 0.05263158, 0.05263158, 0.08421053, 0.        ],
       [0.12631579, 0.0625    , 0.12631579, 0.21052632, 0.        ],
       [0.83157895, 0.12631579, 0.65263158, 0.16842105, 0.        ],
       [0.15789474, 0.06315789, 0.11578947, 0.2       , 0.        ],
       [0.21052632, 0.09473684, 0.31578947, 0.51578947, 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.24210526, 0.16666667, 0.53684211, 0.13684211, 0.        ],
       [0.15789474, 0.09473684, 0.09473684, 0.24210526, 0.        ],
       [0.05263158, 0.08421053, 0.14736842, 0.09473684, 0.        ],
       [0.07216495, 0.05263158, 0.10526316, 0.06315789, 0.        ],
       [0.        , 0.        , 0.        , 0.05263158, 0.        ],
       [0.05263158, 0.05263158, 0.10526316, 0.09473684, 0.        ],
       [0.04210526, 0.03157895, 0.05263158, 0.04210526, 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.04210526, 0.04210526, 0.08421053, 0.07368421, 0.        ],
       [0.0625    , 0.04210526, 0.09473684, 0.05263158, 0.        ]
])
#Below function scores using consensus kozak motif scores
def similarity_score(sequence):
    
    assert len(sequence)==23,'Sequence must be 23 bases long. Codon of interest must be centered, with 10 bases flanking both sides.'
           
 #We need consistency and flexibility:
    sequence = sequence.upper()
    for i in np.arange(len(sequence)):
        if sequence[i] =='U':
            sequence = sequence[0:i]+'T'+sequence[i+1:len(sequence)]
    
    numbers=[0]*len(sequence)    
    for k in np.arange(len(sequence)):
        if sequence[k]=='A':
            numbers[k] = 0
        elif sequence[k]=='T':
            numbers[k] = 1
        elif sequence[k]=='G':
            numbers[k] = 2
        elif sequence[k]=='C':
            numbers[k] = 3
        else:
            numbers[k]=4                
    global score 
    score = 0
    for k in np.arange(len(numbers)):
        score += weights[k][numbers[k]]           
    max_score = np.sum(weights.max(axis=1))       
    score = score/max_score    
    #Final scoring value: we take the maximum possible score 
    #calculated, and return our score divided by the maximum (to normalize from range 0 to 1)    
    return(score)

def getfas(IN):
    global fas 
    fas = {}
    with open(IN, 'r') as in_fa:
        for line in in_fa:
            line = line.strip()
            if line.startswith(">"):
                header = line.split(' ')[0]
                fas[header] = ''
            else:
                fas[header] += line
    in_fa.close()
    return fas

## main
parser = argparse.ArgumentParser()
parser.add_argument('--mod','-m',type=str,required=True,help="choose from: [fa] or [seq]; fasta from file or sequence from stdin")
parser.add_argument('--seq','-s',type=str,required=True,help="fasta file name or sequence")

args = parser.parse_args()

mod = args.mod
seq = args.seq

if mod == "fa":
    fa_file = str(seq)
    outfile = str(seq) + ".KSS"
    getfas(fa_file)
    with open(outfile,'w+') as out:
       for k,v in fas.items():
           similarity_score(v)
           print(k,score,file=out,sep="\t")
    out.close()
elif mod == "seq":
    seq = str(seq)
    if len(seq) == 23:
        print("-"*20,"\n","Koazk Similarity Score is:\n",similarity_score(seq))
    else:
        print("Please input a start codon seq in a length of 23")

else:
    print("Please choose mod from fa or seq!")

