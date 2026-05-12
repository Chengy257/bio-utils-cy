#!/bin/python
import sys
import argparse

dic = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}

def revcomp(seq):
    seq_revcomp = [dic[x] for x in str(seq)[::-1]]  ## reverse complement
    return seq_revcomp

def comp(seq):
    seq_comp = [dic[x] for x in str(seq)]  ## complement
    return seq_comp

def rev(seq):
    seq_rev = str(seq)[::-1]  ## reverse
    return seq_rev


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

def do(mod,seq):
    if mod == "revcomp": 		
        return ''.join(revcomp(seq))
    elif mod == "comp":
        return ''.join(comp(seq))
    elif mod == "rev":       
        return ''.join(rev(seq))

def main(mod,type,seq):
    if type == "fasta":
        fas = getfas(seq)
        for k,v in fas.items():
            print(k,do(mod,v),sep="\n")
    elif type == "stdin":
        seq = str(seq)
        print(do(mod,seq),sep="\n")
			

## main function
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help='Input')
parser.add_argument('-m','--mod',help='Mod. One of revcomp, rev, comp.',default='revcomp')
parser.add_argument('-t','--type',help='Input file type, fasta or stdin.',default='stdin')
args = parser.parse_args()
main(args.mod,args.type,args.input)



# if sys.argv[1] == "-fa":
#     IN = sys.argv[2]
    
#     getfas(IN)
#     for k,v in fas.items():
#         print(k,''.join(revcomp(v)),sep="\n")

# elif sys.argv[1] == "-seq":
#     seqs = str(sys.argv[2])
#     print("-"*30,"\n","reverse complement:\n",''.join(revcomp(seqs)))

# elif sys.argv[1] == "-h":
#     print("Usage: python",__file__,"-seq [sequence]\n       python",__file__,"-fa  [fasta file name] stdout")




##END
