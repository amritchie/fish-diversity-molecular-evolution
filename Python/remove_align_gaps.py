#! /Users/aritchie/anaconda3/bin/python3


## File: remove_align_gaps.py


## Remove columns with some proportion of gaps from alignment
## 
## Usage:
##
## python remove_align_gaps.py align_file.fasta 50 #Remove columns with 50% gaps
##
##

import sys
import os

from argparse import ArgumentParser
from Bio import Seq, SeqRecord, SeqIO, AlignIO

def make_parser():

    parser = ArgumentParser(description =
                            "Remove columns with the designated gap proportion from an alignment")
    parser.add_argument("align_file",
                        help = "The file containing the alignment to process.") 
    parser.add_argument("gap_prop",
                        type = float,
                        help = "Remove columns with this proportion of gaps",
                        default = 100.0)
    parser.add_argument("mode",
                        type = str,
                        help = "What type of gaps to remove. Answers are 'nucleotides', 'missing-codons', 'incomplete-codons'.",
                        choices = ["nucleotides", "missing-codons", "incomplete-codons"],
                        default = "nucleotides")
    return parser


if __name__ == "__main__":

    parser = make_parser()
    args = parser.parse_args()
    ext = os.path.splitext(args.align_file)[-1]

    formatstr = None
    if ext in [".fasta", ".fas"]:
        formatstr = "fasta"
    elif ext in [".phylip", ".phy"]:
        formatstr = "phylip-relaxed"
    elif ext in [".nex", ".nexus"]:
        formatstr = "nexus"

    align = AlignIO.read(args.align_file, formatstr)

    al_ungap = align[:,0:0]
    
    if (args.mode == "nucleotides"):
        for i in range(0, len(align[0])):
            print(100*sum([x in "-?" for x in align[:,i]])/len(align[:,i]))
            if args.gap_prop >= 100*sum([x in "-?" for x in align[:,i]])/len(align[:,i]):
                al_ungap += align[:,i:i+1]
    elif (args.mode == "incomplete-codons"):
        for i in range(0, len(align[0]), 3):
            print(100 * sum(any([x in "-?" for x in align[y,i:i+3]]) for y in range(0, len(align)))/len(align[:,i]))
            if (args.gap_prop >= 100 * sum(any([x in "-?" for x in align[y,i:i+3]]) for y in range(0, len(align)))/len(align[:,i])):
                al_ungap += align[:, i:i+3]
    elif (args.mode == "missing-codons"):
        for i in range(0, len(align[0]), 3):
            print(100 * sum(all([x in "-?" for x in align[y,i:i+3]]) for y in range(0, len(align)))/len(align[:,i]))
            if (args.gap_prop >= 100 * sum(all([x in "-?" for x in align[y,i:i+3]]) for y in range(0, len(align)))/len(align[:,i])):
                al_ungap += align[:, i:i+3]

    direc, fname  = os.path.split(args.align_file)
    out_name = os.path.join(direc, os.path.splitext(fname)[0]) + "_ungapped.fasta"
    
    AlignIO.write(al_ungap, out_name, format="fasta")

    print("Ungapped alignment written to ", out_name, ".")
