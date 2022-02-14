#! /Users/aritchie/anaconda3/bin/python3


## File: bigfish_trim_fillets.py

import sys
import os
import csv
import itertools as it

from argparse import ArgumentParser
from Bio import Seq, SeqRecord, SeqIO, AlignIO

def make_parser():

    parser = ArgumentParser(description =
                            "Remove columns with the designated gap proportion from an alignment")
    parser.add_argument("align_file",
                        help = "The file containing the alignment to process.")
    parser.add_argument("accspairs_file",
                        help = "The file containing the list of accessions, pairs and quartets")
    parser.add_argument("acc_text",
                        help = "The name of the accessions in the accspairs file.")
    parser.add_argument("min_bases",
                        type = int,
                        help = "Do not accept trimmed alignments with fewer than this number of bases.",
                        default = 100)
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
    accs = [x.id for x in align]

    accspairs = {}
    with open(args.accspairs_file) as fl:
        for row in csv.DictReader(fl):
            if row[args.acc_text] in accs:
                if row["pair_no"] in accspairs:
                    accspairs[row["pair_no"]].append(row[args.acc_text])
                else:
                    accspairs[row["pair_no"]] = [row[args.acc_text]]

    alnpairs = [sorted([x for x in range(0, len(align)) if align[x].id in accs]) for p, accs in accspairs.items()]

    print("Pairs:", alnpairs)
    
    missing = [0] * len(alnpairs)
    al_ungap = align[:,0:0]
    # go across columns by codon
    for i in range(0, len(align[0]), 3):
        # sites with every sister in every pair represented by at least one sequence
        accs_miss = [[any([y in "-?" for y in align[x,i:(i+3)]]) for x in p] for p in alnpairs]
        sisters_miss = [all(x) for x in it.chain.from_iterable([(y[0:len(y)//2], y[len(y)//2:len(y)]) for y in accs_miss])]
        if not any(sisters_miss):
            al_ungap += align[:, i:(i+3)]
        else:
            for n,k in enumerate(sisters_miss):
                if k:
                    missing[n//2] += 1

    # go across columns
    if len(al_ungap[0]) < args.min_bases:
        worst = missing.index(max(missing))
        print("!!!Alignment too small, dropping one pair:", alnpairs[worst])
        for i in range(0, len(align[0]), 3):
            al_ungap = align[:,0:0]
            alnpairs = alnpairs[:worst] + alnpairs[(worst+1):]
            accs_miss = [[any([y in "-?" for y in align[x,i:(i+3)]]) for x in p] for p in alnpairs]
            sisters_miss = [all(x) for x in it.chain.from_iterable([(y[0:len(y)//2], y[len(y)//2:len(y)]) for y in accs_miss])]
            if not any(sisters_miss):
                al_ungap += align[:, i:(i+3)]

    r_txt = ""
    if len(al_ungap[0]) < args.min_bases:
        print("!!!Alignment still too small, retaining all gapped positions.")
        al_ungap = align
        r_txt = "_ATTN"

    print("Ungapped row lengths:" + str(len(al_ungap[0])))

    direc, fname  = os.path.split(args.align_file)
    out_name = os.path.join(direc, os.path.splitext(fname)[0]) + r_txt + "_trim.fas"

    AlignIO.write(al_ungap, out_name, format="fasta")

    print("Ungapped alignment written to " + out_name + ".\n")
