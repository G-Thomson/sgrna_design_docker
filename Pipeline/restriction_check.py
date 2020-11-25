#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Restriction import RestrictionBatch, Restriction as re

# import argparse
# parser = argparse.ArgumentParser(description='Check sequences for our cloning restriction sites')
# parser.add_argument('input_seq',
#                        type=str,
#                        help='A input_seq string containing the sgRNA to look for restriction stites in')

# args = parser.parse_args()

def find_RE_sites(input_seq):
    rb = RestrictionBatch(['PacI', 'AscI', 'Esp3I', 'BsaI', 'BpiI'])

    new_seq = Seq(input_seq, IUPACAmbiguousDNA())
    srch = rb.search(new_seq, linear=True)

    count=0
    enz=[]
    for i in range(len(rb)):
        sites=len(list(srch.values())[i])
        count+=sites
        if sites > 0:
            enz.append(list(srch.keys())[i])

    return(count, enz)

# if __name__ == '__main__':
#     input_seq = args.input_seq
#     re = find_RE_sites(input_seq)
#     print(re)