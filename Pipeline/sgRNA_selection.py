#! /usr/bin/env python

import pandas as pd
from functools import reduce
import operator
import os.path

from Bio.SeqUtils import GC
from complexity import calc_complexity
from restriction_check import find_RE_sites

import argparse
parser = argparse.ArgumentParser(description='Select sgRNAs from candidate')
parser.add_argument('sgRNAs',
                       type=str,
                       help='CSV table of candidate sgRNAs')
parser.add_argument('orthogroup',
                       type=str,
                       help='String representing the orthogroup')

args = parser.parse_args()

def set_cover(universe, subsets):
    """Find a family of subsets that covers the universal set"""
    elements = set(e for s in subsets for e in subsets[s])

    # Check the subsets cover the universe
    if elements != universe:
        return None
    covered = set()
    cover = []
  
    while covered != elements:
        subset = max(subsets.keys(), key=lambda s: len(subsets[s] - covered))
        cover.append(subset)
        covered |= subsets[subset]

    return cover

def reformat_sgRNA(dat):
    ref_d = dat.groupby('sgRNA')[['gene']].apply(lambda g: g.values.tolist()).to_dict()

    for key, value in ref_d.items():
        ref_d[key] = set(reduce(operator.add, value))

    return ref_d

def p_table(sgRNAs):
    d_sub = dat[['sgRNA', 'gene']].loc[dat['sgRNA'].isin(sgRNAs)]

    d_sub['gene_number'] = d_sub.groupby('sgRNA').cumcount()+1
    d_sub['gene_number'] = d_sub['gene_number'].apply(lambda x: "gene_" + str(x))
    d_sub = d_sub.pivot(index='sgRNA', columns = 'gene_number').rename_axis(None, axis=0).droplevel(0, axis=1)
    d_sub.columns.name = "sgRNA"
    return d_sub

def get_sgRNA(dat, orthogroup):
    cutoffs = [300,500,800,1200,1500,2000]
    universe = set(dat.gene.unique())

    for c in cutoffs:
        sub_d = dat[dat.min_pos < c]
        set_d = reformat_sgRNA(sub_d)

        cover = set_cover(universe, set_d)

        if(cover == None):
            continue
        else:
            wide_table = p_table(cover) 
            wide_table.insert(loc=0, column='orthogroup', value=orthogroup)

            if os.path.exists('./wide_table.csv'):
                wide_table.to_csv('wide_table.csv', encoding='utf-8', mode = "a", header=False)
            else:
                wide_table.to_csv('wide_table.csv', encoding='utf-8')

            long_table = dat[['sgRNA', 'gene', 'pos_match', 'strand', 'GC']].loc[dat['sgRNA'].isin(cover)] 
            long_table.insert(loc=0, column='orthogroup', value=orthogroup)
            if os.path.exists('./long_table.csv'):
                long_table.to_csv('long_table.csv', encoding='utf-8', mode = "a", header=False, index=False)
            else:
                long_table.to_csv('long_table.csv', encoding='utf-8', index=False)

            break

if __name__ == '__main__':
    # Read data
    dat = pd.read_csv(args.sgRNAs)

    # Data filtering could be done here
    dat['GC'] = dat['sgRNA'].apply(lambda x: GC(x))
    dat['complexity'] = dat['sgRNA'].apply(lambda x: calc_complexity(x))
    dat['re_sites'] = dat['sgRNA'].apply(lambda x: find_RE_sites(x)[0])

    
    filtered_dat = dat.query('40 < GC < 65').query('complexity > 0.5').query('re_sites == 0')

    get_sgRNA(filtered_dat, args.orthogroup)
