#!/usr/bin/env python3

# http://resources.qiagenbioinformatics.com/manuals/clccancerresearchworkbench/200/index.php?manual=How_sequence_complexity_is_calculated.html

import numpy as np 
# import argparse
# parser = argparse.ArgumentParser(description='')
# parser.add_argument('input_seq',
#                     type=str,
#                     help='A input_seq string containing the sgRNA to assess complexity')
# args = parser.parse_args()
  
def unique(list1): 
    x = np.array(list1) 
    return np.unique(x) 
      
def complexity(word, size): 
    l_seq = len(word)

    ngrams = []
    for i in range(l_seq-size+1):
        ngrams.append(word[i : i+size])

    if len(ngrams) < (4 ** size):
        c_score = len(unique(ngrams))/len(ngrams)
    else:
        c_score = len(unique(ngrams))/(4 ** size)

    return c_score

def calc_complexity(word, precision = 8):
    c = 1
    for i in range(1,precision): # seq length + 1
        c = c * complexity(word, i)
    
    return round(c, 3)

# if __name__ == '__main__':
#     input_seq = args.input_seq
#     re = calc_complexity(input_seq)
#     print(re)