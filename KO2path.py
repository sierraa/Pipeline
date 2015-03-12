# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 11:54:09 2015

@author: Sierra Anderson

Quick script to map KO to pathways
"""

from sys import argv
import pandas as pd 

script, initialdata, pathmap, outfile = argv

def main():
    """
    Given a tab-separated file mapping samples to KO abundances, a tab-separated
    file mapping pathways to KOs, script outputs a file mapping sample to pathway
    with unnormalized values naively calculated by adding up the abundance of each
    KO in a sample for each pathways set of KOs. 
    """
    initial = pd.DataFrame.from_csv(path=initialdata, sep='\t')
    if initial.shape[0] > initial.shape[1]:
        initial = initial.T
    pathm = pd.DataFrame.from_csv(path=pathmap, sep='\t')
    if pathm.shape[0] < pathm.shape[1]:
        pathm = pathm.T
    new_index = list()
    for i in pathm.index:
        if i in initial.columns:
            new_index.append(i)

    pathm = pathm.loc[new_index]
    pathm = pathm.T
    
    paths = dict() 
    
    for p in pathm.index: 
        paths[p] = list()
        for k in pathm.columns:
            if pathm.loc[p, k] == 1 and k in initial.columns:
                paths[p].append(k)
    
    df = pd.DataFrame(index=initial.index, columns=paths.keys())
    
    for s in initial.index:
        for p in pathm.index:
            d = sum(initial.loc[s, paths[p]])
            df[p][s] = d
    
    df.to_csv(path_or_buf=outfile, sep="\t")
    
if __name__ == "__main__":
    main()