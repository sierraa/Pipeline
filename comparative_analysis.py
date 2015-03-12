# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 10:12:43 2015

@author: Sierra Anderson

Perform a comparative analysis on metagenomic abundance data.
"""

from sys import argv
script, filename = argv

from check_params import check_params
import mgprofile
import normalization
import pcoa
import enrichment
import area_plot

def main():
    parameters = check_params(filename)
    
    # Create a metagenomic profile
    mp = mgprofile.metagenomic_profile()
    mp.create_profile(parameters)
    
    # Perform normalization
    if parameters["normalization"] != "none":
        normalization.normalize(mp, parameters)
        
    # Enrichment stuff 
    if parameters["enrichment"][0] == "y":
        enrichment.enrichment_test(mp, parameters)
    
    # PCA and PCOA
    if parameters["pca"][0] == "y":
        pcoa.plot(mp, "euclidean", parameters["output_dir"])
    
    if parameters["pcoa"][0] == "y":
        pcoa.plot(mp, parameters["dist_type"], parameters["output_dir"])
        
    if parameters["area_plot"][0] == "y":
        area_plot.plot(mp, parameters)
    
    print("Tests complete.")    

if __name__ == "__main__":
    main()