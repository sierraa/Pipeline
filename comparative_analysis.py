#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 10:57:25 2014

@author: Sierra
"""

from sys import argv
import mgprofile
import ca_functions
import os
script, filename = argv

p = open(filename)
                
def main():
    
    params = ca_functions.check_params(p)
    
    mp = mgprofile.metagenomic_profile()
    mp.create_profile(title=params["name"], collaborator=params["collaborator"], \
    year=params["year"], sequence_type=params["sequence_type"], abundance_data=params["abundance_data"], \
    metadata=params["metadata"])
    
    # perform normalization
    if params["normalization"] != "n":
        ca_functions.normalization(normalization_type=params["normalization"], \
        mg_profile=mp, output_dir=params["output_dir"])
    
    # perform enrichment test
    if params["enrichment"][0] == "y":
        ca_functions.enrichment(enrichment_test=params["enrichment_test"], \
        mg_profile=mp, output_dir=params["output_dir"])
    
    # perform pca
    if params["pca"][0] == 'y':
        ca_functions.pca(params["pca_type"], mp, output_dir=params["output_dir"])
    
    if params["output_dir"] == "current":
        output_dir = os.getcwd()
    else:
        output_dir = params["output_dir"]
    
    print("Tests complete. Output saved at " + output_dir)
                
if __name__ == "__main__":
    main()