# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 10:12:43 2015

@author: Sierra Anderson

Script to perform a comparative analysis on metagenomic abundance data. Run
from the command line. 

usage -- python comparative_analysis.py parameters_file.sh

"""

from sys import argv
script, filename = argv
import os

from check_params import check_params
import mgprofile
import normalization
import pcoa
import enrichment
import area_plot
import to_html

def main():
    parameters = check_params(filename)
    
    # Create a metagenomic profile
    mp = mgprofile.metagenomic_profile()
    mp.create_profile(parameters)
    
    # Perform normalization
    if parameters["normalization"] != "none":
        normalization.normalize(mp, parameters)

    # Enrichment stuff 
    if parameters["enrichment"][0] == "t":
        enrichment.enrichment_test(mp, parameters)

    # PCA and PCOA
    if parameters["pca"][0] == "t":
        pcoa.plot(mp, parameters["output_dir"])

    if parameters["pcoa"][0] == "t":
        pcoa.plot(mp, parameters["output_dir"], dist_type=parameters["dist_type"])

    # Area plot
    if parameters["area_plot"][0] == "t":
        area_plot.plot(mp)

    # HTML page
    if parameters["to_html"][0] == "t":
        page_params = to_html.generate_params(parameters)
        to_html.write_page(page_params, parameters)
    print("Tests complete.")
    if parameters["output_dir"] == "current":
        parameters["output_dir"] = os.getcwd()
    print("Output saved at " + parameters["output_dir"])    

if __name__ == "__main__":
    main()