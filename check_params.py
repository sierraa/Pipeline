# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 09:55:09 2015

@author: Sierra Anderson

Check parameters file for correctness. Exit if parameters file not found. Raise
NameError if crucial data (abundance data and metadata) is absent from parameters file.

"""

import os.path
import sys

def check_params(filename):
    """Check parameters passed to the script for correctness.
    Print message to user if parameters are incorrect.
    Raise NameError if filenames for data are not included or are not found.
    
    filename -- name of the file containing the parameters
    
    Return dictionary containing parameters.
    """
    
    if not os.path.isfile(filename):
        print("Parameters file not found.")
        sys.exit(0)
    
    filename = open(filename)
        
    parameters = dict()
    for line in filename:  
        line = line.split("=")
        line[0] = line[0].lower()
        parameters[line[0]] = line[1].rstrip().lower()   
    
    filename.close()
    
    # Check for essential parameters    
    
    if "metadata" not in parameters:
        raise NameError("Metadata file not found in parameters.")
    elif not os.path.isfile(parameters["metadata"]):
        raise NameError(str(parameters["metadata"]) + " not found.")
        
    if "abundance_data" not in parameters:
        raise NameError("Abundance data not found in parameters.")
    elif not os.path.isfile(parameters["abundance_data"]):
        raise NameError(str(parameters["abundance_data"]) + " not found.")
    
    # Check for non-essential parameters    
    
    if "name" not in parameters:
        parameters["name"]="n/a"
    
    if "year" not in parameters:
        parameters["year"]="n/a"
    
    if "sequence_type" not in parameters:
        parameters["sequence_type"]="n/a"
        
    if "collaborator" not in parameters:
        parameters["collaborator"]="n/a"    
        
    if "metadata_header" not in parameters:
        parameters["metadata_header"]="true"
    
    if "metadata_label" not in parameters:
        parameters["metadata_label"]="n/a"
    
    if "pca" not in parameters:
        parameters["pca"] = "n"
    
    if "pcoa" not in parameters:
        parameters["pcoa"] = "n"
        
    if "pcoa" in parameters and "dist_type" not in parameters:
        parameters["dist_type"] = "braycurtis"
    
    if "normalization" not in parameters:
        parameters["normalization"] = "none"
    
    if "output_dir" not in parameters:
        parameters["output_dir"] = "current"
    
    if "enrichment" not in parameters:
        parameters["enrichment"] = "f"
    elif "enrichment_test" not in parameters:
        parameters["enrichment_test"] = "ttest"
    
    if "pairwise" in parameters and "correction_type" not in parameters:
        parameters["correction_type"] = "bonferroni"
    
    if "multiple_comparisons" not in parameters:
        parameters["multiple_comparisons"] = "false"
    elif "pairwise" not in parameters:
        parameters["pairwise"] = "false"
        
    if "area_plot" not in parameters:
        parameters["area_plot"] = "false"
    if "plot_individual_classes" not in parameters:
        parameters["plot_individual_classes"] = "false"        
        
    if "to_html" not in parameters:
        parameters["to_html"] = "false"
    if "open_page" not in parameters:
        parameters["open_page"] = "false"
    
    return parameters