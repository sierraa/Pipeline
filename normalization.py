# -*- coding: utf-8 -*-
"""
Created on Tue Mar 03 10:43:09 2015

@author: Sierra Anderson

Normalizes the data according to a user specified technique. Updates the abundance
data for the metagenomic profile class and writes normalized data to file. 
"""

def normalize(profile, parameters):
    if parameters["normalization"] == "relative":
        relative_normalization(profile)
        profile.to_file_abundance_data("normalized.tab", parameters["output_dir"])
        
def relative_normalization(profile):
    for s in profile.abundance_data.index: #for each sample
        row_sum = sum(profile.abundance_data.ix[s])
        for i in range(len(profile.abundance_data.ix[s])):
            profile.abundance_data.ix[s, i] = profile.abundance_data.ix[s, i] / row_sum