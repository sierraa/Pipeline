# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 10:18:48 2015

@author: Sierra Anderson

Class instance representing metagenomic profile containing all data.
"""

import pandas as pd

class metagenomic_profile:
    
    def create_profile(self, parameters):
        """
        Initialize an instance of the metagenomic profile class given a 
        dictionary of parameters. 
        """
        # Initialize class attributes
        self.__class__.yr = parameters["year"]
        self.__class__.collab = parameters["collaborator"]
        self.__class__.seq = parameters["sequence_type"]
        self.__class__.name = parameters["name"]
        
        # Get abundance data and metadata (assumes tab-delimited)
        self.__class__.abundance_data = pd.DataFrame.from_csv(path=parameters["abundance_data"], sep='\t')
        
        if parameters["metadata_header"][0] != "f":
            self.__class__.metadata = pd.DataFrame.from_csv(path=parameters["metadata"], sep='\t')
        else:
            self.__class__.metadata = pd.DataFrame.from_csv(path=parameters["metadata"], sep='\t', header=None)
        
        # Only consider one label per instance of class - first label if none specified 
        if parameters["metadata_label"] != "n/a":
            self.__class__.metadata.sort(columns=parameters["metadata_label"], inplace=True)            
            self.__class__.metadata = self.__class__.metadata[parameters["metadata_label"]]
        else:
            self.__class__.metadata.sort(columns=self.__class__.metadata.columns[0], inplace=True)
            self.__class__.metadata = self.__class__.metadata[self.__class__.metadata.columns[0]]
        
        # keys = class labels, values = labels of class members, len(value) = class size
        self.__class__.references = dict()
        
        i = 1
        ref_type = self.__class__.metadata[0]
        self.__class__.references["class_" + str(ref_type)] = list()
        self.__class__.references["class_" + str(ref_type)].append(self.__class__.metadata.index[0])
        
        while i < len(self.__class__.metadata):
            if (self.__class__.metadata[i] == ref_type):
                self.__class__.references["class_" + str(ref_type)].append(self.__class__.metadata.index[i])
                i += 1
            else:
                ref_type = self.__class__.metadata[i]
                self.__class__.references["class_" + str(ref_type)] = list()
                self.__class__.references["class_" + str(ref_type)].append(self.__class__.metadata.index[i])
                i += 1
                
        self.__class__.num_of_classes = len(self.__class__.references)
        self.__class__.total_sample_count = 0
        
        for k in self.__class__.references.keys():
            self.__class__.total_sample_count += len(self.__class__.references[k])
        
        # assumes number of attributes > number of samples 
        # arranges the abundance data matrix to be samples x attributes
        if self.__class__.total_sample_count < len(self.__class__.abundance_data.index):
            self.__class__.abundance_data = self.__class__.abundance_data.T
    
    def to_file_abundance_data(self, filename, output_dir):
        if output_dir != "current":
            filename = output_dir + "//" + filename
        
        self.__class__.abundance_data.to_csv(path_or_buf=filename, sep="\t")
    
    def to_string(self):
        return "Name " + self.__class__.name + "\nCollaborator: " + self.__class__.collab \
        + "\nSequencing type: " + self.__class__.seq + "\n" + self.__class__.yr
        