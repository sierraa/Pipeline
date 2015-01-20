# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 10:15:10 2014

@author: Sierra Anderson
"""

import pandas as pd

class metagenomic_profile:
    
    def create_profile(self, title, collaborator, sequence_type, year, abundance_data=None, metadata=None):
        # Initializes the class attributes for this profile 
        self.__class__.yr = year
        self.__class__.collab = collaborator
        self.__class__.seq = sequence_type
        self.__class__.name = title
        
        # Creates a pandas DataFrame for both sample data and metadata.
        # Assumes files are tab delimited. 
        if (abundance_data != None and metadata != None):
            self.__class__.abundance_data = pd.DataFrame.from_csv(path=abundance_data.rstrip(), sep='\t')
            self.__class__.metadata = pd.DataFrame.from_csv(path=metadata.rstrip(), sep='\t', header=None)
        
    def to_string(self):
        return "Name " + self.__class__.name + "\nCollaborator: " + self.__class__.collab \
        + "\nSequencing type: " + self.__class__.seq + "\n" + self.__class__.yr
    
        