# -*- coding: utf-8 -*-
"""
Created on Tue Jan 06 11:25:20 2015

@author: Sierra
"""

import ca_functions
import mgprofile
import pandas as pd
import numpy as np
import unittest
import random
import os

class test_comparative_analysis(unittest.TestCase): 
    
    def random_profile(self, normalization=True):
        # creates a random instance of the mgprofile with 50 samples in two classes
        # and 100 attributes
        attr_names = []
        for i in range(100):
            attr_names.append("attr" + str(i))
        samp_names = []
        for i in range(25):
            samp_names.append("a" + str(i))
        for i in range(25, 50):
            samp_names.append("b" + str(i))
        samples = []
        for i in range(len(samp_names)):    
            attributes = []
            for x in range(100):
                attributes.append(random.random())
            if normalization:
                attr_sum = sum(attributes)
                for j in range(100):
                    attributes[j] = attributes[j] / attr_sum
            samples.append((samp_names[i], attributes))
        
        f = open('temp.tab','w')
        f.write("attr\t")
        for i in range(len(samp_names)):
            if i != (len(samp_names) -1):
                f.write(samp_names[i] + '\t')
            else:
                f.write(samp_names[i])
        f.write("\n")
        
        for i in range(100):
            f.write(attr_names[i] + "\t")
            for j in range(len(samples)):
                attr_data = samples[j][1]
                if j != (len(samples) - 1):
                    f.write(str(attr_data[i]) + "\t")
                else:
                    f.write(str(attr_data[i]))
            f.write("\n")
                    
        f.close()
        
        f = open('temp2.tab', 'w')
        for s in samp_names:
            f.write(s + "\t")
            if s[0] == "a":
                f.write(str(1))
                f.write("\n")
            else:
                f.write(str(0))
                f.write("\n")
        f.close()
        
        p = mgprofile.metagenomic_profile()
        p.create_profile(title="test",collaborator="test",\
        year="2015", sequence_type="test", abundance_data="temp.tab",\
        metadata="temp2.tab")
        
        os.remove("temp.tab")
        os.remove("temp2.tab")        
        
        return p
        
    def skip_parameter_check(self):
        # checks if necessary files are available for the test_paramater_check
        # method
        testparams = ["goodparams1.sh", "goodparams2.sh", "goodparams3.sh", "goodparams4", \
        "badparams1.sh", "badparams2.sh"]
        try:
            for p in testparams:
                f = open(p, 'r')
                f.close()
        except IOError:
            return True
            
        return False
    
    def temp_file_setup(self):
        # sets up files for test_parameter_check, returns a list of the file names
        # created 
        files = []
        
        a_data = open("tempdata1.tab", 'w')
        a_data.write("\n")
        files.append(a_data.name)
        a_data.close()

        good_params_1 = ["ABUNDANCE_DATA=tempdata1.tab", "METADATA=tempdata1.tab"]
        gp_1 = open("gp1-temp.sh", 'w')
        for line in good_params_1:
            gp_1.write(line)
            gp_1.write("\n")
        files.append(gp_1.name)
        gp_1.close()
        
        good_params_2 = ["ABUNDANCE_DATA=tempdata1.tab", "METADATA=tempdata1.tab", \
        "ENRICHMENT=Yes", "ENRICHMENT_TEST=t-test", "NORMALIZATION=None", "OUTPUT_DIR=Current", \
        "PCA=yes", "PCA_TYPE=pca"]   
        gp_2 = open("gp2-temp.sh", 'w')
        for line in good_params_2:
            gp_2.write(line)
            gp_2.write("\n")
        files.append(gp_2.name)
        gp_2.close()
        
        good_params_3 = ["NAME=Sierra Anderson", "COLLABORATOR=No one in particular", \
        "SEQUENCE_TYPE=Illumina", "YEAR=2015", "ABUNDANCE_DATA=tempdata1.tab", \
        "METADATA=tempdata1.tab", "ENRICHMENT=Yes", "NORMALIZATION=None", "PCA=yes"]
        gp_3 = open("gp3-temp.sh", 'w')
        for line in good_params_3:
            gp_3.write(line)
            gp_3.write("\n")
        files.append(gp_3.name)
        gp_3.close()
        
        good_params_4 = ["METADATA=tempdata1.tab", "ENRICHMENT_TEST=t-test", \
        "COLLABORATOR=No one in particular", "PCA=yes", "SEQUENCE_TYPE=Illumina", \
        "NORMALIZATION=None", "ABUNDANCE_DATA=tempdata1.tab", "ENRICHMENT=Yes", \
        "OUTPUT_DIR=Current", "NAME=Sierra Anderson", "PCA_TYPE=pca", "YEAR=2015"]
        gp_4 = open("gp4-temp.sh", 'w')
        for line in good_params_4:
            gp_4.write(line)
            gp_4.write("\n")
        files.append(gp_4.name)
        gp_4.close()
                
        bad_params_1 = ["NAME=Sierra Anderson", "COLLABORATOR=No one in particular", \
        "SEQUENCE_TYPE=Illumina", "YEAR=2015", "ABUNDANCE_DATA=tempdata1.tab", \
        "ENRICHMENT=Yes", "ENRICHMENT_TEST=t-test", "NORMALIZATION=None", "OUTPUT_DIR=Current", \
        "PCA=yes", "PCA_TYPE=pca"]
        bp_1 = open("bp1-temp.sh", 'w')
        for line in bad_params_1:
            bp_1.write(line)
            bp_1.write("\n")
        files.append(bp_1.name)
        bp_1.close()
        
        bad_params_2 = ["NAME=Sierra Anderson", "COLLABORATOR=No one in particular", \
        "SEQUENCE_TYPE=Illumina", "YEAR=2015", "METADATA=tempdata1.tab", \
        "ENRICHMENT=Yes", "ENRICHMENT_TEST=t-test", "NORMALIZATION=None", "OUTPUT_DIR=Current", \
        "PCA=yes", "PCA_TYPE=pca"]
        bp_2 = open("bp2-temp.sh", 'w')
        for line in bad_params_2:
            bp_2.write(line)
            bp_2.write("\n")
        files.append(bp_2.name)
        bp_2.close()

        return files
        
        
    def temp_file_teardown(self, files):
        # removes the temporary files for test_parameter_check
        for fname in files:
            os.remove(fname)
        
    def test_parameter_check(self):
        # checks that the parameter function accepts good parameters 
        expected = ["name", "collaborator","sequence_type","year","abundance_data", \
        "metadata","enrichment","pca", "normalization", "output_dir"]
        
        optional_info = ["name","year","sequence_type","collaborator"]
        
        defaults = ["output_dir", "enrichment_test", "pca_type"]        
        
        files = self.temp_file_setup()
        
        # parameters with non-essential information missing
        good_parameters_1 = open('gp1-temp.sh','r')
        params = ca_functions.check_params(good_parameters_1)
        for p in expected:
            self.assertTrue(p in params.keys())
        good_parameters_1.close()
        
        good_parameters_2 = open('gp2-temp.sh','r')
        params = ca_functions.check_params(good_parameters_2)
        for p in expected:
            self.assertTrue(p in params)
        for p in optional_info:
            self.assertEqual(params[p], "n/a")
        good_parameters_2.close()
        
        # parameters that should generate default tests/output directory
        good_parameters_3 = open('gp3-temp.sh','r')
        params = ca_functions.check_params(good_parameters_3)
        for p in defaults:
            self.assertTrue(p in params)
        good_parameters_3.close()
        
        # parameters in a different order 
        good_parameters_4 = open('gp4-temp.sh', 'r')
        params = ca_functions.check_params(good_parameters_4)
        for p in expected:
            self.assertTrue(p in params)
        good_parameters_4.close()
        
        # parameters with essential information missing
        bad_parameters_1 = open('bp1-temp.sh','r')
        with self.assertRaises(NameError):
            ca_functions.check_params(bad_parameters_1)
        bad_parameters_1.close()
        
        bad_parameters_2 = open('bp2-temp.sh','r')
        with self.assertRaises(NameError):
            ca_functions.check_params(bad_parameters_2)
        bad_parameters_2.close()
        
        self.temp_file_teardown(files)        
        
    def test_sep_abundances(self):
        # checks that the function that seperates classes works correctly
        class_1 = []
        class_2 = []
        profile = self.random_profile(normalization=False)
        samples = list(profile.abundance_data.columns)
        for s in samples:
            if s[0] == 'a':
                class_1.append(s)
            else:
                class_2.append(s)
        # shuffle the metadata before passing it to the sep_abundances function
        df = profile.metadata
        profile.metadata = df.reindex(np.random.permutation(df.index))
        tup = ca_functions.sep_abundances(profile)
        test_1 = list(tup[0].index)
        test_2 = list(tup[1].index)
        self.assertEqual(sorted(class_1), sorted(test_2))
        self.assertEqual(sorted(class_2), sorted(test_1))
        
    def test_normalization(self):
        # tests that the relative normalization is working correctly 
        profile = self.random_profile(normalization=False)
        ca_functions.normalization(normalization_type="relative", mg_profile=profile, \
        output_dir="current", filename="normal-temp.tab")
        for c in profile.abundance_data.columns:            
            col = list(profile.abundance_data[c])            
            self.assertAlmostEqual(1, sum(col))
            
        profile = self.random_profile()
        for c in profile.abundance_data.columns:
            col = list(profile.abundance_data[c])
            self.assertAlmostEqual(1, sum(col))
            
        os.remove("normal-temp.tab")
        
    def test_enrichment(self):
        # first test the t-test
        
        profile = self.random_profile()            
        ca_functions.enrichment(enrichment_test="t-test", mg_profile=profile, \
        output_dir="current", filename="p-values-temp.tab")
        df = pd.DataFrame.from_csv('p-values-temp.tab', sep="\t")
        pvals = list(df["P-VALUE"])
        avg = sum(pvals) / len(pvals)
        self.assertAlmostEqual(0.5, avg, places=1)
        
        # test the rank-sums test
        avg_pvals = []
        for i in range(10):
            profile = self.random_profile()
            ca_functions.enrichment(enrichment_test="rank-sum", mg_profile=profile, \
            output_dir="current", filename="p-values-temp.tab")
            df = pd.DataFrame.from_csv('p-values-temp.tab', sep="\t")
            pvals = list(df["P-VALUE"])
            avg = sum(pvals) / len(pvals)
            avg_pvals.append(avg)
        
        self.assertAlmostEqual(0.5, sum(avg_pvals) / len(avg_pvals), places=1)
        
        os.remove("p-values-temp.tab")
        
    def test_pca(self):
        # tests the pca for correctness by doing calculations step by step 
        # and comparing to the pca generated by the sklearn.decomposition library
        profile = self.random_profile()        
        shape = profile.abundance_data.shape
        prf = np.array(profile.abundance_data)
        cov_matrix = np.cov(prf)
        eig_val, eig_vec = np.linalg.eig(cov_matrix)
        eig_pairs = [(eig_val[i], eig_vec[:,i]) for i in range(len(eig_val))]
        
        eig_pairs.sort()
        eig_pairs.reverse()
        
        wmatrix = np.hstack((eig_pairs[0][1].reshape(shape[0], 1), \
        eig_pairs[1][1].reshape(shape[0],1)))
        
        trans = wmatrix.T.dot(prf)
        
        test_pca = ca_functions.pca(pca_type="pca", mg_profile=profile, output_dir="current",\
        filename="pca-temp.png")
        
        pc_test_1 = sorted(test_pca[:, 0])
        pc_test_2 = sorted(test_pca[:, 1])
        pc_hand_1 = sorted(trans[0, :].real) # discards the complex part of a number,  
        # since it is discarded when plotting 
        pc_hand_2 = sorted(trans[1, :].real)
        
        for i in range(len(pc_test_1)):
            self.assertAlmostEquals(pc_test_1[i], pc_hand_1[i], places=1)
            self.assertAlmostEquals(pc_test_2[i], pc_hand_2[i], places=1)
        
        os.remove("pca-temp.png")
        
if __name__ == "__main__":
    unittest.main()