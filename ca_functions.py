# -*- coding: utf-8 -*-
"""
Created on Thu Jan 08 10:48:19 2015

@author: Sierra
"""

from scipy import stats
import math
from sklearn.decomposition import PCA
import pandas as pd
from matplotlib import pyplot as plt

def check_params(p):
    # keep this updated as new parameters are added to the file
    
    parameters = dict()
    for line in p:  
        line = line.split("=")
        line[0] = line[0].lower()
        parameters[line[0]] = line[1].rstrip().lower()
        
    p.close()
    
    if "metadata" not in parameters:
        raise NameError("Metadata file not found in parameters.")
        
    if "abundance_data" not in parameters:
        raise NameError("Abundance data not found in parameters.")
        
    if "pca" not in parameters:
        parameters["pca"] = "n"
    elif "pca_type" not in parameters:
        parameters["pca_type"] = "pca"
    
    if "normalization" not in parameters:
        parameters["normalization"] = "none"
    
    if "output_dir" not in parameters:
        parameters["output_dir"] = "current"
    
    if "enrichment" not in parameters:
        parameters["enrichment"] = "n"
    elif "enrichment_test" not in parameters:
        parameters["enrichment_test"] = "default"
    
    if "name" not in parameters:
        parameters["name"]="n/a"
    
    if "year" not in parameters:
        parameters["year"]="n/a"
    
    if "sequence_type" not in parameters:
        parameters["sequence_type"]="n/a"
        
    if "collaborator" not in parameters:
        parameters["collaborator"]="n/a"
    
    return parameters
    
def sep_abundances(mg_profile):
    mg_profile.abundance_data = mg_profile.abundance_data.T
    mg_profile.abundance_data.sort(inplace=True)
    mg_profile.metadata.sort(inplace=True)
    
    bool_vector_1 = []
    for s in list(mg_profile.metadata.index):
        if mg_profile.metadata[1][s] == 1:
            bool_vector_1.append(True)
        else:
            bool_vector_1.append(False)
    
    abundance_exp = mg_profile.abundance_data[bool_vector_1]
    
    bool_vector_2 = []
    for t in bool_vector_1:
        bool_vector_2.append(not t)
    
    abundance_ctrl = mg_profile.abundance_data[bool_vector_2]
    
    mg_profile.abundance_data = mg_profile.abundance_data.T
    return (abundance_ctrl, abundance_exp)
    
def normalization(normalization_type, mg_profile, output_dir, filename="normalized.tab"):
    # Performs a relative normalization for the data
    mg_profile.abundance_data = mg_profile.abundance_data.T 
    if normalization_type == "relative":
        for s in list(mg_profile.abundance_data.index):
            row = list(mg_profile.abundance_data.ix[s])
            row_sum = sum(row)
            for i in range(len(row)):
                mg_profile.abundance_data.ix[s, i] = mg_profile.abundance_data.ix[s, i] / row_sum
        
    mg_profile.abundance_data = mg_profile.abundance_data.T
    
    # write the data to file 
    if output_dir != "current":
        filename = output_dir + "//" + filename
        
    mg_profile.abundance_data.to_csv(path_or_buf=filename, sep="\t")

def enrichment(enrichment_test, mg_profile, output_dir, filename="p-values.tab"):
    abundances = sep_abundances(mg_profile)
    abundance_ctrl = abundances[0]
    abundance_exp = abundances[1]    
    
    p_values = []
    nans = []         
    means = dict()
    
    for s in list(abundance_ctrl.columns):
        ctrl_mean = stats.nanmean(abundance_ctrl[s])
        exp_mean = stats.nanmean(abundance_exp[s])
        means[s] = (ctrl_mean, exp_mean)
    
    if enrichment_test == "t-test" or enrichment_test == "default":
        # do the t-test 
        for s in list(abundance_ctrl.columns):
            directionality = "n/a"
            if means[s][0] > means[s][1]:
                directionality = "control"
            elif means[s][1] > means[s][0]:
                directionality = "experimental"
            p = stats.ttest_ind(abundance_ctrl[s], abundance_exp[s])[1]
            if not math.isnan(p):                
                p_values.append((p, s, directionality))
            else:
                nans.append((p, s))
                
    elif enrichment_test == "rank-sum":
        # do the rank-sum test
        for s in list(abundance_ctrl.columns):
            directionality = "n/a"
            if means[s][0] > means[s][1]:
                directionality = "control"
            elif means[s][1] > means[s][0]:
                directionality = "experimental"
            p = stats.ranksums(abundance_ctrl[s], abundance_exp[s])[1]
            if not math.isnan(p):                
                p_values.append((p, s, directionality))
            else:
                nans.append((p, s))
    
    if (output_dir != "current"):
        filename = output_dir + "\\" + filename
        
    f = open(filename, "w")
    f.write("NAME\tP-VALUE\tDIRECTIONALITY\n")
    p_values.sort()
    
    for tup in p_values:
        f.write(str(tup[1]))
        f.write("\t")
        f.write("%.12f" % tup[0])
        f.write("\t")
        f.write(str(tup[2]))
        f.write("\n")
    
    for n in nans:
        f.write(str(n[1]))
        f.write("\t")
        f.write("n/a")
        f.write("\n")
    
    f.close()
    
def pca(pca_type, mg_profile, output_dir, filename="pca.png"):    
    if pca_type == "pcoa":
        None
        # do the pcoa (to-do)
    else: # do the pca
        abundances = sep_abundances(mg_profile)
        abundance_ctrl = abundances[0]
        abundance_exp = abundances[1]    
        ctrl_count = abundance_ctrl.shape[0]        
        df = pd.concat([abundance_ctrl, abundance_exp])
        
        mypca = PCA(n_components=2)
        sklearnpca = mypca.fit_transform(df)
        plt.plot(sklearnpca[0:ctrl_count,0], sklearnpca[0:ctrl_count,1], 'o', markersize=7, \
        color='blue', alpha=0.5, label="control")
        
        plt.plot(sklearnpca[ctrl_count:,0], sklearnpca[ctrl_count:,1], '^', markersize=7, \
        color='red', alpha=0.5, label="experimental")
        
        plt.legend()
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title('PCA')
        
        if output_dir != "current":
            filename = output_dir + "\\" + filename
        
        plt.savefig(filename)
        
        return sklearnpca