# -*- coding: utf-8 -*-
"""
Created on Wed Mar 04 09:11:17 2015

@author: Sierra Anderson

Perform enrichment tests for this data and write results to file.
"""

from scipy import stats
#from qsturng import qsturng
import math

def enrichment_test(profile, parameters):
    if parameters["multiple_comparisons"][0] == "f":
        if parameters["enrichment_test"] == "ranksums":
            p, n = ranksums(profile)
        else:
            p, n = ttest(profile)
            
        to_file(p, n, parameters["output_dir"], parameters["enrichment_test"] + "_pvals.tab")
    else:
        #do multiple comparisons
        if parameters["pairwise"][0] == "t":
            corrected_pairwise(profile, parameters)
#        else:
#            tukey_hsd(profile, parameters)
 
def ranksums(profile):
    key_list = list(profile.references.keys())
    class_1 = profile.abundance_data.loc[profile.references[key_list[0]]]
    class_2 = profile.abundance_data.loc[profile.references[key_list[1]]]
    
    means = get_means(class_1, class_2)
    
    pvals = list()
    nans = list()
    
    for a in class_1.columns:
        directionality = "n/a"
        if means[a][0] > means[a][1]:
            directionality = key_list[0]
        elif means[a][1] > means[a][0]:
            directionality = key_list[1]
            p = stats.ranksums(class_1[a], class_2[a])[1]
        if not math.isnan(p):                
            pvals.append((p, a, directionality))
        else:
            nans.append((p, a))
    
    return pvals, nans
    
def ttest(profile):
    key_list = list(profile.references.keys())
    class_1 = profile.abundance_data.loc[profile.references[key_list[0]]]
    class_2 = profile.abundance_data.loc[profile.references[key_list[1]]]
    
    means = get_means(class_1, class_2)
    
    pvals = list()
    nans = list()
    
    for a in class_1.columns:
        directionality = "n/a"
        if means[a][0] > means[a][1]:
            directionality = key_list[0]
        elif means[a][1] > means[a][0]:
            directionality = key_list[1]
        p = stats.ttest_ind(class_1[a], class_2[a])[1]
        if not math.isnan(p):                
            pvals.append((p, a, directionality))
        else:
            nans.append((p, a))
    
    return pvals, nans
    
def corrected_pairwise(profile, parameters):
    """
    Perform Bonferroni one-step correction. 
    """
    key_list = profile.references.keys()
    m = len(key_list) - 1
    
    data = dict()
    
    for k in key_list:
        data[k] = profile.abundance_data.loc[profile.references[k]]
        
    for i in range(len(data.keys()) - 1):
        for j in range(i + 1, len(data.keys())):
            fname = str(data.keys()[i]) + "_vs_" + str(data.keys()[j]) + ".tab"
            if parameters["enrichment_test"] == "ranksums":
                p, n = ranksums(profile)
            else:
                p, n = ttest(profile)
            corrected_p = [x * m for x in p]
            to_file(corrected_p, n, parameters["output_dir"], fname)
"""
def tukey_hsd(profile, parameters, filename="tukey.tab", alpha=0.05):
    key_list = profile.references.keys()
    k = len(profile.references)
    
    sig = dict()
    
    for attr in profile.abundance_data.columns:
        data = dict()
        means = dict()        
        for key in key_list:
            l = profile.abundance_data.loc[profile.references[key], attr]
            means[key] = (sum(l)/len(l))
            data[key] = l
        mse, df = calc_mse(data.values())
        if sum(means.values()) == 0:
            continue
        for i in range(len(means) - 1):
            for j in range(i + 1, len(means)):
                s = math.sqrt(abs(mse*2))
                diff = abs(means[means.keys()[i]] - means[means.keys()[j]])
                
                omega = qsturng(1 - alpha, k, df)*(s / math.sqrt(len(data[data.keys()[i]]) + len(data[data.keys()[j]])))
                if diff >= omega:
                    sig[attr] = (omega/diff, (data.keys()[i], data.keys()[j]))
    
    if parameters["output_dir"] != "current":
        filename = parameters["output_dir"] + "\\" + filename
    
    f = open(filename, 'w')        
    f.write("NAME\tDIFFERENCES\tOMEGA/DIFF\n")

    for a in sig.keys():
        f.write(a)
        f.write("\t")
        f.write(str(sig[a][1]))
        f.write("\t")
        f.write(str(sig[a][0]))
        f.write("\n")
 """   
def calc_mse(data):
    """
    Calculate mean standard error
    Return MSE, degrees of freedom
    """
    cm = 0
    n = 0
    
    for l in data:
        cm += sum(l)
        n += len(l)
    
    cm = (cm)**2 / float(n)

    total_ss = 0    
    for l in data:
        for x in l:
            total_ss += x**2
    
    total_ss = total_ss - cm
    
    sst = 0
    for l in data:
        sst = sst + sum(l)**2/float(len(l))
        
    sst = sst - cm
    sse = total_ss - sst
    
    df = n - len(data) # degrees of freedom
    
    return sse / df, df
    
def get_means(df1, df2):
    means = dict()
    for a in df1.columns:
        class1_mean = stats.nanmean(df1[a])
        class2_mean = stats.nanmean(df2[a])
        means[a] = (class1_mean, class2_mean)
    
    return means
        
def to_file(pvals, nans, output_dir, filename):
    
    if (output_dir != "current"):
        filename = output_dir + "\\" + filename
        
    f = open(filename, "w")
    f.write("NAME\tP-VALUE\tDIRECTIONALITY\n")
    pvals.sort()
    
    for tup in pvals:
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

         