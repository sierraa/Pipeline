# -*- coding: utf-8 -*-
"""
Created on Tue Mar 03 11:23:21 2015

@author: Sierra Anderson

Plot Principal Coordinate Analysis (PCoA). Save file to pcoa.png.    
Calculate using methods described in Numerical Ecology (Legendre 1998).     
A Principal Compoenent Analysis (PCA) is a PCoA performed with an 
euclidean distance matrix. This is the default analysis. 

"""
import pandas as pd
import numpy as np
import warnings
from sklearn.metrics.pairwise import pairwise_distances
from matplotlib import pyplot as plt

def plot(profile, output_dir, dist_type="euclidean", filename="pca.png"):
    """
    Generate PCoA or PCA plot for given metagenomic profile abundance data. 
    A PCA is a PCoA with a euclidean distance metric.
    
    profile -- instance of metagenomic profile class 
    output_dir -- directory for image to be saved in 
    dist_type -- distance metric to be used (default="euclidean" [PCA])
    filename -- name of file to be saved. If distance metric is not euclidean, 
    it will be saved as "pcoa_(dist_type).png" (default="pca.png")
    """
    warnings.simplefilter("ignore", np.ComplexWarning) # complex numbers must be cast to real in order to plot
    
    # blue, yellow, red, green, magenta, sea green, orange, purple, lime green,
    # hot pink, cyan, dark red, dark blue, peach, gray, dark green, lavendar
    colors = ["#0000ff", "#ffff00", "#ff0000", "#00ff00", "#ff0066", "#99ff99", "#ff9900", 
         "#660066", "#99ff00", "#ff0099", "#99ffff", "#990000", "#000066", "#ff9966",
         "#c0c0c0", "#006600", "#cc99ff"]
    markers = ['o', 'D', 'v', 'd', '<', 'h' '+', 's', '>', '|', 'p', 'H', '.', 'x', \
    '*', '^', ',', '_']
    
    abundances_to_class = dict()
    for k in profile.references.keys():
        abundances_to_class[k] = profile.abundance_data.loc[profile.references[k]]
        
    #now for each class i have the abundances associated w them
     
    df = pd.concat(list(abundances_to_class.values()))
    
    # effectively sorted now     
    
    dist_matrix = pairwise_distances(df, metric=dist_type)
    
    # 9.20 
    A_matrix = dist_matrix * dist_matrix / -2
        
    n = int(A_matrix.shape[0])  
    a_mean = np.mean(A_matrix)
        
    # 9.21
    ctr_matrix = [[0 for i in range(n)] for j in range(n)]        
        
    for i in range(n):
        for j in range(n):
            s = A_matrix[i][j] - np.mean(A_matrix[i][:]) - np.mean(A_matrix[:][j]) + a_mean
            ctr_matrix[i][j] = s
        
    eig_val, eig_vec = np.linalg.eig(ctr_matrix)

    eig_pairs = [(eig_val[i], eig_vec[:,i]) for i in range(len(eig_val))]

    eig_pairs.sort()
    eig_pairs.reverse()        
        
    PC1 = eig_pairs[0][1]      
    PC2 = eig_pairs[1][1]
    
    my_sum = 0
    for i in range(len(eig_val)):
        my_sum += eig_val[i]
    
    PC1_variance = eig_pairs[0][0]/my_sum * 100
    PC2_variance = eig_pairs[1][0]/my_sum * 100
    plt.clf()
    
    prev= 0
    marker_index = color_index = i = 0
    
    for k in abundances_to_class.keys():
        x = len(profile.references[k])
        plt.plot(PC1[prev:x+prev], PC2[prev:x+prev], markers[marker_index], markersize=5, \
        color=colors[color_index], label=k)
        marker_index = (marker_index + 1) % len(markers)
        color_index = (color_index + 1) % len(colors)          
        prev += x
        i += 1
    
    if dist_type != "euclidean":
        title = "PCoA"
    else:
        title = "PCA"
        
    plt.title(title)
    plt.xlabel('PC1 (' + str(round(PC1_variance, 2)) + '%)')
    plt.ylabel('PC2 (' + str(round(PC2_variance, 2)) + '%)')    
    ax = plt.subplot(111)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1))
        
    if dist_type != "euclidean":    
        filename =  "pcoa_" + dist_type + ".png"
    
    if output_dir != "current":
        filename = output_dir + "\\" + filename
    
    plt.savefig(filename, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=(200))    
    plt.clf()
    
    return PC1, PC2