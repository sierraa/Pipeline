# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 10:45:46 2015

@author: Sierra Anderson

Create an area plot from the abundance data for each class separately or all
classes together and save plot image to output directory. 

"""

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from random import random

MAX_DATA_POINTS = 50

def find_most_abundant(df, row_label):
    sample = df.loc[row_label]
    maximum = 0
    max_name = ''
    for k in sample.index:
        if sample[k] > maximum:
            max_name = k
            maximum = sample[k]
        
    return max_name

def sort_by_most_abundant(df):
    col_label = find_most_abundant(df, df.index[0])
    df.sort(columns=col_label, axis=0, inplace=True)

def generate_colors(n):
    # blue, yellow, red, green, magenta, sea green, orange, purple, lime green,
    # hot pink, cyan, dark red, dark blue, peach, gray, dark green, lavendar
    
    colors = ["#0000ff", "#ffff00", "#ff0000", "#00ff00", "#ff0066", "#99ff99", "#ff9900", 
         "#660066", "#99ff00", "#ff0099", "#99ffff", "#990000", "#000066", "#ff9966",
         "#c0c0c0", "#006600", "#cc99ff"]
    
    if len(colors) >= n:
        return colors
    
    for i in range(n - len(colors)):
        colors.append((random(), random(), random()))
    
    return colors

def sort_for_area_plot(df):
    srs = df.iloc[0]
    srs.sort(ascending=False)
    return list(srs.index)        
         
def plot(profile, filename="area_plot.png", colors=None):
    """Create area plot from abundance data and save figure.
    
    profile -- metagenomic profile instance containing abundance data.
    filename -- name of file to save area plot to (defaut="area_plot.png")
    colors -- colors to be used in the plot, default colors are randomly generated
    
    """
    
    if profile.abundance_data.shape[1] > MAX_DATA_POINTS:
        print("Too many data points to plot area plot.")
        return

    col_label = sort_for_area_plot(profile.abundance_data)

    profile.abundance_data.sort(columns=col_label, axis=0, inplace=True)

    if colors == None:
        colors = generate_colors(len(profile.abundance_data.columns))

    w = 0 # x coordinate to plot the new bar on
    
    prev = dict()
    plt.clf()
    plt.title("Area Plot")
    lgd_labels = dict() # stores information for the plot legend 
    
    for cls in profile.references.keys():
        df = profile.abundance_data.loc[profile.references[cls]]
        sort_by_most_abundant(df)
        
        # change order of columns so most abundant attribute is plotted first
        l = list(df.columns)[::-1]
        df = df[l]
        
        for sample in df.index:
            for i in range(len(df.columns)):
                attr = df.columns[i]                
                if i == 0:
                    prev[sample] = 0
                plt.bar(w, df.loc[sample, attr], linewidth=0, bottom=prev[sample], color=colors[i])
                prev[sample] += df.loc[sample, attr]
                if attr not in lgd_labels.keys():
                    lgd_labels[attr] = mpatches.Patch(color=colors[i], label=attr)
            w += 0.8

    ticks = list()
    ticks.append(0)
    running = 0
    for cls in profile.references.keys():
        running = running + len(profile.references[cls]) * 0.8
        ticks.append(running)
        plt.axvline(x=running, color='black')
        
    plt.xticks(ticks, list(profile.references.keys()))
    plt.xlim(0, len(profile.abundance_data.index) * 0.8)
    plt.ylim(0,1)
    plt.xlabel("Samples")
    plt.ylabel("Abundance")
    lgd = plt.legend(title="Attributes", handles=list(lgd_labels.values()), 
                     loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=8)
    plt.savefig(filename, bbox_extra_artists=(lgd,), bbox_inches='tight', 
                dpi=(400), figsize=(24, 24))