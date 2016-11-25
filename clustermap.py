#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
use seaborn clustermap to plot similarity heatmap
"""

import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as spd
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import os
import sys

def read_pim(pim_f):
    with open(pim_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\n\r') for line in lines] # remove EOF
        lines = [line for line in lines if ':' in line] # remove comment lines
        header = [line.split()[1] for line in lines]
        matrix = [line.split()[2:] for line in lines]

        return header,matrix



def BlueGreenYellow():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                        (0.25, 0.0, 0.0),
                        (0.50, 0.5, 0.5),
                        (1.0, 1.0, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (0.25, 0.75, 0.75),
                       (0.50, 1.0, 1.0),
                       (1.0, 1.0, 1.0)),

             'blue':  ((0.0, 1.0, 1.0),
                         (0.25, 1.0, 1.0),
                         (0.75, 0.0, 0.0),
                         (1.0, 0.0, 0.0))
            }
    ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
    ### modulate between blue and cyan using the last y var in the first green tuple
    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def heatmap(header,scores,file_name,method='average'):

    scores = [map(lambda x: float(x)/100.0,row) for row in matrix]
    scores = np.array(scores)
    distances = [map(lambda x: 1-x,row) for row in scores]
    linkage = sch.linkage(spd.squareform(distances),method=method)
    df = pd.DataFrame(scores,columns=header,index=header)
    # if more than 50 labels, then hide labels
    if len(df.columns) > 20:
        sns_plot = sns.clustermap(df,row_linkage=linkage,col_linkage=linkage,xticklabels='',yticklabels='')
    else:
        sns_plot = sns.clustermap(df,figsize=figsize,row_linkage=linkage,col_linkage=linkage,annot=True,fmt='3.2f')
        plt.setp(sns_plot.ax_heatmap.yaxis.get_majorticklabels(), rotation=20)
        plt.setp(sns_plot.ax_heatmap.xaxis.get_majorticklabels(), rotation=70)
    # plt.yticks(rotation=90)
    sns_plot.savefig(file_name+'.png')
    plt.close('all')



def main():
    pim_f = sys.argv[-1]
    file_path, file_name = os.path.split(pim_f)
    file_name, file_extention = os.path.splitext(file_name)
    header,scores= read_pim(pim_f)
    heatmap(header,scores,file_name)

if __name__ == "__main__":
    main()

