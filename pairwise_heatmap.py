#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example usage: python pairwise_heatmap.py --i /Users/me/pim.txt
"""
### hierarchical_clustering.py
#Copyright 2005-2012 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is furnished
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Changed by Li Tuan, imlituan@gmail.com
#the following code is for making heatmap for visualation of pairwise similarity

#################
### Imports an space-delimited expression matrix and produces a hierarchically clustered heatmap
### the matrix must have lables on the top and right,as the follow shows
#   lables    id1 id2
#   id1       1.0  0.2
#   id2       0.2  1.0
#################

import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import numpy
import string
import time
import sys, os
import getopt
import lt

################# Perform the hierarchical clustering #################

def heatmap(x, row_header, column_header, row_method,
            column_method, row_metric, column_metric,
            color_gradient, filename,max_labels=50,bottom_rotation=310):

    print "\nPerforming hiearchical clustering using %s for columns and %s for rows" % (column_metric,row_metric)

    """
    This below code is based in large part on the protype methods:
    http://old.nabble.com/How-to-plot-heatmap-with-matplotlib--td32534593.html
    http://stackoverflow.com/questions/7664826/how-to-get-flat-clustering-corresponding-to-color-clusters-in-the-dendrogram-cre
    x is an m by m ndarray
    """

    ### Define the color gradient to use based on the provided name
    n = len(x[0]); m = len(x)

    if color_gradient == 'blue_green_yellow':
        cmap=BlueGreenYellow()

    ### Scale the max and min colors so that 0 is white/black
    vmin=max(0, x.min())
    vmax=min(1.0, x.max())
    norm = mpl.colors.Normalize(vmin, vmax) ### Normalize values in x between vmin and vmax

    ### Scale the Matplotlib window size
    default_window_hight = 8.5
    default_window_width = 12
    fig = plt.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
    color_bar_w = 0.015 ### Sufficient size to show

    ## calculate positions for all elements
    # ax1, placement of dendrogram 1, on the left of the heatmap
    #if row_method != None: w1 =
    [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.22,0.2,0.6]   ### The second value controls the position of the matrix relative to the bottom of the view
    width_between_ax1_axr = 0.004
    height_between_ax1_axc = 0.004 ### distance between the top color bar axis and the matrix

    # axr, placement of row side colorbar
    [axr_x, axr_y, axr_w, axr_h] = [0.31,0.1,color_bar_w,0.6] ### second to last controls the width of the side color bar - 0.015 when showing
    axr_x = ax1_x + ax1_w + width_between_ax1_axr
    axr_y = ax1_y; axr_h = ax1_h
    width_between_axr_axm = 0.004

    # axc, placement of column side colorbar
    [axc_x, axc_y, axc_w, axc_h] = [0.4,0.63,0.5,color_bar_w] ### last one controls the hight of the top color bar - 0.015 when showing
    axc_x = axr_x + axr_w + width_between_axr_axm
    axc_y = ax1_y + ax1_h + height_between_ax1_axc
    height_between_axc_ax2 = 0.004

    # axm, placement of heatmap for the data matrix
    [axm_x, axm_y, axm_w, axm_h] = [0.4,0.9,2.5,0.5]
    axm_x = axr_x + axr_w + width_between_axr_axm
    axm_y = ax1_y; axm_h = ax1_h
    axm_w = axc_w

    # ax2, placement of dendrogram 2, on the top of the heatmap
    [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,0.15] ### last one controls hight of the dendrogram
    ax2_x = axr_x + axr_w + width_between_axr_axm
    ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
    ax2_w = axc_w

    # axcb - placement of the color legend
    [axcb_x, axcb_y, axcb_w, axcb_h] = [0.07,0.88,0.18,0.09]

    # transfer similarities to distance, 1-similarity
    def similar_to_dist(x):
        return [[1-xij for xij in xi] for xi in x]
    # Compute and plot top dendrogram
    if column_method != None:
        start_time = time.time()
        # d2 = dist.pdist(x.T)
        d2 = similar_to_dist(x)
        D2 = dist.squareform(d2)
        ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=True)
        Y2 = sch.linkage(D2, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
        Z2 = sch.dendrogram(Y2)
        ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
        ax2.set_xticks([]) ### Hides ticks
        ax2.set_yticks([])
        time_diff = str(round(time.time()-start_time,1))
        print 'Column clustering completed in %s seconds' % time_diff
    else:
        ind2 = ['NA']*len(column_header) ### Used for exporting the flat cluster data

    # Compute and plot left dendrogram.
    if row_method != None:
        start_time = time.time()
        # d1 = dist.pdist(x)
        d1 = similar_to_dist(x)
        D1 = dist.squareform(d1)  # full matrix
        ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True) # frame_on may be False
        Y1 = sch.linkage(D1, method=row_method, metric=row_metric) ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
        # Z1 = sch.dendrogram(Y1, orientation='right')
        Z1 = sch.dendrogram(Y1, orientation='left')
        ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
        ax1.set_xticks([]) ### Hides ticks
        ax1.set_yticks([])
        time_diff = str(round(time.time()-start_time,1))
        print 'Row clustering completed in %s seconds' % time_diff
    else:
        ind1 = ['NA']*len(row_header) ### Used for exporting the flat cluster data

    # Plot distance matrix.
    axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])  # axes for the data matrix
    xt = x
    if column_method != None:
        idx2 = Z2['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
        xt = xt[:,idx2]
        ind2 = ind2[idx2] ### reorder the flat cluster to match the order of the leaves the dendrogram
    if row_method != None:
        idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
        xt = xt[idx1,:]   # xt is transformed x
        ind1 = ind1[idx1] ### reorder the flat cluster to match the order of the leaves the dendrogram
    ### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
    # im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black
    im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, interpolation='nearest') ### norm=norm added to scale coloring of expression with zero = white or black
    # im = axm.imshow(xt, aspect='auto', origin='lower', cmap=cmap, interpolation='nearest') ### norm=norm added to scale coloring of expression with zero = white or black
    axm.set_xticks([]) ### Hides x-ticks
    axm.set_yticks([])

    # Add text
    new_row_header=[]
    new_column_header=[]
    for i in range(x.shape[0]):
        if row_method != None:
            if len(row_header)<max_labels: ### Don't visualize labels when more than max_labels,default is 50
                axm.text(x.shape[1]-0.5, i, '  '+row_header[idx1[i]])
            new_row_header.append(row_header[idx1[i]])
        else:
            if len(row_header)<max_labels:
                axm.text(x.shape[1]-0.5, i, '  '+row_header[i]) ### When not clustering rows
            new_row_header.append(row_header[i])
    for i in range(x.shape[1]):
        if column_method != None:
            if len(column_header)<max_labels:
                axm.text(i, -0.9, ' '+column_header[idx2[i]], rotation=bottom_rotation, verticalalignment="top") # rotation could also be degrees
            new_column_header.append(column_header[idx2[i]])
        else: ### When not clustering columns
            if len(column_header)<max_labels:
                axm.text(i, -0.9, ' '+column_header[i], rotation=bottom_rotation, verticalalignment="top")
            new_column_header.append(column_header[i])

    # Plot colside colors
    # axc --> axes for column side colorbar
    if column_method != None:
        axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
        cmap_c = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
        dc = numpy.array(ind2, dtype=int)
        dc.shape = (1,len(ind2))
        im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
        axc.set_xticks([]) ### Hides ticks
        axc.set_yticks([])

    # Plot rowside colors
    # axr --> axes for row side colorbar
    if row_method != None:
        axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
        dr = numpy.array(ind1, dtype=int)
        dr.shape = (len(ind1),1)
        #print ind1, len(ind1)
        cmap_r = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
        im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
        axr.set_xticks([]) ### Hides ticks
        axr.set_yticks([])

    # Plot color legend
    axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
    cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
    # cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='vertical')
    axcb.set_title("colorkey")

    if '/' in filename:
        dataset_name = string.split(filename,'/')[-1][:-4]
        root_dir = string.join(string.split(filename,'/')[:-1],'/')+'/'
    else:
        dataset_name = string.split(filename,'\\')[-1][:-4]
        root_dir = string.join(string.split(filename,'\\')[:-1],'\\')+'\\'
    filename = root_dir+dataset_name+"_"+'Clustering-%s-hierarchical_%s_%s.pdf' % (dataset_name,column_metric,row_metric)
    # cb.set_label("Differential Expression (log2 fold)")
    cb.set_label("sequence identity")
    exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2)

    ### Render the graphic
    if len(row_header)>50 or len(column_header)>50:
        plt.rcParams['font.size'] = 5
    else:
        plt.rcParams['font.size'] = 8

    plt.savefig(filename)
    print 'Exporting:',filename
    filename = filename[:-3]+'png'
    plt.savefig(filename, dpi=300) #,dpi=200
    #show the image on popup
    # plt.show()

def getColorRange(x):
    """ Determines the range of colors, centered at zero, for normalizing cmap """
    vmax=x.max()
    vmin=x.min()
    if vmax<0 and vmin<0: direction = 'negative'
    elif vmax>0 and vmin>0: direction = 'positive'
    else: direction = 'both'
    if direction == 'both':
        vmax = max([vmax,abs(vmin)])
        vmin = -1*vmax
        return vmax,vmin
    else:
        return vmax,vmin

################# Export the flat cluster data #################

def align_lis_lis(lis_lis):
    """align and  nested list to print a table"""
    lis_lis = [[str(l) for l in lis]
               for lis in lis_lis]  # trans every element to str
    #make all inner lists of the same length
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [lis + (inner_lis_max_len - len(lis)) * [''] for lis in lis_lis]
    #trans list, so that the elements of the same column are in one list
    lis_lis = [[lis[i] for lis in lis_lis] for i in range(inner_lis_max_len)]
    #make element in the same list have the same length
    aligned = []
    for lis in lis_lis:
        width = max([len(l) for l in lis])
        lis = [l + (width - len(l)) * ' ' for l in lis]
        aligned.append(lis)
    #trans list_list to the original list_list
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [[lis[i] for lis in aligned] for i in range(inner_lis_max_len)]
    return lis_lis

def exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2):
    """ Export the clustered results as a text file, only indicating the flat-clusters rather than the tree """

    # transform the matrix, so that origin is on left_bottom
    xt = xt[::-1]
    new_row_header = new_row_header[::-1]

    first_row = [' ','ID '] + new_column_header
    second_row = ['ID ','row_column_clusters_flat'] + map(str,ind2)
    rest_rows = [[new_row_header[i],str(ind1[i])]+map(str,row) for i,row in enumerate(xt)]
    rows = [first_row]+ [second_row]+ rest_rows
    rows = align_lis_lis(rows)
    filename = string.replace(filename,'.pdf','.txt')
    with open(filename,'w') as w_f:
        for row in rows:
            print >> w_f,'\t'.join(row)

    ### Export as CDT file
    first_row = ['UNIQID','NAME','GWEIGHT']+new_column_header
    second_row = ['EWEIGHT',' ',' '] + ['1']*len(new_column_header)
    rest_rows = [[new_row_header[i],str(ind1[i])]+map(str,row) for i,row in enumerate(xt)]
    rows = [first_row]+ [second_row]+ rest_rows
    rows = align_lis_lis(rows)
    filename = string.replace(filename,'.txt','.cdt')
    with open(filename,'w') as w_f:
        for row in rows:
            print >> w_f,'\t'.join(row)

    ### Export as labels file
    filename = string.replace(filename,'.txt','.id')
    with open(filename,'w') as w_f:
        for label in new_column_header:
            print >> w_f,'\t'.join(label)
################# Create Custom Color Gradients #################
#http://matplotlib.sourceforge.net/examples/pylab_examples/custom_cmap.html

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

################# General data import methods #################

def importData(filename):
    start_time = time.time()

    if '/' in filename:
        dataset_name = string.split(filename,'/')[-1][:-4] # unix format path
    else:
        dataset_name = string.split(filename,'\\')[-1][:-4] # windows format path

    with open(filename,'rU') as o_f:
        lines = o_f.readlines()
        lines = [line.strip('\n') for line in lines] # remove EOF
        column_header = lines[0].split()[1:]
        row_header = column_header
        matrix = [line.split()[1:] for line in lines[1:]]
        matrix = [map(float,row) for row in matrix]

    time_diff = str(round(time.time()-start_time,1))
    try:
        print '\n%d rows and %d columns imported for %s in %s seconds...' % (len(matrix),len(column_header),dataset_name,time_diff)
    except Exception:
        print 'No data in input file.'; force_error
    return numpy.array(matrix), column_header, row_header

def main():
    ################  Default Methods ################
    row_method = 'complete'
    column_method = 'complete'
    row_metric = 'euclidean'
    column_metric = 'euclidean'
    color_gradient = 'blue_green_yellow'

    """ Running with cosine or other distance metrics can often produce negative Z scores
        during clustering, so adjustments to the clustering may be required.

    see: http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    see: http://docs.scipy.org/doc/scipy/reference/spatial.distance.htm
    color_gradient = red_white_blue|red_black_sky|red_black_blue|red_black_green|yellow_black_blue|green_white_purple'
    """
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a tab-delimited input expression file in the command-line"
        print "Example: python pairwise_heatmap.py --i /Users/me/pim.txt"
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','row_header','column_method',
                                                    'row_metric','column_metric','color_gradient'])
        for opt, arg in options:
            if opt == '--i': filename=arg
            elif opt == '--row_header': row_header=arg
            elif opt == '--column_method': column_method=arg
            elif opt == '--row_metric': row_metric=arg
            elif opt == '--column_metric': column_metric=arg
            elif opt == '--color_gradient': color_gradient=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

    matrix, column_header, row_header = importData(filename)

    if len(matrix)>0:
        try:
            heatmap(matrix, row_header, column_header, row_method, column_method, row_metric, column_metric, color_gradient, filename,max_labels=50,bottom_rotation=310)
        except Exception:
            print 'Error using %s ... trying euclidean instead' % row_metric
            row_metric = 'euclidean'
            try:
                heatmap(matrix, row_header, column_header, row_method, column_method, row_metric, column_metric, color_gradient, filename,max_labels=50,bottom_rotation=310)
            except IOError:
                print 'Error with clustering encountered'

if __name__ == "__main__":
    main()
