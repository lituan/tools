#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
use mds to plot high-dimentional data
"""
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA
import os
import sys
import scipy.spatial.distance as spd
import numpy as np

def read_pim(pim_f):
    with open(pim_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\n\r') for line in lines] # remove EOF
        lines = [line for line in lines if ':' in line] # remove comment lines
        header = [line.split()[1] for line in lines]
        matrix = [line.split()[2:] for line in lines]
        matrix = [map(lambda x: float(x)/100.0,row) for row in matrix]

        return header,matrix

def plot_mds_3d(labels,matrix,filename):

    scores = np.array(matrix)
    distances = [map(lambda x: 1.0-x,row) for row in scores]
    mds=manifold.MDS(n_components=3,max_iter=30000,eps=1e-12,dissimilarity='precomputed',n_jobs=-1)
    pos = mds.fit(distances)
    pos = model.embedding_
    print 'final stress: ' model.stress_
    # pos_d = euclidean_distances(pos)
    # print 'pos distances: ', pos_d

    fig=plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(pos[:,0],pos[:,1],pos[:,2],s=20,c='r')
    plt.savefig(filename+'.png',dpi=300)
    plt.show()
    plt.close('all')

def plot_mds(labels,matrix,filename):

    scores = np.array(matrix)
    distances = [map(lambda x: 1.0-x,row) for row in scores]
    mds=manifold.MDS(n_components=2,max_iter=30000,eps=1e-12,dissimilarity='precomputed',n_jobs=-1)
    pos = mds.fit(distances)
    pos = model.embedding_
    print 'final stress: ' model.stress_
    # pos_d = euclidean_distances(pos)
    # print 'pos distances: ', pos_d

    fig=plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(pos[:,0],pos[:,1],s=20,c='r')
    plt.savefig(filename+'.png',dpi=300)
    plt.show()
    plt.close('all')

def main():
    header,matrix = read_pim(sys.argv[-1])
    filename = os.path.splitext(os.path.split(sys.argv[-1])[1])[0]
    plot_mds(header,matrix,filename)


if __name__ == "__main__":
    main()
