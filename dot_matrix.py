#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot dot matrix
just a test
"""
import os
import sys
import itertools
import numpy as np
import matplotlib.pyplot as plt

class dot_matrix:
    def __init__(self,seq1,seq2):
        self.vertical_seq = seq1
        self.horizonal_seq = seq2
        self.v_len = len(seq1)
        self.h_len = len(seq2)
        self.dot_matrix = np.zeros((self.v_len,self.h_len))
        self.x = []
        self.y = []
        for i,v in enumerate(seq1):
            for j,h in enumerate(seq2):
                if v == h:
                    self.dot_matrix[i,j] = 1
                    self.x.append(i)
                    self.y.append(j)
    def show_dot_matrix(self):
        dot = u'\u2022'
        print '  ' + ' '.join(self.horizonal_seq)
        for i in range(len(self.dot_matrix)):
            print self.vertical_seq[i],
            for j in range(len(self.dot_matrix[i])):
                if self.dot_matrix[i,j] == 1:
                    print '*',
                else:
                    print ' ',
            print
    def plot_dot_matrix(self):
        fig = plt.figure(figsize=(9,9),dpi=900)
        ax = fig.add_subplot(111,aspect='equal')
        ax.scatter(self.x,self.y,s=1)
        ax.set_xlim(0,self.h_len)
        ax.set_ylim(0,self.v_len)
        ax.yaxis.set_ticks(range(self.v_len)[::5])
        ax.yaxis.set_ticklabels(self.vertical_seq[::5])
        ax.xaxis.set_ticks(range(self.h_len)[::5])
        ax.xaxis.set_ticklabels(self.horizonal_seq[::5])
        fig.savefig('test')
        plt.close('all')

#@lt.log('test')
def test():
    a = 'RKTYYNVRILSDTWDQNRVIHYSGGDLIAVSLSDTWDQNVRNQDWTDRVIHYSGGSNRKIHLLDIIQVVKAILSDTWDQNRVIHYSGGPVEFRGHAGSVRALSDTWDQNRVIHYSGGLFLCE'
    c = 'RKKY'
    b = dot_matrix(a,a)
   # b.show_dot_matrix()
    b.plot_dot_matrix()

test()
