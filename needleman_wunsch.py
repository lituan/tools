#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A needlman_wunsch algorithm for toy
"""
import os
import sys
import numpy as np
from random import randint

from Bio.SubsMat import MatrixInfo

SEQ_MATRIX = MatrixInfo.blosum62
UP = 1
LEFT = 2
DIAG = 0
INDEL = -7


def gap_penalty(i, j):
    return -3 * abs(i - j)


def score_hot(h1, h2):
    score = []
    for h11, h22 in zip(h1, h2):
        if (h11, h22) in SEQ_MATRIX.keys():
            score.append(SEQ_MATRIX[(h11, h22)])
        elif (h22, h11) in SEQ_MATRIX.keys():
            score.append(SEQ_MATRIX[(h22, h11)])
        else:
            score.append(-5)

    return sum(score)

def needleman_wunsch_matrix(hot1, hot2, score_fun):
    """
    calculate the matrix
    """

    m, n = len(hot1), len(hot2)
    score_matrix = np.zeros((m + 1, n + 1))
    pointer_matrix = np.zeros((m + 1, n + 1), dtype=int)

    for i in range(1, m + 1):
        score_matrix[i, 0] = INDEL * i + (-3) * i * (i + 1) * 0.5
    for j in range(1, n + 1):
        score_matrix[0, j] = INDEL * j + (-3) * j * (j + 1) * 0.5

    pointer_matrix[0, 1:] = LEFT
    pointer_matrix[1:, 0] = UP

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            up = score_matrix[i - 1, j] + INDEL + gap_penalty(i, j)
            left = score_matrix[i, j - 1] + INDEL + gap_penalty(i, j)
            diag = score_matrix[i - 1, j - 1] + \
                score_fun(hot1[i - 1], hot2[j - 1])

            max_score = max(up, left, diag)
            if max_score == diag:
                pointer_matrix[i, j] = DIAG
            elif max_score == up:
                pointer_matrix[i, j] = UP
            elif max_score == left:
                pointer_matrix[i, j] = LEFT

            score_matrix[i, j] = max_score
    return score_matrix, pointer_matrix


def needleman_wunsch_trace(hot1, hot2, score_matrix, pointer_matrix):
    align1 = []
    align2 = []
    m, n = len(hot1), len(hot2)

    curr = pointer_matrix[m, n]
    i, j = m, n
    while (i > 0 or j > 0):
        if curr == DIAG:
            align1.append(hot1[i - 1])
            align2.append(hot2[j - 1])
            i -= 1
            j -= 1
        elif curr == LEFT:
            align1.append('---')
            align2.append(hot2[j - 1])
            j -= 1
        elif curr == UP:
            align1.append(hot1[i - 1])
            align2.append('---')
            i -= 1
        curr = pointer_matrix[i, j]

    return align1[::-1], align2[::-1]


def needleman_wunsch(hot1, hot2):
    score_matrix, pointer_matrix = needleman_wunsch_matrix(
        hot1, hot2, score_hot)
    hot1, hot2 = needleman_wunsch_trace(
        hot1, hot2, score_matrix, pointer_matrix)

    hot11 = ''.join(hot1)
    hot22 = ''.join(hot2)
    identity = sum([1 for h1, h2 in zip(hot11, hot22) if h1 == h2])

    return identity,hot1,hot2

# d = ['SN*', 'GRC', 'CEN', 'KAA', 'CQY', 'LRA']
# e = ['QYG', 'SCD', 'GRR', 'CEN', 'KAT', 'CQY', 'LRG']
# print needleman_wunsch(d, e)
