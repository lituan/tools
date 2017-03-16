#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cut sequences according to local alignment with a tomplate sequence
"""

import sys
import os
import operator
import itertools
from random import randint
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as spd
from collections import OrderedDict

from matplotlib.path import Path
from matplotlib import ticker
from matplotlib.patches import PathPatch
from svgpath2mpl import parse_path
from xml.dom import minidom

from wdsp import Wdsp


def align_seq(seq1, seq2):
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62
    gap_open = -10  # usual value
    gap_extend = -0.5  # usual value

    alns = pairwise2.align.localds(seq1, seq2, matrix, gap_open, gap_extend)

    print alns

    # seq1 = alns[0][0]
    # seq2 = alns[0][1]
    # identity = [1 for i, s in enumerate(seq1) if s == seq2[i]]
    # identity = 1.0 * len(identity) / len(seq1)

    # return identity



