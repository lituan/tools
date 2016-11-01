#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot ratio map according to a series of cutoff
"""
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

"""
input ratio matrix:

cutoff      total label1 label2 label3 ...
0           1.0   0.2     0.3   0.5
0.1         0.9   0.18    0.28  0.54
...      ......
0.9         0.1   0.16    0.26  0.58
"""


#"cutoff",  "Total",  "Animalia",  "Plantae",  "Fungi",  "Protista",  "Archaea",  "Bacteria",  "Virus",  "Other"
input_ratio_matrix = [
    [0.0, 1.0000, 0.3470, 0.1209, 0.3272, 0.1389, 0.0010, 0.0640, 0.0009, 0.0000],
    [0.1, 0.9944, 0.3466, 0.1212, 0.3274, 0.1386, 0.0010, 0.0643, 0.0009, 0.0000],
    [0.2, 0.3619, 0.2829, 0.1150, 0.3306, 0.1361, 0.0022, 0.1322, 0.0009, 0.0001],
    [0.3, 0.0427, 0.1096, 0.0179, 0.2639, 0.0841, 0.0047, 0.5195, 0.0000, 0.0004],
    [0.4, 0.0213, 0.0043, 0.0057, 0.2996, 0.0970, 0.0057, 0.5876, 0.0000, 0.0000],
    [0.5, 0.0129, 0.0000, 0.0036, 0.3405, 0.1234, 0.0071, 0.5255, 0.0000, 0.0000],
    [0.6, 0.0076, 0.0000, 0.0020, 0.3327, 0.1443, 0.0100, 0.5110, 0.0000, 0.0000],
    [0.7, 0.0047, 0.0000, 0.0000, 0.3574, 0.1311, 0.0098, 0.5016, 0.0000, 0.0000],
    [0.8, 0.0023, 0.0000, 0.0000, 0.4067, 0.1067, 0.0000, 0.4867, 0.0000, 0.0000],
    [0.9, 0.0005, 0.0000, 0.0000, 0.5000, 0.0313, 0.0000, 0.4688, 0.0000, 0.0000]]


# cutoff    ratio    Animalia    Plantae    Fungi    Protista
# Archaea_Bacteria
input_ratio_matrix = [
    [0.0,    1.0000,    0.3474,    0.1211,    0.3275,    0.1391,   0.0650, ],
    [0.2,    0.3619,    0.2832,    0.1151,    0.3310,   0.1363,    0.1345, ],
    [0.3,    0.0427,    0.1096,    0.0179,    0.2640,   0.0842,    0.5244, ],
    [0.4,    0.0213,    0.0043,    0.0057,    0.2996,   0.0970,    0.5934, ],
    [0.5,    0.0129,    0.0000,    0.0036,    0.3405,   0.1234,    0.5326, ],
    [0.6,    0.0076,    0.0000,    0.0020,    0.3327,   0.1443,    0.5210, ],
    [0.7,    0.0047,    0.0000,    0.0000,    0.3574,   0.1311,    0.5115, ]]


# set color for the labels
LABELS = ["Animalia",   "Plantae",   "Fungi",
          "Protista",   "Archaea_Bacteria"]
COLORS = [[0.98, 0.8, 0.68],
          [0.44, 0.68, 0.27],
          [0.62, 0.81, 0.98],
          [1.0, 1.0, 0.0],
          [0.97, 0.02, 0.91]]

EDGE_COLOR = [0.32, 0.51, 0.69]

# set the coordiates of the inverted triangle
# the bottom of the triangle is parallel to the x-axis
HEIGHT = 1000.0
BASE = 1000.0
LEFT = (0, HEIGHT)
BOTTOM = (BASE / 2.0, 0)
RIGHT = (BASE, HEIGHT)

X_LIM = [0 - BASE * 0.1, BASE + BASE * 0.1]
Y_LIM = [0 - HEIGHT * 0.1, HEIGHT + HEIGHT * 0.1]


def dist(p1, p2):
    return np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# calculate coordiates of the respective ladder in the inversed triangle
# according to the ratio


def get_ladder_coordiate(left, bottom, right, input_ratio_matrix):

    tri_bottom_len = right[0] - left[0]
    ratios = [1 - r[1] for r in input_ratio_matrix[1:]]

    ladder_coordiates = [(left, right)]
    for ratio in ratios:
        ladder_bottom_len = np.sqrt(1.0 - ratio) * tri_bottom_len
        ladder_bottom_y = bottom[
            1] - (1.0 * ladder_bottom_len / tri_bottom_len) * (bottom[1] - right[1])
        ladder_bottom_left = (
            left[0] + 0.5 * (tri_bottom_len - ladder_bottom_len), ladder_bottom_y)
        ladder_bottom_right = (
            right[0] - 0.5 * (tri_bottom_len - ladder_bottom_len), ladder_bottom_y)
        ladder_coordiates.append((ladder_bottom_left, ladder_bottom_right))
    ladder_coordiates.append((bottom, bottom))

    ladder_coordiates = [ladder_coordiates[i:i + 2]
                         for i in range(len(ratios) + 1)]
    ladder_coordiates = [(l[0][0], l[0][1], l[1][0], l[1][1])
                         for l in ladder_coordiates]

    return ladder_coordiates

# calculate coordiates of the respective polytriangle in the inversed ladder according
# to the ratio


def get_ratio_coordiate((ladder_up_left, ladder_up_right, ladder_bottom_left, ladder_bottom_right), ratio):

    # ladder
    if ladder_bottom_left != ladder_bottom_right:
        ladder_up_len = ladder_up_right[0] - ladder_up_left[0]
        ladder_bottom_len = ladder_bottom_right[0] - ladder_bottom_left[0]

        cutoff_ratio_left = 0.5 * \
            (ladder_up_len - ladder_bottom_len) / \
            (ladder_up_len + ladder_bottom_len)
        cutoff_ration_right = 1 - cutoff_ratio_left

        if ratio == 0:
            ratio_up = ratio_bottom = ladder_up_left

        elif ratio < cutoff_ratio_left:
            ratio_up_len = np.sqrt(ratio / cutoff_ratio_left) * \
                0.5 * (ladder_up_len - ladder_bottom_len)
            ratio_up = (ladder_up_left[0] + ratio_up_len, ladder_up_left[1])
            ratio_bottom_y = ladder_up_left[1] + (ladder_bottom_left[1] - ladder_up_left[
                                                  1]) * (2 * ratio_up_len / (ladder_up_len - ladder_bottom_len))
            ratio_bottom = (ladder_up_left[0] + ratio_up_len, ratio_bottom_y)

        elif ratio == cutoff_ratio_left:
            ratio_up_len = 0.5 * (ladder_up_len - ladder_bottom_len)
            ratio_up = (ladder_up_left[0] + ratio_up_len, ladder_up_left[1])
            ratio_bottom = ladder_bottom_left

        elif ratio < cutoff_ration_right:
            ratio_up_len = 0.25 * \
                (2 * ratio * (ladder_up_len + ladder_bottom_len) +
                 ladder_up_len - ladder_bottom_len)
            ratio_up = (ladder_up_left[0] + ratio_up_len, ladder_up_left[1])
            ratio_bottom = (ladder_up_left[0] +
                            ratio_up_len, ladder_bottom_left[1])

        elif ratio == cutoff_ration_right:
            ratio_up_len = 0.5 * (ladder_up_len + ladder_bottom_len)
            ratio_up = (ladder_up_left[0] + ratio_up_len, ladder_up_left[1])
            ratio_bottom = ladder_bottom_right

        elif ratio > cutoff_ration_right:
            ratio_after = 1 - ratio
            ratio_up_len = ladder_up_len - \
                np.sqrt(ratio_after / cutoff_ratio_left) * \
                0.5 * (ladder_up_len - ladder_bottom_len)
            ratio_up = (ladder_up_left[0] + ratio_up_len, ladder_up_left[1])
            ratio_bottom_y = ladder_up_left[1] + (ladder_bottom_left[1] - ladder_up_left[1]) * (
                2 * (ladder_up_len - ratio_up_len) / (ladder_up_len - ladder_bottom_len))
            ratio_bottom = (ladder_up_left[0] + ratio_up_len, ratio_bottom_y)

        elif ratio == 1.0:
            ratio_up = ratio_bottom = ladder_up_right

    # triangle
    elif ladder_bottom_left == ladder_bottom_right:

        right = ladder_up_right
        left = ladder_up_left
        bottom = ladder_bottom_left

        tri_bottom_len = right[0] - left[0]
        if ratio == 0:
            ratio_up = ratio_bottom = left
        elif ratio < 0.5:
            ratio_up_x = left[0] + 0.5 * tri_bottom_len * np.sqrt(2.0*ratio)
            ratio_up = (ratio_up_x, left[1])
            ratio_bottom_y = np.sqrt(2*ratio) * (bottom[1] - left[1]) + left[1]
            ratio_bottom = (ratio_up_x, ratio_bottom_y)
        elif ratio == 0.5:
            ratio_up = (left[0] + 0.5 * tri_bottom_len, left[1])
            ratio_bottom = bottom
        elif ratio < 1.0:
            ratio_up_x = right[0] - 0.5 * tri_bottom_len * np.sqrt(2.0 - 2.0*ratio)
            ratio_up = (ratio_up_x, right[1])
            ratio_bottom_y = np.sqrt(2.0 - 2.0*ratio) * \
                (bottom[1] - left[1]) + left[1]
            ratio_bottom = (ratio_up_x, ratio_bottom_y)
        elif ratio == 1.0:
            ratio_up = ratio_bottom = right

    return (ratio_up, ratio_bottom)


# get polygons in each ladder according to ratio of numbers under each label
def get_polygons(input_ratio_matrix, ladder_coordiates):

    label_cum_ratios = [np.cumsum(r[2:]) for r in input_ratio_matrix]
    label_cum_ratios = [[0 for r in ratios if r <= 0] + [r for r in ratios if 0 <
                                                         r < 1.0] + [1.0 for r in ratios if r >= 1.0] for ratios in label_cum_ratios]

    polygons = []
    for ratio, ladder_coordiate in zip(label_cum_ratios, ladder_coordiates):
        llu, lru, llb, lrb = ladder_coordiate
        ratio_coordiates = [get_ratio_coordiate(
            ladder_coordiate, r) for r in ratio]

        ladder_polygons = []

       # ladder
        if llb[0] != lrb[0]:

            up1, bottom1 = ratio_coordiates[0]
            if up1[0] == llu[0]:
                polygon = (llu, llu, llu)
                ladder_polygons.append(polygon)
            elif up1[0] <= llb[0]:
                polygon = (llu, bottom1, up1)
                ladder_polygons.append(polygon)
            elif up1[0] <= lrb[0]:
                polygon = (llu, llb, bottom1, up1)
                ladder_polygons.append(polygon)
            elif up1[0] < lru[0]:
                polygon = (llu, llb, lrb, bottom1, up1)
                ladder_polygons.append(polygon)
            elif up1[0] == lru[0]:
                polygon = (llu, llb, lrb, lru)
                ladder_polygons.append(polygon)

            for (up1, bottom1), (up2, bottom2) in zip(ratio_coordiates[:-1], ratio_coordiates[1:]):

                # empty polygon
                if up1[0] == up2[0]:
                    polygon = (up1, bottom1, bottom2, up2)
                    ladder_polygons.append(polygon)

                elif up1[0] == llu[0] and up1[0] < up2[0] <= llb[0]:
                    polygon = (llu, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif up1[0] == llu[0] and llb[0] < up2[0] <= lrb[0]:
                    polygon = (up1, llb, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif up1[0] == llu[0] and lrb[0] < up2[0] < lru[0]:
                    polygon = (up1, llb, lrb, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif up1[0] == llu[0] and up2[0] == lru[0]:
                    polygon = (llu, llb, lrb, lru)
                    ladder_polygons.append(polygon)

                elif llu[0] < up1[0] < llb[0] and up1[0] < up2[0] <= llb[0]:
                    polygon = (up1, bottom1, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif llu[0] < up1[0] < llb[0] and llb[0] < up2[0] <= lrb[0]:
                    polygon = (up1, bottom1, llb, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif llu[0] < up1[0] < llb[0] and lrb[0] < up2[0] < lru[0]:
                    polygon = (up1, bottom1, llb, lrb, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif llu[0] < up1[0] < llb[0] and up2[0] == lru[0]:
                    polygon = (up1, bottom1, lrb, lru)
                    ladder_polygons.append(polygon)

                elif up1[0] == llb[0] and up2[0] <= lrb[0]:
                    polygon = (up1, bottom1, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif up1[0] == llb[0] and lrb[0] < up2[0] < lru[0]:
                    polygon = (up1, bottom1, lrb, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif up1[0] == llb[0] and up2[0] == lru[0]:
                    polygon = (up1, bottom1, lrb, up2)
                    ladder_polygons.append(polygon)

                elif llb[0] < up1[0] < lrb[0] and up2[0] <= lrb[0]:
                    polygon = (up1, bottom1, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif llb[0] < up1[0] < lrb[0] and lrb[0] < up2[0] < lru[0]:
                    polygon = (up1, bottom1, lrb, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif llb[0] < up1[0] < lrb[0] and up2[0] == lru[0]:
                    polygon = (up1, bottom1, lrb, lru)
                    ladder_polygons.append(polygon)

                elif up1[0] == lrb[0] and up2[0] < lru[0]:
                    polygon = (up1, bottom1, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif up1[0] == lrb[0] and up2[0] == lru[0]:
                    polygon = (up1, bottom1, lru)
                    ladder_polygons.append(polygon)

                elif lrb[0] < up1[0] < lru[0] and up2[0] < lru[0]:
                    polygon = (up1, bottom1, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif lrb[0] < up1[0] < lru[0] and up2[0] == lru[0]:
                    polygon = (up1, bottom1, lru)
                    ladder_polygons.append(polygon)

        # triangle
        if llb[0] == lrb[0]:

            left, right = llu, lru
            bottom = llb

            up1, bottom1 = ratio_coordiates[0]
            if up1[0] == left[0]:
                polygon = (left, left, left)
                ladder_polygons.append(polygon)
            elif up1[0] <= bottom[0]:
                polygon = (left, bottom1, up1)
                ladder_polygons.append(polygon)
            elif up1[0] < right[0]:
                polygon = (left, bottom, bottom1, up1)
                ladder_polygons.append(polygon)
            elif up1[0] < lru[0]:
                polygon = (llu, llb, lrb, bottom1, up1)
                ladder_polygons.append(polygon)
            elif up1[0] == lru[0]:
                polygon = (llu, llb, lrb, lru)
                ladder_polygons.append(polygon)

            for (up1, bottom1), (up2, bottom2) in zip(ratio_coordiates[:-1], ratio_coordiates[1:]):

                # empty polygon
                if up1[0] == up2[0]:
                    polygon = (up1, bottom1, bottom2, up2)
                    ladder_polygons.append(polygon)

                elif up1[0] == left[0] and up2[0] <= bottom[0]:
                    polygon = (left, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif up1[0] == left[0] and bottom[0] < up2[0] < right[0]:
                    polygon = (left, bottom, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif up1[0] == left[0] and up2[0] == right[0]:
                    polygon = (left, bottom, up2)
                    ladder_polygons.append(polygon)
                elif bottom[0] > up1[0] > left[0] and up2[0] <= bottom[0]:
                    polygon = (up1, bottom1, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif bottom[0] > up1[0] > left[0] and bottom[0] < up2[0] < right[0]:
                    polygon = (up1, bottom1, bottom, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif bottom[0] > up1[0] > left[0] and up2[0] == right[0]:
                    polygon = (up1, bottom1, bottom, right)
                    ladder_polygons.append(polygon)
                elif up1[0] == bottom[0] and up2[0] < right[0]:
                    polygon = (up1, bottom1, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif up1[0] == bottom[0]and right[0] == up2[0]:
                    polygon = (up1, bottom1, up2)
                    ladder_polygons.append(polygon)
                elif bottom[0] < up1[0] and up2[0] < right[0]:
                    polygon = (up1, bottom1, bottom2, up2)
                    ladder_polygons.append(polygon)
                elif bottom[0] < up1[0] and right[0] == up2[0]:
                    polygon = (up1, bottom1, up2)
                    ladder_polygons.append(polygon)
                elif up1[0] == right[0]:
                    polygon = (up1, up1, up1)
                    ladder_polygons.append(polygon)

        polygons.append(ladder_polygons)

    return polygons


def plot_polygons(polygons, label_cutoff, label_cutoff_position):

    columns = len(polygons[0])
    rows = len(polygons)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(X_LIM)
    ax.set_ylim(Y_LIM)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    patches = []
    handles = []
    for pi in polygons:
        for pii in pi:
            polygon = Polygon(pii, True)
            patches.append(polygon)
            handles.append(polygon)
    # set handles for legends, using polygons on first ladder
    handles = handles[:columns]

    p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.9)
    colors = [c for i in range(rows) for c in COLORS]
    p.set_color(colors)
    p.set_edgecolor(EDGE_COLOR)
    ax.add_collection(p)

    # put on legend
    ax.legend(handles, LABELS, loc='lower left', frameon=True)
    legends = ax.get_legend()
    for i in range(columns):
        legends.legendHandles[i].set_color(COLORS[i])

    # put on label_cutoff
    for l, l_p in zip(label_cutoff, label_cutoff_position):
        ax.text(l_p[0] + 5, l_p[1] - 2, l)
    fig.savefig('polygons.png')

    plt.show()


def main():
    ladder_coordiate = get_ladder_coordiate(
        LEFT, BOTTOM, RIGHT, input_ratio_matrix)
    label_cutoff = [l[0] for l in input_ratio_matrix]
    label_cutoff = [str(s) for s in label_cutoff]
    label_cutoff_position = [l[1] for l in ladder_coordiate]
    polygons = get_polygons(input_ratio_matrix, ladder_coordiate)

    plot_polygons(polygons, label_cutoff, label_cutoff_position)

if __name__ == "__main__":
    main()
