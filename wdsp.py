#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
define baseclass Wdsp for dealing with WDSP output file
"""
from collections import OrderedDict


class Wdsp():

    """
    use wdsp output file as input, gives you pros, pro_hotspots, pro_blades
    pro_scores etc.
    """

    def __init__(self, wdsp_f):
        """TODO: to be defined1. """
        self.pros = []
        self.seqs = OrderedDict()
        self.hotspots = OrderedDict()
        self.blades = OrderedDict()
        self.wdsps = OrderedDict()
        self.scores = OrderedDict()
        self.blade_scores = OrderedDict()
        self.repeat_num = OrderedDict()
        self.tetrad_num = OrderedDict()
        self.repeats = OrderedDict()

        self.lines = wdsp_f.readlines()
        self.lines = [line.strip('\n') for line in self.lines]
        self.lines = filter(lambda x: len(x.split()) > 0, self.lines)

        self.get_wdsps()
        self.get_blades()
        self.get_hotspots()

    def get_wdsps(self):
        for line in self.lines:
            words = line.split()
            if words[0] == '>':
                pro_name = words[1]
                # pro_name = words[1].split('|')[2]
                self.pros.append(pro_name)
                self.scores[pro_name] = float(words[2])
                pro_wdsp = [line]
            elif len(words) > 4:
                pro_wdsp.append(line)
                self.wdsps[pro_name] = pro_wdsp

    def get_blades(self):
        for pro, wdsp in self.wdsps.iteritems():
            self.blades[pro] = [line.split()[3:-1] for line in wdsp[1:]]
            self.blade_scores[pro] = [
                float(line.split()[-1]) for line in wdsp]
            self.seqs[pro] = ''.join([''.join(blade)
                                      for blade in self.blades[pro]])
            self.repeat_num[pro] = len(wdsp) - 1
            self.tetrad_num[pro] = len(
                [s for s in self.blade_scores[pro] if s >= 44.0])
            self.repeats[pro] = [''.join(b) for b in self.blades[pro]]

    def get_hotspots(self):
        for pro, blades in self.blades.iteritems():
            hotspot = []
            for blade in blades:
                R1 = blade[2][1]
                R1_2 = blade[1][-1]
                if len(blade[5]) <= 5 and blade[5][1] == 'D':
                    D_1 = blade[5][0]
                elif len(blade[5]) == 3 or len(blade[5]) == 2:
                    D_1 = blade[5][0]
                elif 3 <= len(blade[5]) <= 5 and blade[5][2] == 'D':
                    D_1 = blade[5][1]
                elif 4 <= len(blade[5]) <= 5 and blade[5][3] == 'D':
                    D_1 = blade[5][2]
                elif 5 <= len(blade[5]) <= 5 and blade[5][4] == 'D':
                    Di_1 = blade[5][3]
                elif len(blade[5]) <= 5:
                    D_1 = blade[5][1]
                elif len(blade[5]) <= 7:
                    D_1 = blade[5][0]
                else:
                    D_1 = '*'
                hotspot.append(R1 + R1_2 + D_1)
            self.hotspots[pro] = hotspot


def main():
    with open('test.wdsp') as wdsp_f:
        wdsp = Wdsp(wdsp_f)
        # print wdsp.blade_scores
        print wdsp.scores
        # print wdsp.hotspots
        # print wdsp.seqs
        print wdsp.wdsps
        print wdsp.repeat_num
        print wdsp.tetrad_num

if __name__ == "__main__":
    main()
