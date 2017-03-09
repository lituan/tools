#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
given uniprot accs, query organisms
"""

import sys
import os
import urllib
import urllib2
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3,venn2
from multiprocessing import Pool

def uniprot_query(accs):
    url = 'http://www.uniprot.org/?'
    data = {
        'query':' '.join(accs) ,
        'format':'tab'
        'columns':'id,entry_name,organism,organism-id,reviewed'
    }
    data = urllib.urlencode(data)
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    r = response.readlines()
    lines = set([line.rstrip('\r\n') for line in r])
