#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
use mds to plot high-dimentional data
"""
from matplotlib import pyplot as plt
from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

import numpy as np

n = 3
seed = np.random.RandomState(seed=3)
x=seed.randint(0,20,2*n).astype(np.float)
x=x.reshape((n,2))
x-=x.mean()
similarities=euclidean_distances(x)
print 'original: ',x
print 'similarities: ',similarities

mds=manifold.MDS(n_components=2,max_iter=3000,eps=1e-9,random_state=seed,dissimilarity='precomputed',n_jobs=1)
pos = mds.fit(similarities).embedding_
print 'pos: ',pos
pos *= np.sqrt((x**2).sum())/np.sqrt((pos**2).sum())
clf=PCA(n_components=2)
x=clf.fit_transform(x)

fig=plt.figure(1)
ax=plt.axes([0.,0.,1.,1.])
plt.scatter(pos[:,0],pos[:,1],s=20,c='g')
plt.show()
