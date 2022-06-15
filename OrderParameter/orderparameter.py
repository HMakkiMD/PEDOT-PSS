#!/usr/bin/env python
# coding: utf-8

# In[145]:


import sys
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import statistics
import networkx as nx
import itertools


line_numbers=[]
atomindex=[]
XS=[]
YS=[]
ZS=[]
# director direction is the unit vector of the grain orientation
director = np.array([1,0,0])
# number is the number of atoms in each pi stacking segment
#only the atoms one wants to consider/I consider heavy atoms
number = 9
#number of total monomers in one polymer
polymer=12
#index.ndx is the index file for atoms exist in conjugated structure 
#that one wants to consider/total atoms in index file should be devisible by "number"
indexfilename = input('Enter an index file name: ')
with open(indexfilename) as f:
    file_list = f.readlines()
    for item in file_list[1:]:
        index = item.split()
        for i in range(len(index)):
            line_numbers.append(int(index[i])+1)
if len(line_numbers)%number !=0:
    print('number of atoms in pi-stacking segment does not match the index file')
    sys.exit()
lines = []
structurefilename = input('Enter a gro file name: ')
with open(structurefilename) as f:
    file_list = f.readlines()
    for i, line in enumerate(file_list):
        if i in line_numbers:
            XS.append(float(line[20:28]))
            YS.append(float(line[28:36]))
            ZS.append(float(line[36:44]))
xs=[]
ys=[]
zs=[]
monomernumber=0
CENTROID=[]
NORMAL=[]
pedotdirector=np.array([])
for i in range(int(len(line_numbers)/number/polymer)):
    CENTroid=[]
    CENTROID=[]
    for j in range(polymer):
        xs=[]
        ys=[]
        zs=[]
        for k in range(number):
            xs.append(XS[k+j*number+i*polymer*number])
            ys.append(YS[k+j*number+i*polymer*number])
            zs.append(ZS[k+j*number+i*polymer*number])
        centroid = [sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)]
        CENTroid.append(centroid)
    CENTROID.append(CENTroid)
    datamean = np.squeeze(np.array(CENTROID)).mean(axis=0)
    uu, dd, vv = np.linalg.svd(np.squeeze(np.array(CENTROID)) - datamean)
    pedotdirector=np.append(pedotdirector,vv[0],axis=0)
PEDOTDIRECTOR=pedotdirector[0::3]
orderparameter = np.average((3*PEDOTDIRECTOR**2-1)/2)
orderparameter_std = np.std((3*PEDOTDIRECTOR**2-1)/2)
with open('orderparameter.txt', "a") as f:
    f.write("\n [ "+structurefilename+" ] \n")
    f.writelines("%s\n" % orderparameter)
    f.writelines("%s\n" % orderparameter_std)




