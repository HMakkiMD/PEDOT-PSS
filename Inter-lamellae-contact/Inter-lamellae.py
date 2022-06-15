#!/usr/bin/env python
# coding: utf-8

# In[21]:


import sys
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import statistics
import networkx as nx
import itertools


line_numbers=[]
line_numbers_sp2c=[]
atomindex=[]
XS=[]
YS=[]
ZS=[]
sp2c_coord=[]
# number is the number of atoms in each pi stacking segment
#only the atoms one wants to consider/I consider heavy atoms
number = 9
#number of total monomers in one polymer
polymer=12
#number of SP2 carbon in each monomer
sp2c=4
#cutoff range to see how parallel are normals to the plane of each monomer
parallelrange = 0.04
#cutoff radius for finding possbible CENTROID neighbors
cutoff = 1.5 #nanometer
#cutoff distance for pistacks CENTROID in vertical (normal to their plane) and horizental (parallel to their plane)
verticalcutoff = 0.4 #nanometer
horizentalcutoff = 0.5
#distance between COG of monomers between each network
networkconnectcutoff = 0.3
#index.ndx is the index file for atoms exist in conjugated structure 
#that one wants to consider/total atoms in index file should be devisible by "number"
indexfilename = input('Enter EDOT heavy atoms index file name: ')
with open(indexfilename) as f:
    file_list = f.readlines()
    for item in file_list[1:]:
        index = item.split()
        for i in range(len(index)):
            line_numbers.append(int(index[i])+1)
if len(line_numbers)%number !=0:
    print('number of atoms in pi-stacking segment does not match the index file')
    sys.exit()
indexsp2cfilename = input('Enter EDOT sp2 carbon atoms index file name: ')
with open(indexsp2cfilename) as f:
    file_list = f.readlines()
    for item in file_list[1:]:
        index = item.split()
        for i in range(len(index)):
            line_numbers_sp2c.append(int(index[i])+1)
lines = []
structurefilename = input('Enter a gro file name: ')
with open(structurefilename) as f:
    file_list = f.readlines()
    for i, line in enumerate(file_list):
        if i in line_numbers:
            XS.append(float(line[20:28]))
            YS.append(float(line[28:36]))
            ZS.append(float(line[36:44]))
    for i, line in enumerate(file_list):
        if i in line_numbers_sp2c:
            sp2c_coord.append([float(line[20:28]),float(line[28:36]),float(line[36:44])])
            
xs=[]
ys=[]
zs=[]
monomernumber=0
CENTROID=[]
NORMAL=[]
for i in range(int(len(line_numbers)/number)):
    xs=[]
    ys=[]
    zs=[]
    for j in range(number):
        xs.append(XS[j+i*number])
        ys.append(YS[j+i*number])
        zs.append(ZS[j+i*number])
    tmp_A = []
    tmp_b = []
    for i in range(len(xs)):
        tmp_A.append([xs[i], ys[i], 1])
        tmp_b.append(zs[i])
    b = np.matrix(tmp_b).T
    A = np.matrix(tmp_A)

    # Manual solution
    fit = (A.T * A).I * A.T * b
    centroid = [sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)]
    CENTROID.append(centroid)
    normal = [float(fit[0]), float(fit[1]),-1.0]
    tempabs =  [abs(ele) for ele in normal]
    normal = [float(fit[0])/max(tempabs), float(fit[1])/max(tempabs),-1.0/max(tempabs)]
    NORMAL.append(normal)
PARALLEL=[]
for i in range(0,len(NORMAL)):
    parallel=[]
    for j in range(0,len(NORMAL)):
        if i!=j and (1.0 - parallelrange) < np.dot(NORMAL[i],NORMAL[j])/(np.linalg.norm(NORMAL[i])*np.linalg.norm(NORMAL[j])) < (1.0 + parallelrange):
            parallel.append(j)
    PARALLEL.append(parallel)
DISTANCE=[]
NEIGHBOR=[]
for i in range(0,len(CENTROID)):
    distance=[]
    neighbor=[]
    for j in range(0,len(CENTROID)):

        if i!=j and math.dist(CENTROID[i],CENTROID[j]) < cutoff:
            if i//polymer != j//polymer:
                neighbor.append(j)
                distance.append(math.dist(CENTROID[i],CENTROID[j]))
    NEIGHBOR.append(neighbor)
    DISTANCE.append(distance)
PARALLELANDNEIGHBOR=[]
VERTICALDISTANCE=[]
HORIZENTALDISTANCE=[]
for i in range(0,len(NEIGHBOR)):
    parallelandneighbor=[]
    verticaldistance=[]
    horizentaldistance=[]
    for j in range(0,len(NEIGHBOR[i])):
        for k in range(0,len(PARALLEL[i])):
            if PARALLEL[i][k]==NEIGHBOR[i][j]:
                parallelandneighbor.append(PARALLEL[i][k])
                vector_1 = np.subtract(np.array(CENTROID[i]),np.array(CENTROID[PARALLEL[i][k]]))
                vector_2 = NORMAL[i]
                unit_vector_1 = vector_1 / np. linalg. norm(vector_1)
                unit_vector_2 = vector_2 / np. linalg. norm(vector_2)
                dot_product = np. dot(unit_vector_1, unit_vector_2)
                sine = math.sqrt(1-dot_product**2)
                verticaldistance.append(abs(dot_product)*(np.linalg.norm(np.array(CENTROID[i])-np.array(CENTROID[PARALLEL[i][k]]))))
                horizentaldistance.append(sine*(np.linalg.norm(np.array(CENTROID[i])-np.array(CENTROID[PARALLEL[i][k]]))))
    PARALLELANDNEIGHBOR.append(parallelandneighbor)
    VERTICALDISTANCE.append(verticaldistance)
    HORIZENTALDISTANCE.append(horizentaldistance)
PISTACKWITHREPEAT=[]
for i in range(len(PARALLELANDNEIGHBOR)):
    for j in range(len(PARALLELANDNEIGHBOR[i])):
        PISTACKWITHREPEAT.append([i, PARALLELANDNEIGHBOR[i][j],VERTICALDISTANCE[i][j],HORIZENTALDISTANCE[i][j]])
PISTACKNOLIMIT=[]
for i in range(len(PISTACKWITHREPEAT)-1):
    for j in range(i+1,len(PISTACKWITHREPEAT)):
        if PISTACKWITHREPEAT[i][0]==PISTACKWITHREPEAT[j][1] and PISTACKWITHREPEAT[i][1]==PISTACKWITHREPEAT[j][0]:
            verticaldistance=[PISTACKWITHREPEAT[i][2],PISTACKWITHREPEAT[j][2] ]
            horizentaldistance=[PISTACKWITHREPEAT[i][3],PISTACKWITHREPEAT[j][3] ]
            PISTACKNOLIMIT.append([PISTACKWITHREPEAT[i][0],PISTACKWITHREPEAT[i][1],statistics.mean(verticaldistance),statistics.mean(horizentaldistance)])
PISTACK=[]
for i in range(len(PISTACKNOLIMIT)):
    if PISTACKNOLIMIT[i][2] < verticalcutoff and PISTACKNOLIMIT[i][3] < horizentalcutoff:
        PISTACK.append(PISTACKNOLIMIT[i])
        
atoms=[]
for i in range(0,len(PISTACK)):
        for j in range(0,number):
            atoms.append(line_numbers[PISTACK[i][0]*number+j]-1)
        for k in range(0,number):
            atoms.append(line_numbers[PISTACK[i][1]*number+k]-1)

#writing the vertical and horizental distance of pi stacked monomers

INFO=[]
POLYMERREPEAT=[]
for i in range(len(PISTACK)):
    pistackline = '{0: <12}'.format(PISTACK[i][0])+'{0: <12}'.format(PISTACK[i][1])+'    '+'{:.4f}'.format(PISTACK[i][2])+'    '+'{0:.4f}'.format(PISTACK[i][3])
    INFO.append(pistackline)
    POLYMERREPEAT.append(tuple([PISTACK[i][0]//polymer,PISTACK[i][1]//polymer]))
POLYMER = list(dict.fromkeys(POLYMERREPEAT))

#writing the index file for the polymers who are in a pi stacking netwrok

G = nx.Graph(POLYMER)

NETWORK=list(nx.connected_components(G))

networklist=list(itertools.chain(*NETWORK)) # list of all elements in NETWORK

isolatedchains=[x for x in list(range(int(len(line_numbers)/(number*polymer)))) if x not in networklist]


allnetworks = []
for i in range(len(NETWORK)):
    allnetworks.append(list(NETWORK[i]))
allnetworksnoiso=allnetworks.copy()
for i in range(len(isolatedchains)):
    allnetworks.append([isolatedchains[i]])
networkconnect = [[] for i in range(0, len(allnetworks))]
networkconnectcount = [0 for i in range(0, len(allnetworks))]
connectdistance = [[] for i in range(0, len(allnetworks))]
connectsp2c1 = [[] for i in range(0, len(allnetworks))]
connectsp2c2 = [[] for i in range(0, len(allnetworks))]

for i in range(len(allnetworks)):
    for j in range(len(allnetworks[i])):
        for k in range(len(allnetworks)):
            if k == i:
                continue
            for l in range(len(allnetworks[k])):
                for m in range(polymer*sp2c):
                    for n in range(polymer*sp2c):
                        if math.dist(sp2c_coord[allnetworks[i][j]*polymer*sp2c+m],sp2c_coord[allnetworks[k][l]*polymer*sp2c+n])<networkconnectcutoff:
                            networkconnect[i].append(k)
                            networkconnectcount[i]+=1
                            connectdistance[i].append(math.dist(sp2c_coord[allnetworks[i][j]*polymer*sp2c+m],sp2c_coord[allnetworks[k][l]*polymer*sp2c+n]))
                            connectsp2c1[i].append(allnetworks[i][j]*polymer*sp2c+m)
                            connectsp2c2[i].append(allnetworks[k][l]*polymer*sp2c+n)
#calculate contacts without considering individual chains                            
networkconnectnoiso = [[] for i in range(0, len(allnetworks))]
networkconnectcountnoiso = [0 for i in range(0, len(allnetworks))]
connectdistancenoiso = [[] for i in range(0, len(allnetworks))]
connectsp2c1noiso = [[] for i in range(0, len(allnetworks))]
connectsp2c2noiso = [[] for i in range(0, len(allnetworks))]

for i in range(len(allnetworksnoiso)):
    for j in range(len(allnetworksnoiso[i])):
        for k in range(len(allnetworksnoiso)):
            if k == i:
                continue
            for l in range(len(allnetworksnoiso[k])):
                for m in range(polymer*sp2c):
                    for n in range(polymer*sp2c):
                        if math.dist(sp2c_coord[allnetworksnoiso[i][j]*polymer*sp2c+m],sp2c_coord[allnetworksnoiso[k][l]*polymer*sp2c+n])<networkconnectcutoff:
                            networkconnectnoiso[i].append(k)
                            networkconnectcountnoiso[i]+=1
                            connectdistancenoiso[i].append(math.dist(sp2c_coord[allnetworksnoiso[i][j]*polymer*sp2c+m],sp2c_coord[allnetworksnoiso[k][l]*polymer*sp2c+n]))
                            connectsp2c1noiso[i].append(allnetworksnoiso[i][j]*polymer*sp2c+m)
                            connectsp2c2noiso[i].append(allnetworksnoiso[k][l]*polymer*sp2c+n)


import copy
networkconnectanalysis = copy.deepcopy(networkconnect)
networkconnectanalysiscopy = copy.deepcopy(networkconnect)

for i in range(len(networkconnectanalysis)):
    networkconnectanalysis[i].append(i)
    networkconnectanalysiscopy[i].append(i)
for i in range(len(networkconnectanalysis)):
    for j in range(len(networkconnectanalysis)):
        if i==j:
            continue
        if bool(set(networkconnectanalysiscopy[i]) & set(networkconnectanalysiscopy[j]))==True:
            networkconnectanalysis[i].extend(networkconnectanalysiscopy[j])
networkconnectanalysisTUPLE=[]
for i in range(len(networkconnectanalysis)):
    networkconnectanalysis[i]=list(dict.fromkeys(networkconnectanalysis[i]))
    networkconnectanalysis[i].sort()
    TUPLE=tuple(networkconnectanalysis[i])
    networkconnectanalysisTUPLE.append(TUPLE)

networkconnectanalysisnew = copy.deepcopy(networkconnectanalysis)
networkconnectanalysiscopynew = copy.deepcopy(networkconnectanalysis)

for i in range(len(networkconnectanalysisnew)):
    for j in range(len(networkconnectanalysisnew)):
        if i==j:
            continue
        if bool(set(networkconnectanalysiscopynew[i]) & set(networkconnectanalysiscopynew[j]))==True:
            networkconnectanalysisnew[i].extend(networkconnectanalysiscopynew[j])
networkconnectanalysisTUPLEnew=[]
for i in range(len(networkconnectanalysisnew)):
    networkconnectanalysisnew[i]=list(dict.fromkeys(networkconnectanalysisnew[i]))
    networkconnectanalysisnew[i].sort()
    TUPLEnew=tuple(networkconnectanalysisnew[i])
    networkconnectanalysisTUPLEnew.append(TUPLEnew)
finalnetwork=list(dict.fromkeys(networkconnectanalysisTUPLEnew))

networkconnectanalysisnewnew = copy.deepcopy(networkconnectanalysisnew)
networkconnectanalysiscopynewnew = copy.deepcopy(networkconnectanalysisnew)

for i in range(len(networkconnectanalysisnewnew)):
    for j in range(len(networkconnectanalysisnewnew)):
        if i==j:
            continue
        if bool(set(networkconnectanalysiscopynewnew[i]) & set(networkconnectanalysiscopynewnew[j]))==True:
            networkconnectanalysisnewnew[i].extend(networkconnectanalysiscopynewnew[j])
networkconnectanalysisTUPLEnewnew=[]
for i in range(len(networkconnectanalysisnewnew)):
    networkconnectanalysisnewnew[i]=list(dict.fromkeys(networkconnectanalysisnewnew[i]))
    networkconnectanalysisnewnew[i].sort()
    TUPLEnewnew=tuple(networkconnectanalysisnewnew[i])
    networkconnectanalysisTUPLEnewnew.append(TUPLEnewnew)
finalfinalnetwork=list(dict.fromkeys(networkconnectanalysisTUPLEnewnew))

LARGESTCLUSTER=[]
NUMBEROFCONTACT=[]
NETWORKCONNECTCUTOFF=[]
for i in range(71):
    networkconnect = [[] for i in range(0, len(allnetworks))]
    networkconnectcount = [0 for i in range(0, len(allnetworks))]
    connectdistance = [[] for i in range(0, len(allnetworks))]
    connectsp2c1 = [[] for i in range(0, len(allnetworks))]
    connectsp2c2 = [[] for i in range(0, len(allnetworks))]
    networkconnectcutoff = round(0.30+i/100,2)
    NETWORKCONNECTCUTOFF.append(networkconnectcutoff)
    for i in range(len(allnetworks)):
        for j in range(len(allnetworks[i])):
            for k in range(len(allnetworks)):
                if k == i:
                    continue
                for l in range(len(allnetworks[k])):
                    for m in range(polymer*sp2c):
                        for n in range(polymer*sp2c):
                            if math.dist(sp2c_coord[allnetworks[i][j]*polymer*sp2c+m],sp2c_coord[allnetworks[k][l]*polymer*sp2c+n])<networkconnectcutoff:
                                networkconnect[i].append(k)
                                networkconnectcount[i]+=1
                                connectdistance[i].append(math.dist(sp2c_coord[allnetworks[i][j]*polymer*sp2c+m],sp2c_coord[allnetworks[k][l]*polymer*sp2c+n]))
                                connectsp2c1[i].append(allnetworks[i][j]*polymer*sp2c+m)
                                connectsp2c2[i].append(allnetworks[k][l]*polymer*sp2c+n)
    networkconnectanalysis = copy.deepcopy(networkconnect)
    networkconnectanalysiscopy = copy.deepcopy(networkconnect)

    for i in range(len(networkconnectanalysis)):
        networkconnectanalysis[i].append(i)
        networkconnectanalysiscopy[i].append(i)
    for i in range(len(networkconnectanalysis)):
        for j in range(len(networkconnectanalysis)):
            if i==j:
                continue
            if bool(set(networkconnectanalysiscopy[i]) & set(networkconnectanalysiscopy[j]))==True:
                networkconnectanalysis[i].extend(networkconnectanalysiscopy[j])
    networkconnectanalysisTUPLE=[]
    for i in range(len(networkconnectanalysis)):
        networkconnectanalysis[i]=list(dict.fromkeys(networkconnectanalysis[i]))
        networkconnectanalysis[i].sort()
        TUPLE=tuple(networkconnectanalysis[i])
        networkconnectanalysisTUPLE.append(TUPLE)

    networkconnectanalysisnew = copy.deepcopy(networkconnectanalysis)
    networkconnectanalysiscopynew = copy.deepcopy(networkconnectanalysis)

    for i in range(len(networkconnectanalysisnew)):
        for j in range(len(networkconnectanalysisnew)):
            if i==j:
                continue
            if bool(set(networkconnectanalysiscopynew[i]) & set(networkconnectanalysiscopynew[j]))==True:
                networkconnectanalysisnew[i].extend(networkconnectanalysiscopynew[j])
    networkconnectanalysisTUPLEnew=[]
    for i in range(len(networkconnectanalysisnew)):
        networkconnectanalysisnew[i]=list(dict.fromkeys(networkconnectanalysisnew[i]))
        networkconnectanalysisnew[i].sort()
        TUPLEnew=tuple(networkconnectanalysisnew[i])
        networkconnectanalysisTUPLEnew.append(TUPLEnew)
    finalnetwork=list(dict.fromkeys(networkconnectanalysisTUPLEnew))

    networkconnectanalysisnewnew = copy.deepcopy(networkconnectanalysisnew)
    networkconnectanalysiscopynewnew = copy.deepcopy(networkconnectanalysisnew)

    for i in range(len(networkconnectanalysisnewnew)):
        for j in range(len(networkconnectanalysisnewnew)):
            if i==j:
                continue
            if bool(set(networkconnectanalysiscopynewnew[i]) & set(networkconnectanalysiscopynewnew[j]))==True:
                networkconnectanalysisnewnew[i].extend(networkconnectanalysiscopynewnew[j])
    networkconnectanalysisTUPLEnewnew=[]
    for i in range(len(networkconnectanalysisnewnew)):
        networkconnectanalysisnewnew[i]=list(dict.fromkeys(networkconnectanalysisnewnew[i]))
        networkconnectanalysisnewnew[i].sort()
        TUPLEnewnew=tuple(networkconnectanalysisnewnew[i])
        networkconnectanalysisTUPLEnewnew.append(TUPLEnewnew)
    finalfinalnetwork=list(dict.fromkeys(networkconnectanalysisTUPLEnewnew))
    print(networkconnectcutoff)
    clustersize=[]
    NUMBEROFPEDOT=[]
    for i in range(len(allnetworks)):
        clustersize.append(len(allnetworks[i]))
    for i in range(len(finalfinalnetwork)):
        numberofpedot=[]
        for j in range(len(finalfinalnetwork[i])):
            numberofpedot.append(clustersize[list(finalfinalnetwork[i])[j]])
        NUMBEROFPEDOT.append(numberofpedot)
    sumnumberofpedot=[]
    for i in range(len(NUMBEROFPEDOT)):
        sumnumberofpedot.append(sum(NUMBEROFPEDOT[i]))
    largestcluster=max(sumnumberofpedot)
    LARGESTCLUSTER.append(largestcluster)
    NUMBEROFCONTACT.append(sum(networkconnectcount))
numcontanct = "\n".join("{} {}".format(x, y) for x, y in zip(NETWORKCONNECTCUTOFF, NUMBEROFCONTACT))
largestclust = "\n".join("{} {}".format(x, y) for x, y in zip(NETWORKCONNECTCUTOFF, LARGESTCLUSTER))
with open('NUMBEROFCONTACT.txt', "w") as f:
    f.write(numcontanct)
with open('LARGESTCLUSTER.txt', "w") as f:
    f.write(largestclust)
list_finalfinalnetwork = [list(elem) for elem in finalfinalnetwork]
with open('CONNECTEDNETWORK.txt', "w") as f:
    f.writelines("%s\n" % l for l in list_finalfinalnetwork)


