#!/usr/bin/env python
# coding: utf-8



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
cutoff = 1.0 #nanometer
#cutoff distance for pistacks CENTROID in vertical (normal to their plane) and horizental (parallel to their plane)
verticalcutoff = 0.4 #nanometer
horizentalcutoff = 0.5 #nanometer
#distance between COG of monomers between each network
networkconnectcutoff = 0.3 #nanometer
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
            #lines.append(line.strip())
            #coordinate=line.split()
            XS.append(float(line[20:28]))
            YS.append(float(line[28:36]))
            ZS.append(float(line[36:44]))
    for i, line in enumerate(file_list):
        if i in line_numbers_sp2c:
            #lines.append(line.strip())
            #coordinate=line.split()
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
    #errors = b - A * fit
    #residual = np.linalg.norm(errors)
    centroid = [sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)]
    CENTROID.append(centroid)
    normal = [float(fit[0]), float(fit[1]),-1.0]
    tempabs =  [abs(ele) for ele in normal]
    normal = [float(fit[0])/max(tempabs), float(fit[1])/max(tempabs),-1.0/max(tempabs)]
    NORMAL.append(normal)
    #print(np.dot(normal,normal)/(np.linalg.norm(normal)*np.linalg.norm(normal)))
    #print("solution: %f x + %f y + %f = z" % (fit[0], fit[1], fit[2]))
PARALLEL=[]
for i in range(0,len(NORMAL)):
    parallel=[]
    for j in range(0,len(NORMAL)):
        if i!=j and (1.0 - parallelrange) < np.dot(NORMAL[i],NORMAL[j])/(np.linalg.norm(NORMAL[i])*np.linalg.norm(NORMAL[j])) < (1.0 + parallelrange):
            parallel.append(j)
    PARALLEL.append(parallel)
#print(PARALLEL)
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
#print(NEIGHBOR)
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
#writing the index file for the monomers who have at least one pi stacking interaction
with open('pistack.ndx', "w") as f:
    f.write("[ pi stack ]\n")
    f.writelines("%s\n" % l for l in atoms)
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
with open('pistack.txt', "w") as f:
    f.write("ring1   ring2   vertical distance    horizental distance \n")
    f.writelines("%s\n" % l for l in INFO)
with open('cluster.ndx', "w") as f:
    CLUSTER=[]
for i in range(len(NETWORK)):
    cluster=[]
    for j in range(len(NETWORK[i])):
        for k in range(polymer*number):
            cluster.append(line_numbers[list(NETWORK[i])[j]*number*polymer+k]-1)
    with open('cluster.ndx', "a") as f:
        f.write("\n [ cluster"+str(i)+" ] \n")
        f.writelines("%s\n" % l for l in cluster)
isolatedchains=[x for x in list(range(int(len(line_numbers)/(number*polymer)))) if x not in networklist]
for i in range(len(isolatedchains)):
    isolatedatoms=[]
    for j in range(polymer*number):
        isolatedatoms.append(line_numbers[isolatedchains[i]*number*polymer+j]-1)
    with open('isolatedchains.ndx', "a") as f:
        f.write("\n [ isolatedchains"+str(i)+" ] \n")
        f.writelines("%s\t" % l for l in isolatedatoms) 
with open('network.txt', "a") as f:
    f.write("\n [ Number of pedot networks:   "+str(len(NETWORK))+" ] \n")
    f.writelines("%s\t" % l for l in NETWORK)
    f.write("\n [ Number of isolated chains:   "+str(len(isolatedchains))+" ] \n")
    f.writelines("%s\t" % l for l in isolatedchains)
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
#find number of pi-stack for each polymer couple in each network
rep = [[] for i in range(0, int(len(line_numbers)/polymer/number))]
for i in range(len(PISTACK)):
    rep[PISTACK[i][0]//polymer].append(PISTACK[i][1]//polymer)
repetition=[]
from collections import Counter
for i in range(len(rep)):
    repetition.append(dict(Counter(rep[i])))
pistackrepetition = [[] for i in range(0, len(allnetworksnoiso))]
for i in range(len(allnetworksnoiso)):
    for j in range(len(allnetworksnoiso[i])):
        pistackrepetition[i].append(repetition[allnetworksnoiso[i][j]])
with open('networkconnect2.txt', "w") as f:
    f.write("\n [ cutoff"+str(networkconnectcutoff)+" ] \n")
    f.writelines("%s\t" % l for l in networkconnect) 
with open('information.txt', "w") as f:
    f.write("[averageclustersize(excl isochain) averageclustersize(inc isochain)] \n")
    f.write(str(np.average([len(i) for i in NETWORK]))+"    "+str(np.average([len(i) for i in allnetworks]))+"\n")
    f.write("[varclustersize(excl isochain) varclustersize(inc isochain)] \n")
    f.write(str(np.var([len(i) for i in NETWORK]))+"    "+str(np.var([len(i) for i in allnetworks]))+"\n")
    f.write("[average contact between clusters] \n")
    f.write(str(np.average(networkconnectcount))+'\n')
    f.write("[var contact between clusters] \n")
    f.write(str(np.var(networkconnectcount))+'\n')
    f.write("[number of pistack monomer          number of tunnel position] \n")
    f.write(str(len(PISTACK))+'    '+str(np.sum(networkconnectcount))+'\n')
    f.write("[number of tunnel position - no isolated chain is considered] \n")
    f.write(str(np.sum(networkconnectcountnoiso))+'\n')
    f.write("[repetition of polymer couple in each network no isolated chain] \n")
    f.writelines("%s\t" % l for l in pistackrepetition) 
