#!/usr/bin/env python
# coding: utf-8

# In[29]:


import random
import pdbiniteditor
sequence = ['beg']
seqmid = ['mid']*25
seqmid2 = ['mid2']*8
seqend = ['end']
seqmid.extend(seqmid2)
pdbs=[]
random.shuffle(seqmid)
seqmid.extend(seqend)
sequence.extend(seqmid)
print(sequence)
counter = 0
atomcounter=0
monomerlength=3.200
y_bondcarbon=-1.500
for monomer in sequence:
    if monomer == 'beg':
        with open('pssbeg.pdb') as f:
            file_list = f.readlines()
            for item in file_list:
                    if not item.startswith('END'):
                        coordinate = item.split()
                        coordinate[1]=atomcounter+1
                        coordinate[4]=counter+1
                        if (counter % 2) == 1:
                            coordinate[6]=format(-1*(y_bondcarbon-float(coordinate[6])), '.3f')
                        else:
                            coordinate[6]=format(y_bondcarbon-float(coordinate[6]), '.3f')
                        pdbline = '{0: <6}'.format(coordinate[0])+'{0: >5}'.format(str(coordinate[1]))+' '+'{0: <4}'.format(str(coordinate[2]))+' '+'{0: >3}'.format(str(coordinate[3]))+'{0: >6}'.format(str(coordinate[4]))+'   '+'{0: >8}'.format(str(coordinate[5]))+'{0: >8}'.format(str(coordinate[6]))+'{0: >8}'.format(str(coordinate[7]))+'{0: >6}'.format(str(coordinate[8]))+'{0: >6}'.format(str(coordinate[9]))+'{0: >12}'.format(str(coordinate[10]))
                        pdbs.append(pdbline)
                        atomcounter += 1
            counter += 1
    if monomer == 'mid':
        with open('pssoh.pdb') as f:
            file_list = f.readlines()
            for item in file_list:
                    if not item.startswith('END'):
                        coordinate = item.split()
                        coordinate[1]=atomcounter+1
                        coordinate[4]=counter+1
                        coordinate[5]=format(float(coordinate[5])+monomerlength*counter, '.3f')
                        if (counter % 2) == 1:
                            coordinate[6]=format(-1*(y_bondcarbon-float(coordinate[6])), '.3f')
                        else:
                            coordinate[6]=format(y_bondcarbon-float(coordinate[6]), '.3f')
                        pdbline = '{0: <6}'.format(coordinate[0])+'{0: >5}'.format(str(coordinate[1]))+' '+'{0: <4}'.format(str(coordinate[2]))+' '+'{0: >3}'.format(str(coordinate[3]))+'{0: >6}'.format(str(coordinate[4]))+'   '+'{0: >8}'.format(str(coordinate[5]))+'{0: >8}'.format(str(coordinate[6]))+'{0: >8}'.format(str(coordinate[7]))+'{0: >6}'.format(str(coordinate[8]))+'{0: >6}'.format(str(coordinate[9]))+'{0: >12}'.format(str(coordinate[10]))
                        pdbs.append(pdbline)
                        atomcounter += 1
            counter += 1
    if monomer == 'mid2':
         with open('psso.pdb') as f:
            file_list = f.readlines()
            for item in file_list:
                    if not item.startswith('END'):
                        coordinate = item.split()
                        coordinate[1]=atomcounter+1
                        coordinate[4]=counter+1
                        coordinate[5]=format(float(coordinate[5])+monomerlength*counter, '.3f')
                        if (counter % 2) == 1:
                            coordinate[6]=format(-1*(y_bondcarbon-float(coordinate[6])), '.3f')
                        else:
                            coordinate[6]=format(y_bondcarbon-float(coordinate[6]), '.3f')
                        pdbline = '{0: <6}'.format(coordinate[0])+'{0: >5}'.format(str(coordinate[1]))+' '+'{0: <4}'.format(str(coordinate[2]))+' '+'{0: >3}'.format(str(coordinate[3]))+'{0: >6}'.format(str(coordinate[4]))+'   '+'{0: >8}'.format(str(coordinate[5]))+'{0: >8}'.format(str(coordinate[6]))+'{0: >8}'.format(str(coordinate[7]))+'{0: >6}'.format(str(coordinate[8]))+'{0: >6}'.format(str(coordinate[9]))+'{0: >12}'.format(str(coordinate[10]))
                        pdbs.append(pdbline)
                        atomcounter += 1
            counter += 1
    if monomer == 'end':
        with open('pssend.pdb') as f:
            file_list = f.readlines()
            for item in file_list:
                    if not item.startswith('END'):
                        coordinate = item.split()
                        coordinate[1]=atomcounter+1
                        coordinate[4]=counter+1
                        coordinate[5]=format(float(coordinate[5])+monomerlength*counter, '.3f')
                        if (counter % 2) == 1:
                            coordinate[6]=format(-1*(y_bondcarbon-float(coordinate[6])), '.3f')
                        else:
                            coordinate[6]=format(y_bondcarbon-float(coordinate[6]), '.3f')
                        pdbline = '{0: <6}'.format(coordinate[0])+'{0: >5}'.format(str(coordinate[1]))+' '+'{0: <4}'.format(str(coordinate[2]))+' '+'{0: >3}'.format(str(coordinate[3]))+'{0: >6}'.format(str(coordinate[4]))+'   '+'{0: >8}'.format(str(coordinate[5]))+'{0: >8}'.format(str(coordinate[6]))+'{0: >8}'.format(str(coordinate[7]))+'{0: >6}'.format(str(coordinate[8]))+'{0: >6}'.format(str(coordinate[9]))+'{0: >12}'.format(str(coordinate[10]))
                        pdbs.append(pdbline)
                        atomcounter += 1
            counter += 1       


with open('mypss.pdb', "w") as f:
    f.writelines("%s\n" % l for l in pdbs)
    f.writelines('END\n')


# In[ ]:





# In[ ]:




