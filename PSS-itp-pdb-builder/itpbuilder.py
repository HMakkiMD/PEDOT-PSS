#!/usr/bin/env python
# coding: utf-8

# In[5]:


import random
from pdbbuilder import sequence
#sequence = ['beg']
#seqmid = ['mid']*27
#seqmid2 = ['mid2']*6
#seqend = ['end']
#seqmid.extend(seqmid2)

#random.shuffle(seqmid)
#seqmid.extend(seqend)
#sequence.extend(seqmid)
#print(sequence)
counter = 0
atomcounter=0
atoms=[]
bonds=[]
pairs=[]
angles=[]
dihedrals=[]
impropers=[]
for monomer in sequence:
    if monomer == 'beg':
        with open('pssbeg.itp') as f:
            file_list = f.readlines()
            bondcounter=atomcounter
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomlines = item.split()
                        atomlines[0] = str(atomcounter+1)
                        atomlines[2] = str(counter+1)
                        atomlines[5] = str(atomcounter+1)
                        atomline ='{0: <6}'.format(atomlines[0])+'{0: <6}'.format(atomlines[1])+'{0: <6}'.format(atomlines[2])+'{0: <6}'.format(atomlines[3])+'{0: <6}'.format(atomlines[4])+'{0: <6}'.format(atomlines[5])+'{0: >12}'.format(atomlines[6])+'{0: >12}'.format(atomlines[7])
                        atoms.append(atomline)
                        atomcounter += 1
                if  file_list.index(item) > file_list.index("[ bonds ]\n") and file_list.index(item) < file_list.index("[ pairs ]\n"):
                    if not item.startswith(';'):
                        bondlines = item.split()
                        bondlines[0] = str(int(bondlines[0])+bondcounter)
                        bondlines[1] = str(int(bondlines[1])+bondcounter)
                        bondline ='{0: <6}'.format(bondlines[0])+'{0: <6}'.format(bondlines[1])+'{0: <6}'.format(bondlines[2])
                        bonds.append(bondline)
                if  file_list.index(item) > file_list.index("[ pairs ]\n") and file_list.index(item) < file_list.index("[ angles ]\n"):
                    if not item.startswith(';'):
                        pairlines = item.split()
                        pairlines[0] = str(int(pairlines[0])+bondcounter)
                        pairlines[1] = str(int(pairlines[1])+bondcounter)
                        pairline = '{0: <6}'.format(pairlines[0])+'{0: <6}'.format(pairlines[1])+'{0: <6}'.format(pairlines[2])
                        pairs.append(pairline)
                if  file_list.index(item) > file_list.index("[ angles ]\n") and file_list.index(item) < file_list.index("[ dihedrals ]\n"):
                    if not item.startswith(';'):
                        anglelines = item.split()
                        anglelines[0] = str(int(anglelines[0])+bondcounter)
                        anglelines[1] = str(int(anglelines[1])+bondcounter)
                        anglelines[2] = str(int(anglelines[2])+bondcounter)
                        angleline = '{0: <6}'.format(anglelines[0])+'{0: <6}'.format(anglelines[1])+'{0: <6}'.format(anglelines[2])+'{0: <6}'.format(anglelines[3])
                        angles.append(angleline)
                if  file_list.index(item) > file_list.index("[ dihedrals ]\n") and file_list.index(item) < file_list.index("[ impropers ]\n"):
                    if not item.startswith(';'):
                        dihedrallines = item.split()
                        dihedrallines[0] = str(int(dihedrallines[0])+bondcounter)
                        dihedrallines[1] = str(int(dihedrallines[1])+bondcounter)
                        dihedrallines[2] = str(int(dihedrallines[2])+bondcounter)
                        dihedrallines[3] = str(int(dihedrallines[3])+bondcounter)
                        dihedralline = '{0: <6}'.format(dihedrallines[0])+'{0: <6}'.format(dihedrallines[1])+'{0: <6}'.format(dihedrallines[2])+'{0: <6}'.format(dihedrallines[3])+'{0: <6}'.format(dihedrallines[4])
                        dihedrals.append(dihedralline)
                if  file_list.index(item) > file_list.index("[ impropers ]\n"):
                    if not item.startswith(';'):
                        improperlines = item.split()
                        improperlines[0] = str(int(improperlines[0])+bondcounter)
                        improperlines[1] = str(int(improperlines[1])+bondcounter)
                        improperlines[2] = str(int(improperlines[2])+bondcounter)
                        improperlines[3] = str(int(improperlines[3])+bondcounter)
                        improperline = '{0: <6}'.format(improperlines[0])+'{0: <6}'.format(improperlines[1])+'{0: <6}'.format(improperlines[2])+'{0: <6}'.format(improperlines[3])+'{0: <6}'.format(improperlines[4])
                        impropers.append(improperline)
            counter += 1
    if monomer == 'mid':
        with open('pssoh.itp') as f:
            file_list = f.readlines()
            bondcounter=atomcounter
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomlines = item.split()
                        atomlines[0] = str(atomcounter+1)
                        atomlines[2] = str(counter+1)
                        atomlines[5] = str(atomcounter+1)
                        atomline ='{0: <6}'.format(atomlines[0])+'{0: <6}'.format(atomlines[1])+'{0: <6}'.format(atomlines[2])+'{0: <6}'.format(atomlines[3])+'{0: <6}'.format(atomlines[4])+'{0: <6}'.format(atomlines[5])+'{0: >12}'.format(atomlines[6])+'{0: >12}'.format(atomlines[7])
                        atoms.append(atomline)
                        atomcounter += 1
                if  file_list.index(item) > file_list.index("[ bonds ]\n") and file_list.index(item) < file_list.index("[ pairs ]\n"):
                    if not item.startswith(';'):
                        bondlines = item.split()
                        bondlines[0] = str(int(bondlines[0])+bondcounter)
                        bondlines[1] = str(int(bondlines[1])+bondcounter)
                        bondline ='{0: <6}'.format(bondlines[0])+'{0: <6}'.format(bondlines[1])+'{0: <6}'.format(bondlines[2])
                        bonds.append(bondline)
                if  file_list.index(item) > file_list.index("[ pairs ]\n") and file_list.index(item) < file_list.index("[ angles ]\n"):
                    if not item.startswith(';'):
                        pairlines = item.split()
                        pairlines[0] = str(int(pairlines[0])+bondcounter)
                        pairlines[1] = str(int(pairlines[1])+bondcounter)
                        pairline = '{0: <6}'.format(pairlines[0])+'{0: <6}'.format(pairlines[1])+'{0: <6}'.format(pairlines[2])
                        pairs.append(pairline)
                if  file_list.index(item) > file_list.index("[ angles ]\n") and file_list.index(item) < file_list.index("[ dihedrals ]\n"):
                    if not item.startswith(';'):
                        anglelines = item.split()
                        anglelines[0] = str(int(anglelines[0])+bondcounter)
                        anglelines[1] = str(int(anglelines[1])+bondcounter)
                        anglelines[2] = str(int(anglelines[2])+bondcounter)
                        angleline = '{0: <6}'.format(anglelines[0])+'{0: <6}'.format(anglelines[1])+'{0: <6}'.format(anglelines[2])+'{0: <6}'.format(anglelines[3])
                        angles.append(angleline)
                if  file_list.index(item) > file_list.index("[ dihedrals ]\n") and file_list.index(item) < file_list.index("[ impropers ]\n"):
                    if not item.startswith(';'):
                        dihedrallines = item.split()
                        dihedrallines[0] = str(int(dihedrallines[0])+bondcounter)
                        dihedrallines[1] = str(int(dihedrallines[1])+bondcounter)
                        dihedrallines[2] = str(int(dihedrallines[2])+bondcounter)
                        dihedrallines[3] = str(int(dihedrallines[3])+bondcounter)
                        dihedralline = '{0: <6}'.format(dihedrallines[0])+'{0: <6}'.format(dihedrallines[1])+'{0: <6}'.format(dihedrallines[2])+'{0: <6}'.format(dihedrallines[3])+'{0: <6}'.format(dihedrallines[4])
                        dihedrals.append(dihedralline)
                if  file_list.index(item) > file_list.index("[ impropers ]\n"):
                    if not item.startswith(';'):
                        improperlines = item.split()
                        improperlines[0] = str(int(improperlines[0])+bondcounter)
                        improperlines[1] = str(int(improperlines[1])+bondcounter)
                        improperlines[2] = str(int(improperlines[2])+bondcounter)
                        improperlines[3] = str(int(improperlines[3])+bondcounter)
                        improperline = '{0: <6}'.format(improperlines[0])+'{0: <6}'.format(improperlines[1])+'{0: <6}'.format(improperlines[2])+'{0: <6}'.format(improperlines[3])+'{0: <6}'.format(improperlines[4])
                        impropers.append(improperline)
            counter += 1
            
    if monomer == 'mid2':
        with open('psso.itp') as f:
            file_list = f.readlines()
            bondcounter=atomcounter
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomlines = item.split()
                        atomlines[0] = str(atomcounter+1)
                        atomlines[2] = str(counter+1)
                        atomlines[5] = str(atomcounter+1)
                        atomline ='{0: <6}'.format(atomlines[0])+'{0: <6}'.format(atomlines[1])+'{0: <6}'.format(atomlines[2])+'{0: <6}'.format(atomlines[3])+'{0: <6}'.format(atomlines[4])+'{0: <6}'.format(atomlines[5])+'{0: >12}'.format(atomlines[6])+'{0: >12}'.format(atomlines[7])
                        atoms.append(atomline)
                        atomcounter += 1
                if  file_list.index(item) > file_list.index("[ bonds ]\n") and file_list.index(item) < file_list.index("[ pairs ]\n"):
                    if not item.startswith(';'):
                        bondlines = item.split()
                        bondlines[0] = str(int(bondlines[0])+bondcounter)
                        bondlines[1] = str(int(bondlines[1])+bondcounter)
                        bondline ='{0: <6}'.format(bondlines[0])+'{0: <6}'.format(bondlines[1])+'{0: <6}'.format(bondlines[2])
                        bonds.append(bondline)
                if  file_list.index(item) > file_list.index("[ pairs ]\n") and file_list.index(item) < file_list.index("[ angles ]\n"):
                    if not item.startswith(';'):
                        pairlines = item.split()
                        pairlines[0] = str(int(pairlines[0])+bondcounter)
                        pairlines[1] = str(int(pairlines[1])+bondcounter)
                        pairline = '{0: <6}'.format(pairlines[0])+'{0: <6}'.format(pairlines[1])+'{0: <6}'.format(pairlines[2])
                        pairs.append(pairline)
                if  file_list.index(item) > file_list.index("[ angles ]\n") and file_list.index(item) < file_list.index("[ dihedrals ]\n"):
                    if not item.startswith(';'):
                        anglelines = item.split()
                        anglelines[0] = str(int(anglelines[0])+bondcounter)
                        anglelines[1] = str(int(anglelines[1])+bondcounter)
                        anglelines[2] = str(int(anglelines[2])+bondcounter)
                        angleline = '{0: <6}'.format(anglelines[0])+'{0: <6}'.format(anglelines[1])+'{0: <6}'.format(anglelines[2])+'{0: <6}'.format(anglelines[3])
                        angles.append(angleline)
                if  file_list.index(item) > file_list.index("[ dihedrals ]\n") and file_list.index(item) < file_list.index("[ impropers ]\n"):
                    if not item.startswith(';'):
                        dihedrallines = item.split()
                        dihedrallines[0] = str(int(dihedrallines[0])+bondcounter)
                        dihedrallines[1] = str(int(dihedrallines[1])+bondcounter)
                        dihedrallines[2] = str(int(dihedrallines[2])+bondcounter)
                        dihedrallines[3] = str(int(dihedrallines[3])+bondcounter)
                        dihedralline = '{0: <6}'.format(dihedrallines[0])+'{0: <6}'.format(dihedrallines[1])+'{0: <6}'.format(dihedrallines[2])+'{0: <6}'.format(dihedrallines[3])+'{0: <6}'.format(dihedrallines[4])
                        dihedrals.append(dihedralline)
                if  file_list.index(item) > file_list.index("[ impropers ]\n"):
                    if not item.startswith(';'):
                        improperlines = item.split()
                        improperlines[0] = str(int(improperlines[0])+bondcounter)
                        improperlines[1] = str(int(improperlines[1])+bondcounter)
                        improperlines[2] = str(int(improperlines[2])+bondcounter)
                        improperlines[3] = str(int(improperlines[3])+bondcounter)
                        improperline = '{0: <6}'.format(improperlines[0])+'{0: <6}'.format(improperlines[1])+'{0: <6}'.format(improperlines[2])+'{0: <6}'.format(improperlines[3])+'{0: <6}'.format(improperlines[4])
                        impropers.append(improperline)
            counter += 1
    if monomer == 'end':
        with open('pssend.itp') as f:
            file_list = f.readlines()
            bondcounter=atomcounter
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomlines = item.split()
                        atomlines[0] = str(atomcounter+1)
                        atomlines[2] = str(counter+1)
                        atomlines[5] = str(atomcounter+1)
                        atomline ='{0: <6}'.format(atomlines[0])+'{0: <6}'.format(atomlines[1])+'{0: <6}'.format(atomlines[2])+'{0: <6}'.format(atomlines[3])+'{0: <6}'.format(atomlines[4])+'{0: <6}'.format(atomlines[5])+'{0: >12}'.format(atomlines[6])+'{0: >12}'.format(atomlines[7])
                        atoms.append(atomline)
                        atomcounter += 1
                if  file_list.index(item) > file_list.index("[ bonds ]\n") and file_list.index(item) < file_list.index("[ pairs ]\n"):
                    if not item.startswith(';'):
                        bondlines = item.split()
                        bondlines[0] = str(int(bondlines[0])+bondcounter)
                        bondlines[1] = str(int(bondlines[1])+bondcounter)
                        bondline ='{0: <6}'.format(bondlines[0])+'{0: <6}'.format(bondlines[1])+'{0: <6}'.format(bondlines[2])
                        bonds.append(bondline)
                if  file_list.index(item) > file_list.index("[ pairs ]\n") and file_list.index(item) < file_list.index("[ angles ]\n"):
                    if not item.startswith(';'):
                        pairlines = item.split()
                        pairlines[0] = str(int(pairlines[0])+bondcounter)
                        pairlines[1] = str(int(pairlines[1])+bondcounter)
                        pairline = '{0: <6}'.format(pairlines[0])+'{0: <6}'.format(pairlines[1])+'{0: <6}'.format(pairlines[2])
                        pairs.append(pairline)
                if  file_list.index(item) > file_list.index("[ angles ]\n") and file_list.index(item) < file_list.index("[ dihedrals ]\n"):
                    if not item.startswith(';'):
                        anglelines = item.split()
                        anglelines[0] = str(int(anglelines[0])+bondcounter)
                        anglelines[1] = str(int(anglelines[1])+bondcounter)
                        anglelines[2] = str(int(anglelines[2])+bondcounter)
                        angleline = '{0: <6}'.format(anglelines[0])+'{0: <6}'.format(anglelines[1])+'{0: <6}'.format(anglelines[2])+'{0: <6}'.format(anglelines[3])
                        angles.append(angleline)
                if  file_list.index(item) > file_list.index("[ dihedrals ]\n") and file_list.index(item) < file_list.index("[ impropers ]\n"):
                    if not item.startswith(';'):
                        dihedrallines = item.split()
                        dihedrallines[0] = str(int(dihedrallines[0])+bondcounter)
                        dihedrallines[1] = str(int(dihedrallines[1])+bondcounter)
                        dihedrallines[2] = str(int(dihedrallines[2])+bondcounter)
                        dihedrallines[3] = str(int(dihedrallines[3])+bondcounter)
                        dihedralline = '{0: <6}'.format(dihedrallines[0])+'{0: <6}'.format(dihedrallines[1])+'{0: <6}'.format(dihedrallines[2])+'{0: <6}'.format(dihedrallines[3])+'{0: <6}'.format(dihedrallines[4])
                        dihedrals.append(dihedralline)
                if  file_list.index(item) > file_list.index("[ impropers ]\n"):
                    if not item.startswith(';'):
                        improperlines = item.split()
                        improperlines[0] = str(int(improperlines[0])+bondcounter)
                        improperlines[1] = str(int(improperlines[1])+bondcounter)
                        improperlines[2] = str(int(improperlines[2])+bondcounter)
                        improperlines[3] = str(int(improperlines[3])+bondcounter)
                        improperline = '{0: <6}'.format(improperlines[0])+'{0: <6}'.format(improperlines[1])+'{0: <6}'.format(improperlines[2])+'{0: <6}'.format(improperlines[3])+'{0: <6}'.format(improperlines[4])
                        impropers.append(improperline)
            counter += 1
#new bonds
begmonomer = 4
mid1monomer = 3
mid1monomer2=0
mid2monomer = 3
mid2monomer2=0
endmonomer = 0
atomcounter=0
addibonds=['none']*3
addbond=[]
from itertools import islice
for i, monomer in islice(enumerate(sequence),0,len(sequence)-1):
    if sequence[i] == 'beg' and sequence[i+1]== 'mid':
        with open('pssbeg.itp') as f:
            file_list = f.readlines()
            addibonds[0]=str(begmonomer+atomcounter+1)
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomcounter += 1
            addibonds[1]=str(mid1monomer2+atomcounter+1)
            addibonds[2]='1'
            addibondline = '{0: <6}'.format(addibonds[0])+'{0: <6}'.format(addibonds[1])+'{0: <6}'.format(addibonds[2])
            addbond.append(addibondline)
    if sequence[i] == 'beg' and sequence[i+1]== 'mid2':
        with open('pssbeg.itp') as f:
            file_list = f.readlines()
            addibonds[0]=str(begmonomer+atomcounter+1)
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomcounter += 1
            addibonds[1]=str(mid2monomer2+atomcounter+1)
            addibonds[2]='1'
            addibondline = '{0: <6}'.format(addibonds[0])+'{0: <6}'.format(addibonds[1])+'{0: <6}'.format(addibonds[2])
            addbond.append(addibondline)
    if sequence[i] == 'beg' and sequence[i+1]== 'end':
        with open('pssbeg.itp') as f:
            file_list = f.readlines()
            addibonds[0]=str(begmonomer+atomcounter+1)
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomcounter += 1
            addibonds[1]=str(endmonomer+atomcounter+1)
            addibonds[2]='1'
            addibondline = '{0: <6}'.format(addibonds[0])+'{0: <6}'.format(addibonds[1])+'{0: <6}'.format(addibonds[2])
            addbond.append(addibondline)
    if sequence[i] == 'mid' and sequence[i+1]== 'mid':
        with open('pssoh.itp') as f:
            file_list = f.readlines()
            addibonds[0]=str(mid1monomer+atomcounter+1)
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomcounter += 1
            addibonds[1]=str(mid1monomer2+atomcounter+1)
            addibonds[2]='1'
            addibondline = '{0: <6}'.format(addibonds[0])+'{0: <6}'.format(addibonds[1])+'{0: <6}'.format(addibonds[2])
            addbond.append(addibondline)
    if sequence[i] == 'mid' and sequence[i+1]== 'mid2':
        with open('pssoh.itp') as f:
            file_list = f.readlines()
            addibonds[0]=str(mid1monomer+atomcounter+1)
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomcounter += 1
            addibonds[1]=str(mid2monomer2+atomcounter+1)
            addibonds[2]='1'
            addibondline = '{0: <6}'.format(addibonds[0])+'{0: <6}'.format(addibonds[1])+'{0: <6}'.format(addibonds[2])
            addbond.append(addibondline)
    if sequence[i] == 'mid' and sequence[i+1]== 'end':
        with open('pssoh.itp') as f:
            file_list = f.readlines()
            addibonds[0]=str(mid1monomer+atomcounter+1)
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomcounter += 1
            addibonds[1]=str(endmonomer+atomcounter+1)
            addibonds[2]='1'
            addibondline = '{0: <6}'.format(addibonds[0])+'{0: <6}'.format(addibonds[1])+'{0: <6}'.format(addibonds[2])
            addbond.append(addibondline)
    if sequence[i] == 'mid2' and sequence[i+1]== 'mid':
        with open('psso.itp') as f:
            file_list = f.readlines()
            addibonds[0]=str(mid2monomer+atomcounter+1)
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomcounter += 1
            addibonds[1]=str(mid1monomer2+atomcounter+1)
            addibonds[2]='1'
            addibondline = '{0: <6}'.format(addibonds[0])+'{0: <6}'.format(addibonds[1])+'{0: <6}'.format(addibonds[2])
            addbond.append(addibondline)
    if sequence[i] == 'mid2' and sequence[i+1]== 'mid2':
        with open('psso.itp') as f:
            file_list = f.readlines()
            addibonds[0]=str(mid2monomer+atomcounter+1)
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomcounter += 1
            addibonds[1]=str(mid2monomer2+atomcounter+1)
            addibonds[2]='1'
            addibondline = '{0: <6}'.format(addibonds[0])+'{0: <6}'.format(addibonds[1])+'{0: <6}'.format(addibonds[2])
            addbond.append(addibondline)
    if sequence[i] == 'mid2' and sequence[i+1]== 'end':
        with open('psso.itp') as f:
            file_list = f.readlines()
            addibonds[0]=str(mid2monomer+atomcounter+1)
            for item in file_list:
                if  file_list.index(item) > file_list.index("[ atoms ]\n") and file_list.index(item) < file_list.index("[ bonds ]\n"):
                    if not item.startswith(';'):
                        atomcounter += 1
            addibonds[1]=str(endmonomer+atomcounter+1)
            addibonds[2]='1'
            addibondline = '{0: <6}'.format(addibonds[0])+'{0: <6}'.format(addibonds[1])+'{0: <6}'.format(addibonds[2])
            addbond.append(addibondline)
addangle=[]
addanglelines=['none']*4
adddihedral=[]
adddihedrallines=['none']*5
#new angles and some dihedrals
bonds.extend(addbond)
for each in addbond:
    addbondline = each.split()
    for item in bonds:
        bondlines = item.split()
        if addbondline[0]==bondlines[0]:
            if  bondlines[1]!=addbondline[1]:
                addanglelines[0]=str(bondlines[1])
                addanglelines[1]=str(bondlines[0])
                addanglelines[2]=str(addbondline[1])
                addanglelines[3]='1'
                addangleline= '{0: <6}'.format(addanglelines[0])+'{0: <6}'.format(addanglelines[1])+'{0: <6}'.format(addanglelines[2])+'{0: <6}'.format(addanglelines[3])
                addangle.append(addangleline)
            for seconditem in bonds:
                secondbondlines = seconditem.split()
                if bondlines[1]==secondbondlines[0] and secondbondlines[1]!=bondlines[0]and secondbondlines[0]!=addbondline[1]:
                            adddihedrallines[0]=str(secondbondlines[1])
                            adddihedrallines[1]=str(secondbondlines[0])
                            adddihedrallines[2]=str(bondlines[0])
                            adddihedrallines[3]=str(addbondline[1])
                            adddihedrallines[4]='9'
                            adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                            adddihedral.append(adddihedralline)
                if bondlines[1]==secondbondlines[1] and secondbondlines[0]!=bondlines[0] and secondbondlines[1]!=addbondline[1]:
                            adddihedrallines[0]=str(secondbondlines[0])
                            adddihedrallines[1]=str(secondbondlines[1])
                            adddihedrallines[2]=str(bondlines[0])
                            adddihedrallines[3]=str(addbondline[1])
                            adddihedrallines[4]='9'
                            adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                            adddihedral.append(adddihedralline)
        if addbondline[0]==bondlines[1]:
            if  bondlines[0]!=addbondline[1]:
                addanglelines[0]=str(bondlines[0])
                addanglelines[1]=str(bondlines[1])
                addanglelines[2]=str(addbondline[1])
                addanglelines[3]='1'
                addangleline= '{0: <6}'.format(addanglelines[0])+'{0: <6}'.format(addanglelines[1])+'{0: <6}'.format(addanglelines[2])+'{0: <6}'.format(addanglelines[3])
                addangle.append(addangleline)
            for seconditem in bonds:
                secondbondlines = seconditem.split()
                if bondlines[0]==secondbondlines[0] and secondbondlines[1]!=bondlines[1] and secondbondlines[0]!=addbondline[1]:
                            adddihedrallines[0]=str(secondbondlines[1])
                            adddihedrallines[1]=str(secondbondlines[0])
                            adddihedrallines[2]=str(bondlines[1])
                            adddihedrallines[3]=str(addbondline[1])
                            adddihedrallines[4]='9'
                            adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                            adddihedral.append(adddihedralline)
                if bondlines[0]==secondbondlines[1] and secondbondlines[0]!=bondlines[1] and secondbondlines[1]!=addbondline[1]:
                            adddihedrallines[0]=str(secondbondlines[0])
                            adddihedrallines[1]=str(secondbondlines[1])
                            adddihedrallines[2]=str(bondlines[1])
                            adddihedrallines[3]=str(addbondline[1])
                            adddihedrallines[4]='9'
                            adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                            adddihedral.append(adddihedralline)
        if addbondline[1]==bondlines[0]:
            if  bondlines[1]!=addbondline[0]:
                addanglelines[0]=str(bondlines[1])
                addanglelines[1]=str(bondlines[0])
                addanglelines[2]=str(addbondline[0])
                addanglelines[3]='1'
                addangleline= '{0: <6}'.format(addanglelines[0])+'{0: <6}'.format(addanglelines[1])+'{0: <6}'.format(addanglelines[2])+'{0: <6}'.format(addanglelines[3])
                addangle.append(addangleline)
            for seconditem in bonds:
                secondbondlines = seconditem.split()
                if bondlines[1]==secondbondlines[0] and secondbondlines[1]!=bondlines[0] and secondbondlines[0]!=addbondline[0]:
                            adddihedrallines[0]=str(secondbondlines[1])
                            adddihedrallines[1]=str(secondbondlines[0])
                            adddihedrallines[2]=str(bondlines[0])
                            adddihedrallines[3]=str(addbondline[0])
                            adddihedrallines[4]='9'
                            adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                            adddihedral.append(adddihedralline)
                if bondlines[1]==secondbondlines[1] and secondbondlines[0]!=bondlines[0] and secondbondlines[1]!=addbondline[0]:
                            adddihedrallines[0]=str(secondbondlines[0])
                            adddihedrallines[1]=str(secondbondlines[1])
                            adddihedrallines[2]=str(bondlines[0])
                            adddihedrallines[3]=str(addbondline[0])
                            adddihedrallines[4]='9'
                            adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                            adddihedral.append(adddihedralline)
        if addbondline[1]==bondlines[1]:
            if bondlines[0]!=addbondline[0]:
                addanglelines[0]=str(bondlines[0])
                addanglelines[1]=str(bondlines[1])
                addanglelines[2]=str(addbondline[0])
                addanglelines[3]='1'
                addangleline= '{0: <6}'.format(addanglelines[0])+'{0: <6}'.format(addanglelines[1])+'{0: <6}'.format(addanglelines[2])+'{0: <6}'.format(addanglelines[3])
                addangle.append(addangleline)
            for seconditem in bonds:
                secondbondlines = seconditem.split()
                if bondlines[0]==secondbondlines[0] and secondbondlines[1]!=bondlines[1] and secondbondlines[0]!=addbondline[0]:
                            adddihedrallines[0]=str(secondbondlines[1])
                            adddihedrallines[1]=str(secondbondlines[0])
                            adddihedrallines[2]=str(bondlines[1])
                            adddihedrallines[3]=str(addbondline[0])
                            adddihedrallines[4]='9'
                            adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                            adddihedral.append(adddihedralline)
                if bondlines[0]==secondbondlines[1] and secondbondlines[0]!=bondlines[1] and secondbondlines[1]!=addbondline[0]:
                            adddihedrallines[0]=str(secondbondlines[0])
                            adddihedrallines[1]=str(secondbondlines[1])
                            adddihedrallines[2]=str(bondlines[1])
                            adddihedrallines[3]=str(addbondline[0])
                            adddihedrallines[4]='9'
                            adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                            adddihedral.append(adddihedralline)
#add additional dihedrals
for each in addbond:
    addbondline = each.split()
    for item in bonds:
        bondlines = item.split()
        if addbondline[0]==bondlines[0] and bondlines[1]!=addbondline[1]:
            for seconditem in bonds:
                secondbondlines = seconditem.split()
                if addbondline[1]==secondbondlines[0] and secondbondlines[1]!=addbondline[0]:
                    adddihedrallines[0]=str(secondbondlines[1])
                    adddihedrallines[1]=str(addbondline[1])
                    adddihedrallines[2]=str(addbondline[0])
                    adddihedrallines[3]=str(bondlines[1])
                    adddihedrallines[4]='9'
                    adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                    adddihedral.append(adddihedralline)
                if addbondline[1]==secondbondlines[1] and secondbondlines[0]!=addbondline[0]:
                    adddihedrallines[0]=str(secondbondlines[0])
                    adddihedrallines[1]=str(addbondline[1])
                    adddihedrallines[2]=str(addbondline[0])
                    adddihedrallines[3]=str(bondlines[1])
                    adddihedrallines[4]='9'
                    adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                    adddihedral.append(adddihedralline)
        if addbondline[0]==bondlines[1] and bondlines[0]!=addbondline[1]:
            for seconditem in bonds:
                secondbondlines = seconditem.split()
                if addbondline[1]==secondbondlines[0] and secondbondlines[1]!=addbondline[0]:
                    adddihedrallines[0]=str(secondbondlines[1])
                    adddihedrallines[1]=str(addbondline[1])
                    adddihedrallines[2]=str(addbondline[0])
                    adddihedrallines[3]=str(bondlines[0])
                    adddihedrallines[4]='9'
                    adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                    adddihedral.append(adddihedralline)
                if addbondline[1]==secondbondlines[1] and secondbondlines[0]!=addbondline[0]:
                    adddihedrallines[0]=str(secondbondlines[0])
                    adddihedrallines[1]=str(addbondline[1])
                    adddihedrallines[2]=str(addbondline[0])
                    adddihedrallines[3]=str(bondlines[0])
                    adddihedrallines[4]='9'
                    adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                    adddihedral.append(adddihedralline)
        if addbondline[1]==bondlines[0] and bondlines[1]!=addbondline[0]:
            for seconditem in bonds:
                secondbondlines = seconditem.split()
                if addbondline[0]==secondbondlines[0] and secondbondlines[1]!=addbondline[1]:
                    adddihedrallines[0]=str(secondbondlines[1])
                    adddihedrallines[1]=str(addbondline[0])
                    adddihedrallines[2]=str(addbondline[1])
                    adddihedrallines[3]=str(bondlines[1])
                    adddihedrallines[4]='9'
                    adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                    adddihedral.append(adddihedralline)
                if addbondline[0]==secondbondlines[1] and secondbondlines[0]!=addbondline[1]:
                    adddihedrallines[0]=str(secondbondlines[0])
                    adddihedrallines[1]=str(addbondline[0])
                    adddihedrallines[2]=str(addbondline[1])
                    adddihedrallines[3]=str(bondlines[1])
                    adddihedrallines[4]='9'
                    adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                    adddihedral.append(adddihedralline)
        if addbondline[1]==bondlines[1] and bondlines[0]!=addbondline[0]:
            for seconditem in bonds:
                secondbondlines = seconditem.split()
                if addbondline[0]==secondbondlines[0] and secondbondlines[1]!=addbondline[1]:
                    adddihedrallines[0]=str(secondbondlines[1])
                    adddihedrallines[1]=str(addbondline[0])
                    adddihedrallines[2]=str(addbondline[1])
                    adddihedrallines[3]=str(bondlines[0])
                    adddihedrallines[4]='9'
                    adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                    adddihedral.append(adddihedralline)
                if addbondline[0]==secondbondlines[1] and secondbondlines[0]!=addbondline[1]:
                    adddihedrallines[0]=str(secondbondlines[0])
                    adddihedrallines[1]=str(addbondline[0])
                    adddihedrallines[2]=str(addbondline[1])
                    adddihedrallines[3]=str(bondlines[0])
                    adddihedrallines[4]='9'
                    adddihedralline='{0: <6}'.format(adddihedrallines[0])+'{0: <6}'.format(adddihedrallines[1])+'{0: <6}'.format(adddihedrallines[2])+'{0: <6}'.format(adddihedrallines[3])+'{0: <6}'.format(adddihedrallines[4])
                    adddihedral.append(adddihedralline)

angle1lines=['none']*5
toremove=[]
import itertools
#removing repitition of angles
addangle=list(dict.fromkeys(addangle))
for each in addangle:
    angle1 = each.split()
    index = addangle.index(each)
    for line in addangle[index:]:
        angle2 = line.split()
        if angle1[0]==angle2[2] and angle1[1]==angle2[1] and angle1[2]==angle2[0]:
                toremove.append(index)            

def delete_multiple_element(list_object, indices):
    indices = sorted(indices, reverse=True)
    for idx in indices:
        if idx < len(list_object):
            list_object.pop(idx)
delete_multiple_element(addangle, toremove)


dihed1lines=['none']*5
toremove=[]


#removing repitition of dihedrals
adddihedral=list(dict.fromkeys(adddihedral))
for each in adddihedral:
    dihed1 = each.split()
    index = adddihedral.index(each)
    for line in adddihedral[index:]:
        dihed2 = line.split()
        if dihed1[0]==dihed2[3] and dihed1[1]==dihed2[2] and dihed1[2]==dihed2[1] and dihed1[3]==dihed2[0]:
                toremove.append(index)            

def delete_multiple_element(list_object, indices):
    indices = sorted(indices, reverse=True)
    for idx in indices:
        if idx < len(list_object):
            list_object.pop(idx)
delete_multiple_element(adddihedral, toremove)

addpair=[]
addpairlines=['none']*3
#add new pairs
for each in adddihedral:
    dihed1 = each.split()
    addpairlines[0]=dihed1[0]
    addpairlines[1]=dihed1[3]
    addpairlines[2]='1'
    addpairline='{0: <6}'.format(addpairlines[0])+'{0: <6}'.format(addpairlines[1])+'{0: <6}'.format(addpairlines[2])
    addpair.append(addpairline)


with open('mypss.itp', "w") as f:
    f.write("[ moleculetype ]\n; Name            nrexcl\n PSS               3\n")
    f.write("[ atoms ]\n")
    f.writelines("%s\n" % l for l in atoms)
    f.write("[ bonds ]\n")
    f.writelines("%s\n" % l for l in bonds)
    f.write("[ pairs ]\n")
    f.writelines("%s\n" % l for l in pairs)
    f.writelines("%s\n" % l for l in addpair)
    f.write("[ angles ]\n")
    f.writelines("%s\n" % l for l in angles)
    f.writelines("%s\n" % l for l in addangle)
    f.write("[ dihedrals ]\n")
    f.writelines("%s\n" % l for l in dihedrals)
    f.writelines("%s\n" % l for l in adddihedral)
    f.write("[ dihedrals ]\n")
    f.writelines("%s\n" % l for l in impropers)


# In[ ]:





# In[ ]:




