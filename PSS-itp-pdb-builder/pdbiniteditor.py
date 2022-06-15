#!/usr/bin/env python
# coding: utf-8

# In[1]:


x=[]
y=[]
z=[]
coord=[]
#pdb cleaner
with open('pssbeg-init.pdb') as f:
    file_list = f.readlines()
    for item in file_list:
        if not item.startswith(';'):
            coordinate = item.split()
            x.append(float(coordinate[5]))
            y.append(float(coordinate[6]))
            z.append(float(coordinate[7]))
with open('pssbeg-init.pdb') as f:
    file_list = f.readlines()
    for i, item in enumerate(file_list):
        if not item.startswith(';'):
            coordinate = item.split()
            coordinate[4]='{0: >3}'.format(str(coordinate[3]))
            coordinate[5]=str(format((x[i]-max(x)), '.3f'))
            coordinate[6]=str(format((y[i]-max(y)), '.3f'))
            coordinate[7]=str(format((z[i]-max(z)), '.3f'))
            coordinateline = '{0: <6}'.format(coordinate[0])+'{0: >5}'.format(str(coordinate[1]))+' '+'{0: <4}'.format(str(coordinate[2]))+' '+'{0: >3}'.format(str(coordinate[3]))+'{0: >6}'.format(str(coordinate[4]))+'   '+'{0: >8}'.format(str(coordinate[5]))+'{0: >8}'.format(str(coordinate[6]))+'{0: >8}'.format(str(coordinate[7]))+'{0: >6}'.format(str(coordinate[8]))+'{0: >6}'.format(str(coordinate[9]))+'{0: >12}'.format(str(coordinate[10]))
            coord.append(coordinateline)
with open('pssbeg.pdb', 'w') as f:
    f.writelines("%s\n" % l for l in coord)
    f.writelines('END\n')
    
x=[]
y=[]
z=[]
coord=[]
#pdb cleaner
with open('psso-init.pdb') as f:
    file_list = f.readlines()
    for item in file_list:
        if not item.startswith(';'):
            coordinate = item.split()
            x.append(float(coordinate[5]))
            y.append(float(coordinate[6]))
            z.append(float(coordinate[7]))
with open('psso-init.pdb') as f:
    file_list = f.readlines()
    for i, item in enumerate(file_list):
        if not item.startswith(';'):
            coordinate = item.split()
            coordinate[4]='{0: >3}'.format(str(coordinate[3]))
            coordinate[5]=str(format((x[i]-max(x)), '.3f'))
            coordinate[6]=str(format((y[i]-max(y)), '.3f'))
            coordinate[7]=str(format((z[i]-max(z)), '.3f'))
            coordinateline = '{0: <6}'.format(coordinate[0])+'{0: >5}'.format(str(coordinate[1]))+' '+'{0: <4}'.format(str(coordinate[2]))+' '+'{0: >3}'.format(str(coordinate[3]))+'{0: >6}'.format(str(coordinate[4]))+'   '+'{0: >8}'.format(str(coordinate[5]))+'{0: >8}'.format(str(coordinate[6]))+'{0: >8}'.format(str(coordinate[7]))+'{0: >6}'.format(str(coordinate[8]))+'{0: >6}'.format(str(coordinate[9]))+'{0: >12}'.format(str(coordinate[10]))
            coord.append(coordinateline)
with open('psso.pdb', 'w') as f:
    f.writelines("%s\n" % l for l in coord)
    f.writelines('END\n')

x=[]
y=[]
z=[]
coord=[]
#pdb cleaner
with open('pssoh-init.pdb') as f:
    file_list = f.readlines()
    for item in file_list:
        if not item.startswith(';'):
            coordinate = item.split()
            x.append(float(coordinate[5]))
            y.append(float(coordinate[6]))
            z.append(float(coordinate[7]))
with open('pssoh-init.pdb') as f:
    file_list = f.readlines()
    for i, item in enumerate(file_list):
        if not item.startswith(';'):
            coordinate = item.split()
            coordinate[4]='{0: >3}'.format(str(coordinate[3]))
            coordinate[5]=str(format((x[i]-max(x)), '.3f'))
            coordinate[6]=str(format((y[i]-max(y)), '.3f'))
            coordinate[7]=str(format((z[i]-max(z)), '.3f'))
            coordinateline = '{0: <6}'.format(coordinate[0])+'{0: >5}'.format(str(coordinate[1]))+' '+'{0: <4}'.format(str(coordinate[2]))+' '+'{0: >3}'.format(str(coordinate[3]))+'{0: >6}'.format(str(coordinate[4]))+'   '+'{0: >8}'.format(str(coordinate[5]))+'{0: >8}'.format(str(coordinate[6]))+'{0: >8}'.format(str(coordinate[7]))+'{0: >6}'.format(str(coordinate[8]))+'{0: >6}'.format(str(coordinate[9]))+'{0: >12}'.format(str(coordinate[10]))
            coord.append(coordinateline)
with open('pssoh.pdb', 'w') as f:
    f.writelines("%s\n" % l for l in coord)
    f.writelines('END\n')

x=[]
y=[]
z=[]
coord=[]
#pdb cleaner
with open('pssend-init.pdb') as f:
    file_list = f.readlines()
    for item in file_list:
        if not item.startswith(';'):
            coordinate = item.split()
            x.append(float(coordinate[5]))
            y.append(float(coordinate[6]))
            z.append(float(coordinate[7]))
with open('pssend-init.pdb') as f:
    file_list = f.readlines()
    for i, item in enumerate(file_list):
        if not item.startswith(';'):
            coordinate = item.split()
            coordinate[4]='{0: >3}'.format(str(coordinate[3]))
            coordinate[5]=str(format((x[i]-max(x)), '.3f'))
            coordinate[6]=str(format((y[i]-max(y)), '.3f'))
            coordinate[7]=str(format((z[i]-max(z)), '.3f'))
            coordinateline = '{0: <6}'.format(coordinate[0])+'{0: >5}'.format(str(coordinate[1]))+' '+'{0: <4}'.format(str(coordinate[2]))+' '+'{0: >3}'.format(str(coordinate[3]))+'{0: >6}'.format(str(coordinate[4]))+'   '+'{0: >8}'.format(str(coordinate[5]))+'{0: >8}'.format(str(coordinate[6]))+'{0: >8}'.format(str(coordinate[7]))+'{0: >6}'.format(str(coordinate[8]))+'{0: >6}'.format(str(coordinate[9]))+'{0: >12}'.format(str(coordinate[10]))
            coord.append(coordinateline)
with open('pssend.pdb', 'w') as f:
    f.writelines("%s\n" % l for l in coord)
    f.writelines('END\n')


# In[ ]:





# In[ ]:




