# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 10:47:31 2019

@author: Mattia
"""

import h5py

a=2
b=6

c=a+b
d=a*b

hf = h5py.File('./provetta.h5','w') #crea il file
hf.create_dataset('c', data=int(c))

hf.create_dataset('d', data=int(d))
hf.close() 
  
#%%

e=c*d

hf = h5py.File('./provetta.h5','r+') #lo legge e ci puo aggiungere roba
hf.create_dataset('e', data=int(e))
hf.close() 

#%%

h=h5py.File('./provetta.h5','r')
list(h) #ti dice che c'è dentro il file


