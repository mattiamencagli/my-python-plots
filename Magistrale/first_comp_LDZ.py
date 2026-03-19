# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 13:50:34 2019

@author: Mattia
"""

import numpy as np ; import time ; from scipy.io import FortranFile ; import funzioni_matti as matti

start=time.time()

bx = [] ; by = [] ; bz = [] ; vx = [] ; vy = [] ; vz = [] ; pe = [] ; rh = [] ; btot = []

tempi=1
t_out=np.zeros(tempi+1,dtype=float)

directory='/Users/Mattia/Desktop/risultati LDZ/2048^2 (init_turb. 24tempi con tmax=12, res=0)/'

for i in range(tempi+1):
    if i<10: c = str(i) ; a = '0'+c
    else: a = str(i)
    f=FortranFile(directory+'out%s.dat' %a)
    
    data0=f.read_reals(dtype=np.float32)
    t_out[i]=data0[0]
    
    data1=f.read_reals(dtype=np.float32)
    nx=int(data1[0]) ; ny=int(data1[1]) ; nz=int(data1[2]) ; nv=int(data1[3])
    
    data2=f.read_reals(dtype=np.float32)
    ixmin=data2[0:3]-1 ; ixmax=data2[3:6]
    ix1=int(ixmin[0]) ; ix2=int(ixmax[0])
    iy1=int(ixmin[1]) ; iy2=int(ixmax[1])
    iz1=int(ixmin[2]) ; iz2=int(ixmax[2])
    
    data3=f.read_reals(dtype=np.float32)
    arr=np.transpose( np.reshape(data3,(nx,ny,nz,nv), order='F'), (2,1,0,3) ) #cosi lo rimetto z-y-x

    rh.append( arr[:,:,:,0] )
    vx.append( arr[:,:,:,1] )
    vy.append( arr[:,:,:,2] )
    vz.append( arr[:,:,:,3] )
    pe.append( arr[:,:,:,4] )
    
    #5,6,7,8 ? tre campi elettrici e l'entropia?
    bx.append( arr[:,:,:,9] )
    by.append( arr[:,:,:,10] )
    bz.append( arr[:,:,:,11] )
    btot.append( bx[i]**2 + by[i]**2 + bz[i]**2 )
    
    f.close()

end=time.time()

matti.timing(start,end)