# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 11:40:49 2019

@author: Mattia
"""

print('START')

import sys ; import numpy as np ; from scipy.io import FortranFile 
import time ; import h5py ; import funzioni_matti_CINECA as matti

start=time.time()

#dir_lavoro='/Users/Mattia/Desktop/risultati LDZ/256^3 (init_turb. va=0.6 tmax=40)/'
#dir_save='/'

dir_lavoro='/marconi_scratch/userexternal/mmencagl/LDZ3D/'
dir_save='/marconi/home/userexternal/mmencagl/risultati_analisi/256^3(va=0.6,tmax=40)post/'
#dir_save='./'

print('dir work : ',dir_lavoro)
print('dir save : ',dir_save)

N=1024 ; print('numero processori usati : ',N)        #quanti processori ho usato

n=256   ; print('dimensione x,y : ',n)
nnz=256 ; print('dimensione  z  : ',nnz)

d=1                              #ogni quanti tempi leggo?

ti=int(sys.argv[1])     ; print('ti :  ',ti)     #tempo iniziale
tf=int(sys.argv[2])     ; print('tf :  ',tf)     #tempo finale

tempi = int((tf-ti)/d) 
time.sleep(6)

for i in range(tempi+1):
    start1=time.time()
    array=np.zeros((nnz,n,n,12),dtype=np.float32) #ha bisogno di annullarlo ogni volta per motivi poco chiari.....
    
    if d*i+ti<1000: a = str(d*i+ti)
    if d*i+ti<100:  c = str(d*i+ti) ; a = '0'+c
    if d*i+ti<10:   c = str(d*i+ti) ; a = '00'+c
    
    print('---LOADING--- out%s.dat' %a)
    f=FortranFile(dir_lavoro+'out%s.dat' %a)
    
    data0=f.read_reals(dtype=np.float32)
    t_out = data0[0]
     
    data01=f.read_reals(dtype=np.float32)
    nx=int(data01[0]) ; ny=int(data01[1]) ; nz=int(data01[2]) ; nv=int(data01[3])
    
    for p in range(N): #poi uso 2p+1 in modo da andare da 3 a 33 prendendo solo i dispari; oppure prendo i 2p per andare da 2 a 32 prendendo solo i pari
        data02=f.read_reals(dtype=np.float32)
        ixmin=data02[0:3]-1 ; ixmax=data02[3:6]
        ix1=int(ixmin[0]); ix2=int(ixmax[0])
        iy1=int(ixmin[1]); iy2=int(ixmax[1])
        iz1=int(ixmin[2]); iz2=int(ixmax[2])
        data03=f.read_reals(dtype=np.float32)
        arr=np.transpose( np.reshape(data03,(nx,ny,nz,nv), order='F'), (2,1,0,3) ) #cosi lo rimetto z-y-x
        for j in range(12):
            array[iz1:iz2,iy1:iy2,ix1:ix2,j] = arr[:,:,:,j]
      
    startload=time.time()
    rh = array[:,:,:,0] ; print('---LOAD--- rh')
    vx = array[:,:,:,1] ; print('---LOAD--- vx')
    vy = array[:,:,:,2] ; print('---LOAD--- vy')
    vz = array[:,:,:,3] ; print('---LOAD--- vz')
    pe = array[:,:,:,4] ; print('---LOAD--- pe')
    bx = array[:,:,:,9] ;  print('---LOAD--- bx')
    by = array[:,:,:,10] ; print('---LOAD--- by')
    bz = array[:,:,:,11] ; print('---LOAD--- bz')
    endload=time.time()
    matti.timing(startload,endload)
	
    f.close()

    Jx , Jy , Jz = matti.Curl_migliorissimo(bx,by,bz)
    Ax , Ay , Az = matti.Potenziale_Vettore(Jx,Jy,Jz)
    
    helicityJ = np.sum(bx*Jx + by*Jy + bz*Jz)
    helicityB = np.sum(bx*Ax + by*Ay + bz*Az)

    del(Jx , Jy , Jz, Ax , Ay , Az)

    w = rh + 4*pe + bx*bx + by*by + bz*bz
    
    del(rh , pe)
    
    chw = matti.Media( (vx*bx + vy*by + vz*bz)/(np.sqrt(w)) )
    normw = matti.Media( vx*vx + vy*vy + vz*vz + (bx*bx + by*by + bz*bz)/(w) )
    ch = matti.Media( vx*bx + vy*by + vz*bz )
    norm = matti.Media( vx*vx + vy*vy + vz*vz  +  bx*bx + by*by + bz*bz )
    
    del(bx , by , bz , vx , vy , vz , w)
    
    chw = 2*chw/normw
    ch  = 2*ch/norm

    hf = h5py.File(dir_save+'out_analysis%s.h5' %a,'r+') 
    
    hf.create_dataset('helicityJ', data=np.float32(helicityJ))
    hf.create_dataset('helicityB', data=np.float32(helicityB))
    hf.create_dataset('chw_giusto', data=np.float32(chw))
    hf.create_dataset('ch_giusto', data=np.float32(ch))
    
    hf.close()   
      
    
    print('#################### out_analysis%s.h5 ####################' %a)
    
    end1=time.time()
    
    matti.timing(start1,end1)
    print('***********************************************************')
    matti.timing(start,end1)
    
print()
print('  -------------------------')
print('  |   END  OF  ANALYSIS   |')
print('  -------------------------')
end=time.time()
matti.timing(start,end)
