# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:31:17 2019

@author: Mattia
"""

print('START')

import numpy as np ; from scipy.io import FortranFile ; from numba import jit
import time ; import h5py ; import funzioni_matti_CINECA as matti

start=time.time()

cartella='4096^2 (init_turb. va=0.65 tmax=5)'
dir_lavoro='/Users/Mattia/Desktop/risultati_LDZ/'+cartella+'/'
dir_save='/Users/Mattia/Desktop/risultati_LDZ/'+cartella+'/analisi/'
del cartella

d=1         #ogni quanti tempi leggo?
ti=0       #tempo iniziale
tf=20       #tempo finale

tempi = int((tf-ti)/d) 

#dir_lavoro='/marconi_scratch/userexternal/mmencagl/LDZ3D/'
#dir_lavoro='/marconi/home/userexternal/mmencagl/test_mom/4096^2 (init_turb. va=0.65 tmax=5)/'
#dir_save='/marconi/home/userexternal/mmencagl/risultati_analisi/4096^2(va=0.65,tmax=5)/'
#dir_save='./'
print()
print('dir work : ',dir_lavoro)
print('dir save : ',dir_save)
print()
N=256 ; print('numero processori usati : ',N)        #quanti processori ho usato
print()
n=4096   ; print('dimensione x,y : ',n)
nnz=1 ; print('dimensione  z  : ',nnz)
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
       
    rh = array[:,:,:,0] ; print('---LOAD--- rh')
    vx = array[:,:,:,1] ; print('---LOAD--- vx')
    vy = array[:,:,:,2] ; print('---LOAD--- vy')
    vz = array[:,:,:,3] ; print('---LOAD--- vz')
    pe = array[:,:,:,4] ; print('---LOAD--- pe')
    #5,6,7,8 ? tre campi elettrici e l'entropia?
    bx = array[:,:,:,9] ;  print('---LOAD--- bx')
    by = array[:,:,:,10] ; print('---LOAD--- by')
    bz = array[:,:,:,11] ; print('---LOAD--- bz')
    
    f.close()

    rms_rh = matti.Varianza(rh) ; rms_pe = matti.Varianza(pe)
    min_pe=np.amin(pe) ; max_pe=np.amax(pe)
    min_rh=np.amin(rh) ; max_rh=np.amax(rh)
    rms_rh = matti.Varianza(rh) ; rms_pe = matti.Varianza(pe)

    ''' cross helicity ############################################## '''
    w = rh + 4*pe + bx*bx + by*by + bz*bz
    
    chw = matti.Media( (vx*bx + vy*by + vz*bz)/(np.sqrt(w)) )
    normw = matti.Media( (vx*vx + vy*vy + vz*vz) + (bx*bx + by*by + bz*bz)/(w) )
    ch = matti.Media( (vx*bx + vy*by + vz*bz) )
    norm = matti.Media( (vx*vx + vy*vy + vz*vz) + (bx*bx + by*by + bz*bz) )
    chw = 2*chw/normw
    ch  = 2*ch/norm
    
    del(data0,data01,data02,data03,ix1,ix2,iy1,iy2,iz1,iz2,arr,ixmin,ixmax,array,p,nx,ny,nz,nv,pe,rh,w)
    
    ''' B e V ############################################## '''
    @jit(cache=True)
    def Sp_Fourier(fx,fy,fz):
        start=time.time()
        FfM2 = (np.absolute(np.fft.rfftn(fx)))**2 + (np.absolute(np.fft.rfftn(fy)))**2 + (np.absolute(np.fft.rfftn(fz)))**2 
        end=time.time()
        print()
        print('END - Sp_Fourier')
        matti.timing(start,end)
        return FfM2

    BfM2 = Sp_Fourier(bx,by,bz)
    VfM2 = Sp_Fourier(vx,vy,vz)
    
    spettro_B_rad , k = matti.radial_spectrum_2D(BfM2)
    spettro_V_rad , k = matti.radial_spectrum_2D(VfM2)
    
    rms_bx = matti.Varianza(bx) ; rms_by = matti.Varianza(by) ; rms_bz = matti.Varianza(bz)
    rms_vx = matti.Varianza(vx) ; rms_vy = matti.Varianza(vy) ; rms_vz = matti.Varianza(vz)
    
    max_bx=(np.amax(bx)) ; max_by=(np.amax(by)) ; max_bz=(np.amax(bz))
    max_vx=(np.amax(vx)) ; max_vy=(np.amax(vy)) ; max_vz=(np.amax(vz))
    
    ''' J e vor e helicity ############################################## '''
    
    Jx , Jy , Jz = matti.Curl_migliorissimo(bx,by,bz)
    vorx , vory, vorz = matti.Curl_migliorissimo(vx,vy,vz)
    
    max_Jx=(np.amax(Jx)) ; max_Jy=(np.amax(Jy)) ; max_Jz=(np.amax(Jz))
    #max_vorx=(np.amax(vorx)) ; max_vory=(np.amax(vory)) ; max_Jz=(np.amax(vorz))
    
    rms_Jx = matti.Varianza(Jx) ; rms_Jy = matti.Varianza(Jy) ; rms_Jz = matti.Varianza(Jz)
    rms_vorx = matti.Varianza(vorx) ; rms_vory = matti.Varianza(vory) ; rms_vorz = matti.Varianza(vorz)
    
    Ax , Ay , Az = matti.Potenziale_Vettore(Jx,Jy,Jz)
    
    helicityJ = np.sum(bx*Jx + by*Jy + bz*Jz)
    helicityB = np.sum(bx*Ax + by*Ay + bz*Az)
    
    ''' stampo Hdf5 ############################################## '''
    hf = h5py.File(dir_save+'out_analysis%s.h5' %a,'w')
    
    hf.create_dataset('spettro_B_rad', data=np.float32(spettro_B_rad))
    hf.create_dataset('spettro_V_rad', data=np.float32(spettro_V_rad))
    hf.create_dataset('max_Jx', data=np.float32(max_Jx))
    hf.create_dataset('max_Jy', data=np.float32(max_Jy))
    hf.create_dataset('max_Jz', data=np.float32(max_Jz))
    hf.create_dataset('max_bx', data=np.float32(max_bx))
    hf.create_dataset('max_by', data=np.float32(max_by))
    hf.create_dataset('max_bz', data=np.float32(max_bz))
    hf.create_dataset('max_vx', data=np.float32(max_vx))
    hf.create_dataset('max_vy', data=np.float32(max_vy))
    hf.create_dataset('max_vz', data=np.float32(max_vz))
    hf.create_dataset('min_rh', data=np.float32(min_rh))  
    hf.create_dataset('min_pe', data=np.float32(min_pe))  
    hf.create_dataset('max_rh', data=np.float32(max_rh))  
    hf.create_dataset('max_pe', data=np.float32(max_pe))  
    hf.create_dataset('rms_rh', data=np.float32(rms_rh))
    hf.create_dataset('rms_pe', data=np.float32(rms_pe))
    hf.create_dataset('rms_bx', data=np.float32(rms_bx))
    hf.create_dataset('rms_by', data=np.float32(rms_by))
    hf.create_dataset('rms_bz', data=np.float32(rms_bz))
    hf.create_dataset('rms_vx', data=np.float32(rms_vx))
    hf.create_dataset('rms_vy', data=np.float32(rms_vy))
    hf.create_dataset('rms_vz', data=np.float32(rms_vz))
    hf.create_dataset('rms_Jx', data=np.float32(rms_Jx))
    hf.create_dataset('rms_Jy', data=np.float32(rms_Jy))
    hf.create_dataset('rms_Jz', data=np.float32(rms_Jz))
    hf.create_dataset('rms_vorx', data=np.float32(rms_vorx))
    hf.create_dataset('rms_vory', data=np.float32(rms_vory))
    hf.create_dataset('rms_vorz', data=np.float32(rms_vorz))
    hf.create_dataset('t_out', data=np.float32(t_out))
    hf.create_dataset('k', data=np.float32(k))
    
    hf.create_dataset('helicityJ', data=np.float32(helicityJ))
    hf.create_dataset('helicityB', data=np.float32(helicityB))
    hf.create_dataset('chw', data=np.float32(chw))
    hf.create_dataset('ch', data=np.float32(ch))

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
