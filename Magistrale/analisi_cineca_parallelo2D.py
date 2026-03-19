print('START')

import sys ; import numpy as np ; from scipy.io import FortranFile #import matplotlib as mplt 
import time ; import h5py ; import funzioni_matti_CINECA as matti

start=time.time()

dir_lavoro='/Users/Mattia/Desktop/risultati LDZ/2048^2 (init_turb. va=0.5 tmax=36)/'
dir_save='/Users/Mattia/Desktop/risultati LDZ/2048^2 (init_turb. va=0.5 tmax=36)/analisi/'

#dir_lavoro='/marconi_scratch/userexternal/mmencagl/LDZ3D/'
#dir_lavoro='/marconi/home/userexternal/mmencagl/test_mom/4096^2 (init_turb. va=0.65 tmax=5)/'
#dir_save='/marconi/home/userexternal/mmencagl/risultati_analisi/4096^2(va=0.65,tmax=5)/'
#dir_save='./'

print('dir work : ',dir_lavoro)
print('dir save : ',dir_save)

N=256 ; print('numero processori usati : ',N)        #quanti processori ho usato

n=2048   ; print('dimensione x,y : ',n)
nnz=1 ; print('dimensione  z  : ',nnz)

array=np.zeros((nnz,n,n,12),dtype=np.float32) #ha bisogno di annullarlo ogni volta per motivi poco chiari.....

a = '005'#sys.argv[1] 

print('---LOADING--- out%s.dat' %a)

time.sleep(3)

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
#btot.append( np.sqrt(bx[i]**2 + by[i]**2 + (bz[i]-1.)**2) )
#vtot.append( np.sqrt(vx[i]**2 + vy[i]**2 + (vz[i]-0.)**2) )
#gamma.append( (1.0-vtot[i])**(-1/2) )

f.close()

del(data0,data01,data02,data03,ix1,ix2,iy1,iy2,iz1,iz2,arr,ixmin,ixmax,array,n,p,nx,ny,nz,nv,nnz,N)

min_pe=np.amin(pe) ; max_pe=np.amax(pe)
min_rh=np.amin(rh) ; max_rh=np.amax(rh)

del(pe,rh)

''' B e V ############################################## '''
print('building : BfM2 and VfM2')

BfM2 = (np.absolute(np.fft.rfftn(bx)))**2 + (np.absolute(np.fft.rfftn(by)))**2 + (np.absolute(np.fft.rfftn(bz)))**2 
VfM2 = (np.absolute(np.fft.rfftn(vx)))**2 + (np.absolute(np.fft.rfftn(vy)))**2 + (np.absolute(np.fft.rfftn(vz)))**2 

spettro_B_rad , k = matti.radial_spectrum_2D(BfM2)
spettro_V_rad , k = matti.radial_spectrum_2D(VfM2)

rms_bx = matti.Varianza(bx) ; rms_by = matti.Varianza(by) ; rms_bz = matti.Varianza(bz)
rms_vx = matti.Varianza(vx) ; rms_vy = matti.Varianza(vy) ; rms_vz = matti.Varianza(vz)

max_bx=(np.amax(bx)) ; max_by=(np.amax(by)) ; max_bz=(np.amax(bz))
max_vx=(np.amax(vx)) ; max_vy=(np.amax(vy)) ; max_vz=(np.amax(vz))

''' J e vor ############################################## '''

Jx , Jy , Jz = matti.Curl_migliorissimo(bx,by,bz)
vorx , vory, vorz = matti.Curl_migliorissimo(vx,vy,vz)

max_Jx=(np.amax(Jx)) ; max_Jy=(np.amax(Jy)) ; max_Jz=(np.amax(Jz))
#max_vorx=(np.amax(vorx)) ; max_vory=(np.amax(vory)) ; max_vorz=(np.amax(vorz))

rms_Jx = matti.Varianza(Jx) ; rms_Jy = matti.Varianza(Jy) ; rms_Jz = matti.Varianza(Jz)
rms_vorx = matti.Varianza(vorx) ; rms_vory = matti.Varianza(vory) ; rms_vorz = matti.Varianza(vorz)

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
hf.close()   
print('#################### out_analysis%s.h5 ####################' %a)

end=time.time()

print()
print('  -------------------------')
print('  |   END  OF  ANALYSIS   |')
print('  -------------------------')
matti.timing(start,end)
