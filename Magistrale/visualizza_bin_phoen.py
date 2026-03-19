# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 12:13:45 2019

@author: Mattia
"""
import sys
import numpy as np ; from scipy.io import FortranFile ; import pylab as plt ; import matplotlib.cm as cm ; import time

start=time.time()
	
directory='../512^3(init_turb,va=0.6,tmax=40)/'      ; N=1024 ; n=512 ; nnz=512

array=np.zeros((nnz,n,n,12),dtype=np.float32)
	
a=sys.argv[1]
b=int(sys.argv[2])
if b==1 : 
    save=True
else:
    save=False

print(' --- LOAD --- out%s.dat' %a)

f=FortranFile(directory+'out%s.dat' %a)

data0=f.read_reals(dtype=np.float32)
t_out=data0[0]
 
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

   
#rh = ( array[:,:,:,0] )   ; he=rh ; cosa='Densità' ; segno=1
#vz = ( array[:,:,:,3] )  ; he=rh ; cosa='Velocità z' ; segno=-1
# pe = ( array[:,:,:,4] )  ; he=pe ; cosa='Pressione' ; segno=1
# bz = ( array[:,:,:,11] ) ; b0=1.3 ; he= bz-b0 ; cosa='Campo Magnetico z (senza B0)' ; segno=-1
# by = ( array[:,:,:,10] ) he= bz-b0 ; cosa='Campo Magnetico y' ; segno=-1
# bx = ( array[:,:,:,9] ) he= bz-b0 ; cosa='Campo Magnetico x' ; segno=-1

# vx = ( array[:,:,:,1] ) ; vy = ( array[:,:,:,2] ) ; vz = ( array[:,:,:,3] ) 
# vtot= np.sqrt( vx**2 + vy**2 + vz**2 )  ; he=vtot ; cosa='Modulo Velocità' ; segno=1

# bx = ( array[:,:,:,9] ) ; by = ( array[:,:,:,10] ) ; bz = ( array[:,:,:,11] )
# btot= np.sqrt( bx**2 + by**2 + bz**2 )  ; he=btot ; cosa='Modulo Campo Magnetico' ; segno=1

# bx = ( array[:,:,:,9] ) ; by = ( array[:,:,:,10] ) ; bz = ( array[:,:,:,11] ) ; b0=1.3
# btot= np.sqrt( bx**2 + by**2 + (bz-b0)**2 )  ; he=btot ; cosa='Modulo Campo Magnetico (senza B0)' ; segno=1

bx = ( array[:,:,:,9] ) ; by = ( array[:,:,:,10] ) ; b0=1.3 ; bz = ( array[:,:,:,11] - b0) 
vx = ( array[:,:,:,1] ) ; vy = ( array[:,:,:,2] ) ; vz = ( array[:,:,:,3] ) 
ch = 2*( bx*vx + by*vy + vz*bz )/( bx**2 + by**2 + bz**2 + vx**2 + vy**2 + vz**2 )
he=ch ; cosa='Cross Helicity (senza B0)' ; segno=-1

del(data0,data01,data02,data03,ix1,ix2,iy1,iy2,iz1,iz2,arr,ixmin,ixmax,array,n,p,nx,ny,nz,nv,nnz,N,a)
f.close()

# plt.figure(0)
# y=0
# vmax=np.amax(he[:,y,:]) ; vmin=np.amin(he[:,y,:]) ; vrif=max(vmax,np.abs(vmin))
# if segno==1 : plt.contourf(he[:,y,:],levels=np.linspace(vmin,vmax,255),cmap=cm.seismic)
# elif segno==-1 :plt.contourf(he[:,y,:],levels=np.linspace(-vrif,vrif,255),cmap=cm.seismic)
# title=cosa+'. Piano zx, y=%i'%y
# plt.title(title)
# plt.colorbar()
# if (save==True): plt.savefig(directory+title+'.png', format='png',dpi=1000,bbox_inches='tight')

# plt.figure(1)
# y=256
# vmax=np.amax(he[:,y,:]) ; vmin=np.amin(he[:,y,:]) ; vrif=max(vmax,np.abs(vmin))
# if segno==1 : plt.contourf(he[:,y,:],levels=np.linspace(vmin,vmax,255),cmap=cm.seismic)
# elif segno==-1 :plt.contourf(he[:,y,:],levels=np.linspace(-vrif,vrif,255),cmap=cm.seismic)
# title=cosa+'. Piano zx, y=%i'%y
# plt.title(title)
# plt.colorbar()
# if (save==True): plt.savefig(directory+title+'.png', format='png',dpi=1000,bbox_inches='tight')

# plt.figure(2)
# x=0
# vmax=np.amax(he[:,:,x]) ; vmin=np.amin(he[:,:,x]) ; vrif=max(vmax,np.abs(vmin))
# if segno==1 : plt.contourf(he[:,:,x],levels=np.linspace(vmin,vmax,255),cmap=cm.seismic)
# elif segno==-1 :plt.contourf(he[:,:,x],levels=np.linspace(-vrif,vrif,255),cmap=cm.seismic)
# title=cosa+'. Piano zy, x=%i'%x
# plt.title(title)
# plt.colorbar()
# if (save==True): plt.savefig(directory+title+'.png', format='png',dpi=1000,bbox_inches='tight')

# plt.figure(3)
# x=256
# vmax=np.amax(he[:,:,x]) ; vmin=np.amin(he[:,:,x]) ; vrif=max(vmax,np.abs(vmin))
# if segno==1 : plt.contourf(he[:,:,x],levels=np.linspace(vmin,vmax,255),cmap=cm.seismic)
# elif segno==-1 :plt.contourf(he[:,:,x],levels=np.linspace(-vrif,vrif,255),cmap=cm.seismic)
# title=cosa+'. Piano zy, x=%i'%x
# plt.title(title)
# plt.colorbar()
# if (save==True): plt.savefig(directory+title+'.png', format='png',dpi=1000,bbox_inches='tight')

# bx=vx ; by=vy ; bz=vz ; btot=vtot 
plt.figure(4)
y=0 ; x=420
# plottobz=np.concatenate((bz[429:512,y,x],bz[0:41,y,x]),axis=None)
# plottoby=np.concatenate((by[429:512,y,x],by[0:41,y,x]),axis=None)
# plottobx=np.concatenate((bx[429:512,y,x],bx[0:41,y,x]),axis=None)
# plottobtot=np.concatenate((btot[429:512,y,x],btot[0:41,y,x]),axis=None)
plottoch=np.concatenate((ch[429:512,y,x],ch[0:41,y,x]),axis=None)
# z=np.concatenate((np.linspace(429,511,511-429+1),np.linspace(0,40,41)),axis=None)
z=np.linspace(1,512-429+41,512-429+41)
# plt.plot(z,plottobz-b0,'xkcd:blue',linewidth=0.8,label='bz')
# plt.plot(z,plottoby,'xkcd:red',linewidth=0.8,label='by')
# plt.plot(z,plottobx,'xkcd:green',linewidth=0.8,label='bx')
# plt.plot(z,plottobtot,'xkcd:black',linewidth=0.8,label='btot')
plt.plot(z,plottoch,'xkcd:black',linewidth=0.8,label='ch')
title='CH senza b0. Taglio z=429:41, y=%i, x=%i'%(y,x)
plt.title(title)
plt.legend()
if (save==True): plt.savefig(directory+title+'.png', format='png',dpi=300,bbox_inches='tight')
end=time.time()

plt.show()

print('  time  :  ',int(end-start), 'sec')

