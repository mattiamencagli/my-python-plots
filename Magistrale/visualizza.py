# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 10:42:08 2019

@author: Mattia
"""

import h5py ; import matplotlib.cm as cm ; import pylab as plt ; import numpy as np ; import matplotlib as mplt

directory='/Users/Mattia/Desktop/Risultati L.Franci/256^3(init_onda_2D. LF. 256tempi)(Lz=64Lx)/'

bx=[]; by=[]; bz=[]; vx=[]; vy=[]; vz=[]; pe=[]; rh=[]
tempi = 12
t_out=np.zeros(tempi+1,dtype=float)
for i in range(tempi+1):
    if 2*i<10: c = str(2*i) ; a = '0'+c
    else: a = str(2*i)
    h=h5py.File(directory + 'out%s0.h5' %a)
    t_out[i]=h.get('time')[...]
    bx.append( h.get('bx')[...] )
    by.append( h.get('by')[...] )
    bz.append( h.get('bz')[...] )
    vx.append( h.get('vx')[...] )
    vy.append( h.get('vy')[...] )
    vz.append( h.get('vz')[...] )
    rh.append( h.get('rh')[...] )
    pe.append( h.get('pe')[...] )

n = len( bz[0][0,:,0] )

#%% ESEMPIO PER DARE IL RANGE CHE SI VUOLE ALLA COLORBAR
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
#    plt.contourf(bx[i][0,:,:],np.arange(-1., 1.,.05),cmap=cm.YlGnBu,vmin=-1.,vmax=1.)
    plt.contourf(bx[i][0,:,:],levels=np.linspace(-1., 1.,255),cmap=cm.YlGnBu)
    plt.title('bx%s' %a)
    plt.colorbar()
#%% bx
for i in range(9,10):
    a = str(i)
    plt.figure(i)
    plt.contourf(bx[i][0,:,:],255,cmap=cm.YlGnBu)
    plt.title('bx%s' %a)
    plt.colorbar()

#%% by
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(by[i][0,:,:],255,cmap=cm.YlGnBu)
    plt.title('by%s' %a)
    plt.colorbar()
    
#%% bz
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(bz[i][0,:,:],200,cmap=cm.YlGnBu)
    plt.title('bz%s' %a)
    plt.colorbar()

#levels=np.linspace(np.amin(bz[i][0,:,:]),np.amax(bz[i][0,:,:]),255)

#%% b perp
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(np.sqrt(by[i][0,:,:]**2+bx[i][0,:,:]**2),255,cmap=cm.YlGnBu)
    plt.title('bz%s' %a)
    plt.colorbar()

 #%% vx
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(vx[i][0,:,:],255,cmap=cm.YlOrRd)
    plt.title('vx%s' %a)
    plt.colorbar()
    
#%% vy
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(vy[i][0,:,:],255,cmap=cm.YlOrRd)
    plt.title('vy%s' %a)
    plt.colorbar()
    
#%% vz
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(vz[i][0,:,:],255,cmap=cm.YlOrRd)
    plt.title('vz%s' %a)
    plt.colorbar()
    
#%% pe
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(pe[i][0,:,:],255,cmap=cm.Spectral)
    plt.title('pe%s' %a)
    plt.colorbar()
    
#%% rh
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(rh[i][0,:,:],255,cmap=cm.Spectral_r)
    plt.title('rh%s' %a)
    plt.colorbar()
