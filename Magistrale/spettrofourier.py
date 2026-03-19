# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 17:42:35 2019

@author: Mattia
"""

import numpy as np; import h5py; import pylab
import matplotlib.cm as cm; from matplotlib import ticker

#antiex = [];ex=[]

bxf=[];byf=[];bzf=[];vxf=[];vyf=[];vzf=[];rhf=[];pef=[];antibxf=[]

bxfM2=[];byfM2=[];bzfM2=[];vxfM2=[];vyfM2=[];vzfM2=[];rhfM2=[];pefM2=[];BfM2=[];zero=[]

tempi=12
for i in range(tempi+1):
    if i<10: c = str(i) ; a = '0'+c
    else: a = str(i)
    bxf.append( np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('bx')[...]) )
    byf.append( np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('by')[...]) )
    bzf.append( np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('bz')[...]) )
    vxf.append( np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('vx')[...]) )
    vyf.append( np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('vy')[...]) )
    vzf.append( np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('vz')[...]) )
    rhf.append( np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('rh')[...]) )
    pef.append( np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('pe')[...]) )
    
#    ex.append( np.zeros((128,128,65),dtype=complex) )
#    ex[i][3,2,0]=(208138.03 + 25671.379j)
#    antiex.append( np.fft.irfftn( ex[i] ) )
#    antibxf.append( np.fft.irfftn( bxf[i] ) )
    
    bxfM2.append( (np.absolute(bxf[i]))**2 )
    byfM2.append( (np.absolute(byf[i]))**2 )
    bzfM2.append( (np.absolute(bzf[i]))**2 )
    vxfM2.append( (np.absolute(vxf[i]))**2 )
    vyfM2.append( (np.absolute(vyf[i]))**2 )
    vzfM2.append( (np.absolute(vzf[i]))**2 )
    rhfM2.append( (np.absolute(rhf[i]))**2 )
    pefM2.append( (np.absolute(pef[i]))**2 )
    BfM2.append( bxfM2[i]+byfM2[i]+bzfM2[i] )

ny = len( bzf[0][0,:,0] )

#%%
#pylab.figure(0)
#pylab.contourf(antiex[0][:,:,0],255,cmap=cm.YlGnBu)
#pylab.title('antiex')
#pylab.colorbar()
#pylab.figure(1)
#pylab.contourf(antibxf[0][:,:,0],255,cmap=cm.YlGnBu)
#pylab.title('antibxf')
#pylab.colorbar()

#%% Bf modulo
for i in range(tempi+1):
    a = str(i)
    pylab.figure(i)
    #pylab.contourf(BfM2[i][0,:,:],255,cmap=cm.YlGnBu,locator=ticker.LogLocator())
    pylab.contourf(BfM2[i][0,:,:],255,cmap=cm.YlGnBu,norm = LogNorm())
    pylab.title('Bf%s modulo quadro' %a)
    pylab.colorbar()
    
#%% bxf modulo
for i in range(1):
    a = str(i)
    pylab.figure(i)
    pylab.contourf(bxfM2[i][:,:,0],255,cmap=cm.flag)
    pylab.title('bxf%s modulo' %a)
    pylab.colorbar()
    
#%% byf modulo
for i in range(tempi+1):
    a = str(i)
    pylab.figure(i)
    pylab.contourf(byfM2[i][0,:,:],255,cmap=cm.YlGnBu)
    pylab.title('byf%s modulo' %a)
    pylab.colorbar()
 
#%% bzf modulo
for i in range(1):
    a = str(i)
    pylab.figure(i)
    pylab.contourf(bzfM2[i][0,:,:],255,cmap=cm.YlGnBu)
    pylab.title('bzf%s modulo' %a)
    pylab.colorbar()
    
#%% vxf modulo
for i in range(tempi+1):
    a = str(i)
    pylab.figure(i)
    pylab.contourf(vxfM2[i][0,:,:],255,cmap=cm.YlOrRd)
    pylab.title('vxf%s modulo' %a)
    pylab.colorbar()
    
#%% vyf modulo
for i in range(tempi+1):
    a = str(i)
    pylab.figure(i)
    pylab.contourf(vyfM2[i][0,:,:],255,cmap=cm.YlOrRd)
    pylab.title('vyf%s modulo' %a)
    pylab.colorbar()
    
#%% vzf modulo
for i in range(tempi+1):
    a = str(i)
    pylab.figure(i)
    pylab.contourf(vzfM2[i][0,:,:],255,cmap=cm.YlOrRd)
    pylab.title('vzf%s modulo' %a)
    pylab.colorbar()

#%% rhf modulo
for i in range(tempi+1):
    a = str(i)
    pylab.figure(i)
    pylab.contourf(rhfM2[i][0,:,:],255,cmap=cm.Spectral)
    pylab.title('rhf%s modulo' %a)
    pylab.colorbar()
    
#%% pef modulo
for i in range(tempi+1):
    a = str(i)
    pylab.figure(i)
    pylab.contourf(pefM2[i][0,:,:],255,cmap=cm.Spectral)
    pylab.title('pef%s modulo' %a)
    pylab.colorbar()