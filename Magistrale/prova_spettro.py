# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 09:35:40 2019

@author: Mattia
"""

import numpy as np ; import h5py ; import math ; import pylab as plt ; import matplotlib.cm as cm

bzfM2=[]; byfM2=[]; bxfM2=[]; BfM2=[]       #MODULO QUADRO DEL B TOTALE (come somma dei moduli quadri di bx, by e bz)

for i in range(9):
    if i<10: c = str(i) ; a = '0'+c
    else: a = str(i)
    bzfM2.append( (np.absolute(np.fft.rfft2((h5py.File('out0%s.h5' %a)).get('bz')[...])))**2 )
    byfM2.append( (np.absolute(np.fft.rfft2((h5py.File('out0%s.h5' %a)).get('by')[...])))**2 )
    bxfM2.append( (np.absolute(np.fft.rfft2((h5py.File('out0%s.h5' %a)).get('bx')[...])))**2 )
    BfM2.append( bxfM2[i]+byfM2[i]+bzfM2[i] )

ny = len( bzfM2[0][0,:,0] ) ; nx = len( bzfM2[0][0,0,:] ) ; nz = len( bzfM2[0][:,0,0] )

#for i in range(9):
#    BfM2.append( np.zeros((nz,ny,nx))+1. )       #matrice di 1 per provare

#%%    questa è la mia griglia di kx e ky. dovranno partire da 1 sennò il primo k avrà modulo zero...
kx = np.zeros(nx, dtype=int)
ky = np.zeros(ny, dtype=int)

for i in range(0,nx):
    kx[i]=i
    ky[i]=i
for i in range(nx,ny):
    ky[i]=-ny+i    #serve per come mi fa le trasformate di Fourier python
    
#ora definisco k[y,x] (ricordando che python negli argomenti mette prima y e poi x)
k = np.zeros((ny,nx),dtype=float)
for j in range(ny):
    for i in range(nx):
        k[j,i] = math.sqrt(kx[i]**2+ky[j]**2)
 
#%% visualizzo k e phi per controllare cos'ho combinato      
#    plt.figure(1)
#    plt.contourf(k,255,cmap=cm.flag)
#    plt.title('modulo k')
#    plt.colorbar()

#%%   #prendere la diagonale fa un pò schifo.. quindi prendo la "base"
k_base = np.zeros(nx,dtype=float)
for i in range(0,nx):
    k_base[i] =  k[0,i]
    
#%%
spettro_B = []
for i in range(9):
    spettro_B.append( np.zeros(nx,dtype=float) )

for t in range(9):
    for i in range(0,nx):
        for j in range(0,ny):  
                l = np.argmin(np.abs((k_base-k[j,i])))
                spettro_B[t][l] += BfM2[t][0,j,i]
            
spettro_B_finale = []
for t in range(9):
    spettro_B_finale.append( np.delete(spettro_B[t],[nx-1] ))
    
#%% Bf modulo
for i in range(9):
    a = str(i)
    plt.figure(i)
    plt.loglog(k_base[0:nx-1],spettro_B_finale[i])
#    plt.ylim((10**(-9),10**9))
'''#########################################################################'''
#%%
    
spettro_f = []
for i in range(9):
    spettro_f.append( np.zeros(nx+1,dtype=float) )
    
for t in range(9):
    for i in range(0,nx):
        for j in range(0,ny):      
            l = np.argmin(np.abs(k_base-k[j,i]))
            dist = k_base[l]-k[j,i]
            if dist>0 :
                spettro_f[t][l] += BfM2[t][0,j,i]*(1-dist)
                spettro_f[t][l-1] += BfM2[t][0,j,i]*(dist)
            elif dist<0 :
                spettro_f[t][l] += BfM2[t][0,j,i]*(1-np.abs(dist))
                spettro_f[t][l+1] += BfM2[t][0,j,i]*(np.abs(dist))
            else:# (k_base - k[j,i])=0 :
                spettro_f[t][l] += BfM2[t][0,j,i]
                    
spettro_B_finale_2 = []
for t in range(9):
    spettro_B_finale_2.append(np.delete(spettro_f[t],[nx-1,nx])) 
#%% Bf modulo
for i in range(9):
    a = str(i)
    plt.figure(i)
    plt.loglog(k_base[0:nx-1],spettro_B_finale_2[i])
    plt.ylim((10**(-9),10**9))
    
##%% provo a fare un fit di k_base[i] per ricavarne la funzione
#    
#from scipy import optimize  
#
#def func_k(x, m, q):
#    return m*x + q  #fitto con una retta
#
#params, params_covariance = optimize.curve_fit(func_k, kx , k_base)
#
##m = params[0]      #per la mia retta
##q = params[1]
#
#k_doppio=np.zeros(n+1,float)        #raddoppio il numero di punti
#for i in range(2,n+3):
#    k_doppio[i-2]= i/2
#
#funzione_k_base_2 = params[0]*k_doppio + params[1]
#
#k_quadruplo=np.zeros(2*n+1,float)       #quadruplo il numero di punti
#for i in range(3,2*n+5):
#    k_quadruplo[i-4]= i/4
#
#funzione_k_base_4 = params[0]*k_quadruplo + params[1]
#
##%%
#spettro_B_interp_2 = []
#for i in range(9):
#    spettro_B_interp_2.append( np.zeros(int(n+1),dtype=float) )
#
#for t in range(9):
#    for i in range(0,nx):
#        for j in range(0,n):      
#            l = np.argmin(np.abs((funzione_k_base_2-k[j,i])))
#            spettro_B_interp_2[t][l] += BfM2[t][0,j,i]
##            spettro_B_sotto[t][l] += sotto_vero_rib[t][j,i]
#            
#spettro_B_finale_interp_2 = []
#for t in range(9):
#    spettro_B_finale_interp_2.append( np.delete(spettro_B_interp_2[t],[int(n)] ))
##%%
#for i in range(9):
#    a = str(i)
#    plt.figure(i)
#    plt.loglog(funzione_k_base_2[0:int(n)],spettro_B_finale_interp_2[i])
#
##%%
#spettro_B_interp_4 = []
#for i in range(9):
#    spettro_B_interp_4.append( np.zeros(int(2*n+1),dtype=float) )
#
#for t in range(9):
#    for i in range(0,nx):
#        for j in range(0,n):      
#            l = np.argmin(np.abs((funzione_k_base_4-k[j,i])))
#            spettro_B_interp_4[t][l] += BfM2[t][0,j,i]
##            spettro_B_sotto[t][l] += sotto_vero_rib[t][j,i]
#            
#spettro_B_finale_interp_4 = []
#for t in range(9):
#    spettro_B_finale_interp_4.append( np.delete(spettro_B_interp_4[t],[int(2*n)] ))
##%%
#for i in range(9):
#    a = str(i)
#    plt.figure(i)
#    plt.loglog(funzione_k_base_4[0:int(2*n)],spettro_B_finale_interp_4[i])
