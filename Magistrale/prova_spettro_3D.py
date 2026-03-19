# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:25:48 2019

@author: Mattia
"""

import numpy as np ; import h5py ; import math ; import pylab as plt ; import matplotlib.cm as cm 
import matplotlib as mplt
bzfM2=[]; byfM2=[]; bxfM2=[]; BfM2=[]       #MODULO QUADRO DEL B TOTALE (come somma dei moduli quadri di bx, by e bz)
tempi=12
for i in range(tempi+1):
    if i<10: c = str(i) ; a = '0'+c
    else: a = str(i)
    bzfM2.append( (np.absolute(np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('bz')[...])))**2 )
    byfM2.append( (np.absolute(np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('by')[...])))**2 )
    bxfM2.append( (np.absolute(np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('bx')[...])))**2 )
    BfM2.append( bxfM2[i]+byfM2[i]+bzfM2[i] )

nx = len( bzfM2[0][0,0,:] ) ; ny = len( bzfM2[0][0,:,0] ) ; nz = len( bzfM2[0][:,0,0] )

#%%
kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)

#for i in range(0,nx):
#    kx[i]=i
#    ky[i]=i
#    kz[i]=i
#for i in range(nx,ny):
#    ky[i]=-ny+i   
#for i in range(nx,nz):
#    kz[i]=-ny+i
for i in range(0,nx):
    kx[i]=i ; ky[i]=i 
    if nz>nx: kz[i]=i
for i in range(nx,ny):
    ky[i]=-ny+i   
if nz>nx:
    for i in range(nx,nz):
        kz[i]=-ny+i  
#ora definisco k[y,x] (ricordando che python negli argomenti mette prima y e poi x)
k_perp = np.zeros((nz,ny),dtype=float)
for j in range(nz):
    for i in range(ny):
        k_perp[j,i] = math.sqrt(kz[j]**2+ky[i]**2)

k_base = np.zeros(nx,dtype=float)
for i in range(0,nx):
    k_base[i] =  k_perp[0,i]
#%% 
compt_x= [] ; compt_y = [] ; compt_z =[]
for t in range(tempi+1):
    if t<10: c = str(t) ; a = '0'+c
    else: a = str(t)
    compt_x.append(np.zeros((nz,ny,nx),dtype=complex))
    compt_y.append(np.zeros((nz,ny,nx),dtype=complex))
    compt_z.append(np.zeros((nz,ny,nx),dtype=complex))

J_x= [] ; J_y = [] ; J_z =[]
bxt= [] ; byt = [] ; bzt =[]
  
for t in range(tempi+1):
    if t<10: c = str(t) ; a = '0'+c
    else: a = str(t)
    for w in range(nz):
        for j in range(ny):
            for i in range(nx):
                bxt.append(np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('bx')[...]))
                byt.append(np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('by')[...])) 
                bzt.append(np.fft.rfftn((h5py.File('out0%s.h5' %a)).get('bz')[...]))
                compt_x[t][w,j,i] = ( ky[j]*bzt[t][w,j,i] - kz[w]*byt[t][w,j,i] )
                compt_y[t][w,j,i] = ( kz[w]*bxt[t][w,j,i] - kx[i]*bzt[t][w,j,i] )
                compt_z[t][w,j,i] = ( kx[i]*byt[t][w,j,i] - ky[j]*bxt[t][w,j,i] )
for t in range(tempi+1):    
    J_x.append(np.fft.irfftn(compt_x[t]))
    J_y.append(np.fft.irfftn(compt_y[t]))
    J_z.append(np.fft.irfftn(compt_z[t]))


#%%
#plt.figure(1)
#plt.contourf(k_perp,255,cmap=cm.flag)
#plt.title('modulo k_perp')
#plt.colorbar()    

#%%
spettro_B_gyro_pesi = []
for i in range(tempi+1):
    spettro_B_gyro_pesi.append( np.zeros((nx+1,nx),dtype=float) )

for t in range(tempi+1):
    for w in range(0,nx):
        for j in range(0,nz):
            for i in range(0,ny): 
                l = np.argmin(np.abs((k_base-k_perp[j,i])))
                dist = k_base[l]-k_perp[j,i]
                if dist>0 :
                    spettro_B_gyro_pesi[t][l,w] += BfM2[t][j,i,w]*(1-dist)
                    spettro_B_gyro_pesi[t][l-1,w] += BfM2[t][j,i,w]*(dist)
                elif dist<0 :
                    spettro_B_gyro_pesi[t][l,w] += BfM2[t][j,i,w]*(1-np.abs(dist))
                    spettro_B_gyro_pesi[t][l+1,w] += BfM2[t][j,i,w]*(np.abs(dist))
                else:# (k_base - k[j,i])=0 :
                    spettro_B_gyro_pesi[t][l,w] += BfM2[t][j,i,w]

spettro_B_gyro_pesi_finale = []
for t in range(tempi+1):
    spettro_B_gyro_pesi_finale.append( np.delete(spettro_B_gyro_pesi[t],nx,axis=0)) 

#spettro_B_gyro_pesi_finale = []
#for t in range(tempi+1):
#    spettro_B_gyro_pesi_finale.append( np.delete(spettro_B_gyro_pesi_fin[t],nx-1,axis=0)) 

#%% Bf modulo
prova_andamento=(1)*(kx)**(3/2)
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)#=np.logspace(-13, 13, 200)
#    plt.contourf(spettro_B_gyro_pesi_finale[i],cmap=cm.YlGnBu,norm=LogNorm())
    plt.contourf(spettro_B_gyro_pesi_finale[i],cmap=cm.jet,levels=np.logspace(-1, 13, 255),locator=mplt.ticker.LogLocator())
    mplt.pyplot.yscale('log') ; mplt.pyplot.xscale('log') ; plt.ylim((1,64)) ; plt.xlim((1,64))
    plt.xlabel('k parallel') ; plt.ylabel('k perpendicular') ; plt.colorbar() ; plt.title('Spettro B (con pesi)')

#    plt.plot(kx,prova_andamento,'k')
#%%
    
k_mod = np.zeros((nx,nx),dtype=float)
for i in range(nx):
    for j in range(nx):
            k_mod[j,i] = math.sqrt(kx[j]**2+kx[i]**2)

spettro_B_sferico_pesi = []
for i in range(tempi+1):
    spettro_B_sferico_pesi.append( np.zeros(nx+1,dtype=float) )
    
for t in range(tempi+1):
    for i in range(0,nx):
        for j in range(0,nx):      
            l = np.argmin(np.abs(k_base-k_mod[j,i]))
            dist = k_base[l]-k_mod[j,i]
            if dist>0 :
                spettro_B_sferico_pesi[t][l] += spettro_B_gyro_pesi_finale[t][j,i]*(1-dist)
                spettro_B_sferico_pesi[t][l-1] += spettro_B_gyro_pesi_finale[t][j,i]*(dist)
            elif dist<0 :
                spettro_B_sferico_pesi[t][l] += spettro_B_gyro_pesi_finale[t][j,i]*(1-np.abs(dist))
                spettro_B_sferico_pesi[t][l+1] += spettro_B_gyro_pesi_finale[t][j,i]*(np.abs(dist))
            else:# (k_base - k[j,i])=0 :
                spettro_B_sferico_pesi[t][l] += spettro_B_gyro_pesi_finale[t][j,i]
                 
spettro_B_sferico_finale = []
for t in range(tempi+1):
    spettro_B_sferico_finale.append(np.delete(spettro_B_sferico_pesi[t],[nx-1,nx])) 
#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.loglog(k_base[0:nx-1],spettro_B_sferico_finale[i])
#    plt.ylim((10**(-9),10**9))
#%% COSì FO LO SPETTRO SFERICO DA ZERO... CI VUOLE DI PIù MA VIENE PRATICAMENTE IDENTICO A QUELLO SOPRA        
#k_sfera = np.zeros((nz,ny,nx),dtype=float)
#for i in range(nx):
#    for j in range(ny):
#        for w in range(nz):
#            k_sfera[w,j,i] = math.sqrt(kz[w]**2+ky[j]**2+kx[i]**2)
#            
#            
#spettro_B_sferico_pesi_1 = []
#for i in range(tempi+1):
#    spettro_B_sferico_pesi_1.append( np.zeros(nx+1,dtype=float) )
#
#for t in range(tempi+1):
#    for w in range(0,nx):
#        for j in range(0,nz):
#            for i in range(0,ny): 
#                l = np.argmin(np.abs((k_base-k_sfera[j,i,w])))
#                dist = k_base[l]-k_sfera[j,i,w]
#                if dist>0 :
#                    spettro_B_sferico_pesi_1[t][l] += BfM2[t][j,i,w]*(1-dist)
#                    spettro_B_sferico_pesi_1[t][l-1] += BfM2[t][j,i,w]*(dist)
#                elif dist<0 :
#                    spettro_B_sferico_pesi_1[t][l] += BfM2[t][j,i,w]*(1-np.abs(dist))
#                    spettro_B_sferico_pesi_1[t][l+1] += BfM2[t][j,i,w]*(np.abs(dist))
#                else:# (k_base - k[j,i])=0 :
#                    spettro_B_sferico_pesi_1[t][l] += BfM2[t][j,i,w]
#
#spettro_B_sferico_pesi_1_finale = []
#for t in range(tempi+1):
#    spettro_B_sferico_pesi_1_finale.append( np.delete(spettro_B_sferico_pesi_1[t],[nx,nx-1]))            
# #%% Bf modulo
#for i in range(tempi+1):
#    a = str(i)
#    plt.figure(i)
#    plt.loglog(k_base[0:nx-1],spettro_B_sferico_pesi_1_finale[i])
#    plt.loglog(k_base[0:nx-1],spettro_B_sferico_finale[i])  
#%%           
            
spettro_B_perp_rid = [] ; spettro_B_parall_rid = []
for i in range(tempi+1):
    spettro_B_perp_rid.append( np.zeros(nx,dtype=float) )
    spettro_B_parall_rid.append( np.zeros(nx,dtype=float) )
    
for t in range(tempi+1):
    for i in range(0,nx):
        for j in range(0,nx):      
                spettro_B_perp_rid[t][j] += spettro_B_gyro_pesi_finale[t][j,i] 
                spettro_B_parall_rid[t][i] += spettro_B_gyro_pesi_finale[t][j,i] 
                
spettro_B_perp_rid_finale = [] ; spettro_B_parall_rid_finale = []
for t in range(tempi+1):
    spettro_B_perp_rid_finale.append(np.delete(spettro_B_perp_rid[t],[nx-1]))     
    spettro_B_parall_rid_finale.append(np.delete(spettro_B_parall_rid[t],[nx-1]))
#%%    
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.loglog(k_base[0:nx-1],spettro_B_perp_rid_finale[i])    
#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.loglog(k_base[0:nx-1],spettro_B_parall_rid_finale[i])     
#%% CONFRONTO
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.loglog(k_base[0:nx-1],spettro_B_perp_rid_finale[i])  
    plt.loglog(k_base[0:nx-1],spettro_B_parall_rid_finale[i])         
    plt.loglog(k_base[0:nx-1],spettro_B_sferico_finale[i])
    plt.title("blue(K_perp) - orange(K_parallel) - green(K_isotropic)")
    plt.ylim((10**-1,10**11))
    