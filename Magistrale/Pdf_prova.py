# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 10:19:43 2019

@author: Mattia
"""

import h5py ; import pylab as plt ; import numpy as np ; import time ; from scipy import stats

bx=[]; by=[]; bz=[]; 
tempi = 8

for i in range(tempi+1):
    if i<10: c = str(i) ; a = '0'+c
    else: a = str(i)
    bx.append( (h5py.File('../out0%s.h5' %a)).get('bx')[...] )
    by.append( (h5py.File('../out0%s.h5' %a)).get('by')[...] )
    bz.append( (h5py.File('../out0%s.h5' %a)).get('bz')[...] )

nx = len( bx[0][0,0,:] ) ; ny = len( bx[0][0,:,0] ) ; nz = len( bx[0][:,0,0] )

#%% 
''' pdf 2D, solo alla scala più piccola (easy easy...) '''
start=time.time()

d_bx = []
for t in range(tempi+1):
    d_bx.append([])

for t in range(tempi+1):
    for i in range(nx-1):    #distanze orizontali , metodo veloce 
        d_bx[t].append( bx[t][0,:,i+1] - bx[t][0,:,i] )  
    d_bx[t].append( bx[t][0,:,0] - bx[t][0,:,nx-1] )
    
#    d_bx[t].append( bx[t][0,0,:] - bx[t][0,ny-1,:] )  #distanze verticali , metodo veloce
#    for j in range(ny-1):
#        d_bx[t].append( bx[t][0,j+1,:] - bx[t][0,j,:] ) 

    d_bx[t] = np.asarray(d_bx[t]).reshape(-1)

    
end=time.time()
print('time = ',end-start)
#%%
for t in range(tempi+1):
    plt.figure(t)
    plt.hist(d_bx[t],bins=100,log=True)
    plt.title('distribuzione delta bx')
#%% 
hist=[] ; base_hist=[]
for t in range(tempi+1):
    hist.append(np.histogram(d_bx[t],bins=100)[0])
    base_hist.append(np.histogram(d_bx[t],bins=100)[1])
    
#%%
mom_bx_py = np.zeros((tempi+1,5),dtype=float)
for t in range(tempi+1):
    for j in range(5):
        mom_bx_py[t,j] = stats.moment(d_bx[t],moment=j)

kurtosis = mom_bx_py[:,4]/((mom_bx_py[:,2])**2) - 3 #; print('kurtosis=',kurtosis)
skewness = mom_bx_py[:,3]/((mom_bx_py[:,2])**1.5)   #; print('skewness=',skewness)

#%%
''' pdf 2D a tutte le scale. solo orizzontale (ancora abbastanza easy)(su matrici troppo grosse non funge... ha bisogno di troppa ram)'''
start=time.time()

d_bx_or = []
for t in range(tempi+1):
    d_bx_or.append([])
    for l in range(1,int(nx/2)+1):
        d_bx_or[t].append([])

for t in range(tempi+1):
    for l in range(1,10):
        for i in range(nx):
            if i+l < nx :
                d_bx_or[t][l-1].append( bx[t][0,:,i+l] - bx[t][0,:,i] )  
            elif i+l >= nx :
                d_bx_or[t][l-1].append( bx[t][0,:,i+l-nx] - bx[t][0,:,i] )  
        d_bx_or[t][l-1] = np.asarray(d_bx_or[t][l-1]).reshape(-1)

#mom = np.zeros((tempi+1,int(nx/2),5),dtype=float)
#for t in range(tempi+1):
#    for l in range(int(nx/2)):
#        for j in range(5):
#            mom[t,l,j] = stats.moment(d_bx_or[t][l],moment=j)
#        mom[t,l,:] = mom[t,l,:]/mom[t,l,0]
#
#K = mom[:,:,4]/((mom[:,:,2])**2) - 3 
#S = mom[:,:,3]/((mom[:,:,2])**1.5)


end=time.time()
print('time = ',end-start)    
#%%
''' pdf 2D a tutte le scale. solo orizzontale (ancora abbastanza easy)(questo ci mette due secoli a girare)'''
start=time.time()

d_bx_or = np.zeros((tempi+1,int(nx/2),int(nx*nx)),dtype=float)

for t in range(tempi+1):
    for l in range(1,int(nx/2)+1):
        a=0
        for i in range(nx):
            if i+l < nx :
                a=np.append( bx[t][0,:,i+l] - bx[t][0,:,i] , a )  
            elif i+l >= nx :
                a=np.append( bx[t][0,:,i+l-nx] - bx[t][0,:,i] , a )  
        a=np.delete(a,[int(nx*nx)])
        d_bx_or[t,l-1,:] = a
#        d_bx_or[t,l-1] = np.asarray(d_bx_or[t][l-1]).reshape(-1)
#
#mom = np.zeros((tempi+1,int(nx/2),5),dtype=float)
#for t in range(tempi+1):
#    for l in range(int(nx/2)):
#        for j in range(5):
#            mom[t,l,j] = stats.moment(d_bx_or[t][l],moment=j)
#        mom[t,l,:] = mom[t,l,:]/mom[t,l,0]
#
#K = mom[:,:,4]/((mom[:,:,2])**2) - 3 
#S = mom[:,:,3]/((mom[:,:,2])**1.5)


end=time.time()
print('time = ',end-start)    


#%%   

   
        
        