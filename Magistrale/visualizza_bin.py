# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 12:13:45 2019

@author: Mattia
"""

import numpy as np ; from scipy.io import FortranFile ; import pylab as plt ; import matplotlib.cm as cm ; import matplotlib as mplt ; import time
import funzioni_matti as matti ; import funzioni_matti_CINECA as matti_C
start=time.time()

bx=[]; by=[]; bz=[]; btot=[]; vx=[]; vy=[]; vz=[]; pe=[]; rh=[] ; w=[]; vtot=[]; gamma=[]; se=[]; ex=[]; ey=[]; ez=[]; ch=[]; chw=[] ; cs=[]; va=[] ; MACs=[] ;MACa=[]
I=[] ; Q=[] ; U=[] ; js=[] ; pezzoB=[] ; D=[] ; qy=[] ; qz=[] ; cos2chi=[] ; sin2chi=[] ; Itf=[] ; Imed=[]
#directory='/Users/Mattia/Desktop/risultati_LDZ/256^3 (init_turb. va=0.6 tmax=40)/'                   ; N=1024 ; n=256  ; nnz=256

#directory='/Users/Mattia/Desktop/risultati_LDZ/2048^2 (init_turb. va=0.5 tmax=36)/'                  ; N=256  ; n=2048 ; nnz=1  ; ti=24 ; tf=24 ; B0=1.0
#directory='/Users/Mattia/Desktop/risultati_LDZ/2048^2 (init_turb. va=0.87 tmax=4.5)/'                ; N=16   ; n=2048 ; nnz=1 ; ti=16 ; tf=16 ; B0=3.0
#directory='/Users/Mattia/Desktop/risultati_LDZ/2048^2 (init_turb. va=0.75 tmax=6)/'                  ; N=256  ; n=2048 ; nnz=1 ; ti=23 ; tf=23 ; B0=2.0
directory='/Users/Mattia/Desktop/risultati_LDZ/2048^2 (init_turb. va=0.6 tmax=25)/'                  ; N=256  ; n=2048 ; nnz=1 ; ti=51 ; tf=51 ; B0=1.3
#directory='/Users/Mattia/Desktop/risultati_LDZ/4096^2 (init_turb. va=0.65 tmax=5)/'                  ; N=256  ; n=4096 ; nnz=1 ; ti=20 ; tf=20 ; B0=1.5

#directory='/Users/Mattia/Desktop/risultati_LDZ/2048^2 (init_turb_rms. va=0.57 tmax=20)/'             ; N=1024 ; n=2048 ; nnz=1 ; ti=85 ; tf=85 ; B0=3.0

#directory='/Users/Mattia/Documents/Tesi_Mag/Echo-LDZ/'                                               ; N=4    ; n=128  ; nnz=1


d=13        #ogni quanti tempi leggo?
ti=0        #tempo iniziale
tf=39      #tempo finale

save=True

tempi = int((tf-ti)/d) 
t_out=np.zeros(tempi+1,dtype=float)

for i in range(tempi+1):

    array=np.zeros((nnz,n,n,12),dtype=np.float32) #ha bisogno di annullarlo ogni volta per motivi poco chiari.....
    
    if d*i+ti<1000: a = str(d*i+ti)
    if d*i+ti<100:  c = str(d*i+ti) ; a = '0'+c
    if d*i+ti<10:   c = str(d*i+ti) ; a = '00'+c
    
#    if d*i+ti<100:  a = str(d*i+ti) 
#    if d*i+ti<10:   c = str(d*i+ti) ; a = '0'+c
    
    print(i,' - out%s.dat' %a)
    f=FortranFile(directory+'out%s.dat' %a)
    
    data0=f.read_reals(dtype=np.float32)
    t_out[i]=data0[0]
     
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

   
#    rh.append( array[:,:,:,0] )
#    vx.append( array[:,:,:,1] )
#    vy.append( array[:,:,:,2] )
#    vz.append( array[:,:,:,3] )
#    pe.append( array[:,:,:,4] )
#    se.append( array[:,:,:,5] )
#    ex.append( array[:,:,:,6] )
#    ey.append( array[:,:,:,7] )
#    ez.append( array[:,:,:,8] )
    bx.append( array[:,:,:,9] )
    by.append( array[:,:,:,10] )
#    bz.append( array[:,:,:,11] )     ############################################### RICORDA IL B0
#    btot.append( np.sqrt(bx[i]**2 + by[i]**2 + bz[i]**2) )
#    vtot.append( np.sqrt(vx[i]**2 + vy[i]**2 + vz[i]**2) )
#    gamma.append( (1.0-vtot[i]**2)**(-1/2) )
    
#    cs.append( np.sqrt( 4*pe[i]/(3*(rh[i]+4*pe[i])) ) )
#    va.append( B0/(np.sqrt(B0**2 + rh[i] + 4*pe[i])) )
#    MACs.append(vtot[i]/cs[i])
#    MACa.append(vtot[i]/va[i])
#    del(rh,pe,cs,va)
 
    
#    ch.append(2*( bx[i]*vx[i] + by[i]*vy[i] + vz[i]*bz[i] )/( bx[i]**2 + by[i]**2 + bz[i]**2 + vx[i]**2 + vy[i]**2 + vz[i]**2 ))

#    w.append( rh[i] + 4*pe[i] + bx[i]*bx[i] + by[i]*by[i] + bz[i]*bz[i] )
##    chw.append( (2/np.sqrt(w[i]))*( bx[i]*vx[i] + by[i]*vy[i] + vz[i]*bz[i] )/( ((bx[i]*bx[i] + by[i]*by[i] + bz[i]*bz[i])/w[i]) + vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i] ) )
#    chw.append( (2/np.sqrt(w[i]))*( bx[i]*vx[i] + by[i]*vy[i] + vz[i]*bz[i] ))
#    norm = matti_C.Media( ((bx[i]*bx[i] + by[i]*by[i] + bz[i]*bz[i])/w[i]) + vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i] ) 
#    chw[i]=chw[i]/norm
#    del(w,rh,pe,bx,by,bz,vx,vy,vz,norm)
    
#    '''BUCCIA'''
#    D.append( 1.0/(gamma[i]*(1-vx[i])))
#    pezzoB.append( (1/gamma[i])*np.sqrt(btot[i]**2 - (bx[i]*D[i])**2 +2*gamma[i]*D[i]*bx[i]*( vx[i]*bx[i] + vy[i]*by[i] + vz[i]*bz[i]) ) )
#    js.append( pe[i]*(pezzoB[i]**(1.5))*D[i]**(2.5) )
#    qy.append( (1-vx[i])*by[i]+vy[i]*bx[i] )
#    qz.append( (1-vx[i])*bz[i]+vz[i]*bx[i] )
#    cos2chi.append( (qy[i]**2-qz[i]**2)/(qy[i]**2+qz[i]**2) )
#    sin2chi.append( -(qy[i]*qz[i])/(qy[i]**2+qz[i]**2) )
#    I.append(np.zeros(len(js[i][0,0,:]),dtype=float))
#    Q.append(np.zeros(len(js[i][0,0,:]),dtype=float))
#    U.append(np.zeros(len(js[i][0,0,:]),dtype=float))
#    for j in range(len(js[i][0,:,0])):
#        I[i][:] +=  js[i][0,:,j] 
#        Q[i][:] +=  js[i][0,:,j]*cos2chi[i][0,:,j]
#        U[i][:] +=  js[i][0,:,j]*sin2chi[i][0,:,j]
#    Imed.append(np.average(I[i]))
#    Itf.append( (np.absolute((np.fft.rfftn(I[i]-Imed[i]))**2) ))
    
    f.close()    
#
#dexdt=[] ; deydt=[] ; dezdt=[]
#for i in range(tempi):
#    dexdt.append( (ex[i+1]-ex[i])/(t_out[i+1]-t_out[i]) )
#    deydt.append( (ey[i+1]-ey[i])/(t_out[i+1]-t_out[i]) )
#    dezdt.append( (ez[i+1]-ez[i])/(t_out[i+1]-t_out[i]) )

#divV=matti.Divergence(vx,vy,vz,tempi)
#rotVx , rotVy , rotVz =matti.Curl_migliorissimo(vx,vy,vz,tempi)
#rotV=[]
#for i in range(tempi+1):
#    rotV.append(np.sqrt(rotVx[i]**2+rotVy[i]**2+rotVz[i]**2))

del(data0,data01,data02,data03,ix1,ix2,iy1,iy2,iz1,iz2,arr,ixmin,ixmax,array,n,p,nx,ny,nz,nv,nnz,N)

end=time.time()
matti.timing(start,end)

#%% MAC
for i in range(tempi+1):
    vmax=np.amax(MACs[i])
    norm = mplt.colors.Normalize(0,vmax)
    if vmax>1 : colors= [[norm(0), "darkblue"],[norm(1.0), "white"],[norm(1.0), "red"],[norm(vmax), "yellow"]]
    elif vmax<1 : colors= [[norm(0), "darkblue"],[norm(vmax), "white"]]
    cmap = mplt.colors.LinearSegmentedColormap.from_list("", colors)
    a = str(i)
    plt.figure(23)
    plt.contourf(MACs[i][0,:,:],255,cmap=cmap)    #levels=np.linspace(-1.0,1.0,255)
    title='MAC Sonoro'
    plt.title(title)
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')
for i in range(tempi+1):
    vmax=np.amax(MACa[i])
    norm = mplt.colors.Normalize(0,vmax)
    if vmax>1 :colors= [[norm(0), "darkblue"],[norm(1.0), "white"],[norm(1.0), "red"],[norm(vmax), "yellow"]]
    elif vmax<1 : colors= [[norm(0), "darkblue"],[norm(vmax), "white"]]
    cmap = mplt.colors.LinearSegmentedColormap.from_list("", colors)
    a = str(i)
    plt.figure(24)
    plt.contourf(MACa[i][0,:,:],255,cmap=cmap)    #levels=np.linspace(-1.0,1.0,255)
    title='MAC Alfvénico'
    plt.title(title)
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')
    
#%% comprimibilita 
for i in range(tempi +1):
    vmax=np.amax(divV[i]) ; vmin=np.amin(divV[i])  ; vrif=max(vmax,np.abs(vmin))
    norm = mplt.colors.Normalize(vmin,vmax)
    colors= [[norm(vmin), "darkblue"],[norm(vmin/2), "blue"],[norm(0.0), "white"],[norm(vmax/2), "red"],[norm(vmax), "darkred"]]
    cmap = mplt.colors.LinearSegmentedColormap.from_list("", colors)
    a = str(i)
    plt.figure(i)
    plt.contourf(divV[i][0,:,:],255,cmap=cmap)    #levels=np.linspace(-1.0,1.0,255)
    title='Comprimibilita'
    plt.title(r' $\nabla\cdot v$')
    plt.colorbar()
#    if (save==True): plt.savefig(directory+'immagini/'+title+'300.png', format='png',dpi=300,bbox_inches='tight')

#%% rotore v
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(rotV[i][0,:,:],100,cmap=cm.seismic_r)    #levels=np.linspace(-1.0,1.0,255)
    title='Rotore V'
    plt.title(title+r'$\nabla\times v$')
    plt.colorbar()
#    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')

#%% gamma
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(gamma[i][0,:,:],255,cmap=cm.YlOrRd_r)    #levels=np.linspace(-1.0,1.0,255)
    title='Lorentz Factor'
    plt.title(title+' $\gamma$')
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')
    
#%% CH
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(chw[i][0,:,:],levels=np.linspace(-1.0,1.0,255),cmap=cm.seismic)    #levels=np.linspace(-1.0,1.0,255)
    title='CrossHelicityDensity_zero'
    plt.title('Relativistic cross helicity density')
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')
    
#%% bx
for i in range(1):
    a = str(i)
    plt.figure(i)
    plt.contourf(btot[i][0,:,:],255,cmap=cm.YlGnBu_r)
    title='Bx iniziale'
    plt.title(title)
    plt.colorbar()
#    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')

#%% by
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(by[i][0,:,:],255,cmap=cm.YlGnBu)
    plt.title('by%s' %a)
    plt.colorbar()
    
#%% bz
''' ATTENTO AL B0   °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'''
for i in range(tempi+1):
    vmax=np.amax(bz[i]) ; vmin=np.amin(bz[i])
    norm = mplt.colors.Normalize(vmin,vmax)
    colors= [[norm(vmin), "darkgreen"],[norm(vmin/2), "green"],[norm(0.0), "white"],[norm(vmax/2), "orange"],[norm(vmax), "#fe0002"]]
    cmap = mplt.colors.LinearSegmentedColormap.from_list("", colors)
    a = str(i)
    plt.figure(i)
#    plt.contourf(bz[i][0,:,:],255,cmap=cm.YlGnBu_r) 
    plt.contourf(bz[i][0,:,:],255,cmap=cmap) 
    title='B_z-col2'
    plt.title(title)
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')
    
#%% b perp
h=0
for i in [2]:
    a = str(i)
    plt.figure(100+i)
    plt.contourf(np.sqrt(bx[i][h,:,:]**2+by[i][h,:,:]**2),levels=np.linspace(0,1.7,255),cmap=cm.YlGnBu_r)
#    plt.xlim(250,750) ; plt.ylim(300,800) 
    title='B_perpendicolare t=6.50'
#    plt.title('Perpendicular magnetic field map')
    plt.title('t=6.50')
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')
    
#%% v x
for i in range(tempi+1):
    a = str(i)
    plt.figure(8)
    plt.plot(vx[i][0,1733,:])#,155,cmap=cm.hot_r)
    title='Vx iniziale'
    plt.title(title)
#    plt.colorbar()
#    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')
    

#%% v z
for i in range(1):
    a = str(i)
    plt.figure(i)
    plt.contourf(vx[i][0,:,:],255,cmap=cm.YlOrRd_r)
    title='V_z'
    plt.title(title)
    plt.colorbar()
#    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')
    
#%% v perp
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(np.sqrt(vx[i][0,:,:]**2+vy[i][0,:,:]**2),255,cmap=cm.YlOrRd_r)
    title='V_perpendicolare'
    plt.title('Perpendicular velocity field map')
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')
    
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
    plt.figure(30)
    plt.contourf(rh[i][0,:,:],155,cmap=cm.hot_r)
#    plt.xlim(1690,1750) ; plt.ylim(990,1050) 
    title='Densità'# particolare'
    plt.title(title)
    plt.colorbar()
#    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')
    
#%% dedt perp
for i in range(tempi):
    a = str(i)
    plt.figure(13)
    plt.contourf(np.sqrt(dexdt[i][0,:,:]**2+deydt[i][0,:,:]**2),255,cmap=cm.Spectral)
    title='deperpdt'
    plt.title(title)
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')
    
#%% dedt z
for i in [0]:
    vmax=np.amax(dezdt[i]) ; vmin=np.amin(dezdt[i]) ; vrif=max(vmax,np.abs(vmin))
    a = str(i)
    plt.figure(14)
    plt.contourf(dezdt[i][0,:,:],levels=np.linspace(-vrif,vrif,255),cmap=cm.seismic_r)
    title='dezdt'
    plt.title(r'$\partial E_z/\partial t$')
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')
    
#%%
# salvare un immagine vettoriale con python
plt.savefig('/Users/Mattia/Desktop/risultati LDZ/4096^2 (init_turb. va=0.65 tmax=5.25)/nome.eps', format='eps',dpi=300,bbox_inches='tight')
#%%
'''SPETTRO B 3D ################################################################################################################'''
start=time.time()
BfM2=[]
for i in range(tempi+1):
    BfM2.append( (np.absolute(np.fft.rfftn(np.transpose(bx[i],(2,1,0)))))**2 + (np.absolute(np.fft.rfftn(np.transpose(by[i],(2,1,0)))))**2 + (np.absolute(np.fft.rfftn(np.transpose(bz[i],(2,1,0)))))**2 )
end=time.time()
matti.timing(start,end)
#%%
spettro_B_gyro , kparal, kperp = matti.gyrotropic_spectrum_3D(BfM2,tempi)
#%%
nx=len(spettro_B_gyro[0][0,:])
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(spettro_B_gyro[i],cmap=cm.jet,levels=np.logspace(-1, 14, 255),locator=mplt.ticker.LogLocator())
    mplt.pyplot.yscale('log') ; mplt.pyplot.xscale('log') ; plt.ylim((1,nx-1)) ; plt.xlim((1,nx-1))
    plt.ylabel('k parallel') ; plt.xlabel('k perpendicular') ; plt.colorbar() ; plt.title('Spettro B gyrotropico (%s)' %a)
#%%
spettro_B_iso, spettro_B_perp, spettro_B_parall = matti.spettri_ridotti(BfM2,spettro_B_gyro,tempi)
#%%
alpha=5/3
nx=len(spettro_B_perp[0]) 

kolmogorov=np.zeros(nx,dtype=float)
for i in range(nx):
    kolmogorov[i]=4*10**13

for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.loglog(kperp[0:nx],kolmogorov)
    plt.loglog(kperp[0:nx],spettro_B_perp[i]*(kperp[0:nx])**(alpha),label='Spettro perpendicolare')
    plt.loglog(kparal[0:nx],spettro_B_parall[i]*(kparal[0:nx])**(alpha),label='Spettro parallelo')
    plt.loglog(kperp[0:nx],spettro_B_iso[i]*(kperp[0:nx])**(alpha),label='Spettro isotropico')
    mplt.pyplot.legend()
    plt.title('Spettro gyrotropic B, alpha=%f'%alpha)
    plt.ylim((10**7,10**15))

    
#%%
'''LINEE CAMPO B 2D ################################################################################################################'''    
Az_x , Az_y , Az = matti.linee_campo_B_2D(bx,by,tempi) 
#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(1000)
    plt.contour(Az[i][0,:,:],255,cmap=cm.Spectral_r)
    plt.title('Az %s' %a)
    plt.colorbar()

#%%
'''RMS ################################################################################################################'''
rms_bx = matti.Varianza(bx,tempi) ; rms_by = matti.Varianza(by,tempi) ; rms_bz = matti.Varianza(bz,tempi) ; rms_B = rms_by + rms_bz + rms_bx
rms_vx = matti.Varianza(vx,tempi) ; rms_vy = matti.Varianza(vy,tempi) ; rms_vz = matti.Varianza(vz,tempi) ; rms_V = rms_vy + rms_vz + rms_vx
#%%
'''SPETTRO B 2D ################################################################################################################'''
BfM2=[]
for i in range(tempi+1):
    BfM2.append( (np.absolute(np.fft.rfftn(bx[i])))**2 + (np.absolute(np.fft.rfftn(by[i])))**2 + (np.absolute(np.fft.rfftn(bz[i])))**2 )
spettro_b , k = matti.radial_spectrum_2D(BfM2,tempi)
#%% visualizza BfM2   
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(BfM2[i][0,:,:],levels=np.logspace(-2,7,25),cmap=cm.Spectral_r,locator=mplt.ticker.LogLocator())
    plt.title('BfM2 (%s)' %a)
    plt.colorbar()
#%% visualizza spettro 2D
alpha=1
kolmogorov=np.zeros(len(spettro_b[0]),dtype=float)
for i in range(len(spettro_b[0])):
    kolmogorov[i]=5*10**14
    
plt.figure(789)
#plt.loglog(k[0:len(spettro_b[0])],kolmogorov)
for i in range(tempi+1):
    a = str(i)
    #plt.figure(i)
    plt.loglog(k[0:len(spettro_b[0])],spettro_b[i]*k[0:len(spettro_b[0])]**(alpha),label='%s'%a)
    plt.title('spettro 2D B, alpha=%f'%alpha)
plt.legend()
if (save==True): plt.savefig('/Users/Mattia/Desktop/risultati LDZ/4096^2 (init_turb. va=0.65 tmax=5)/immagini/spettroB2.png', format='png',dpi=300,bbox_inches='tight')

#%%
'''SPETTRO V 2D ################################################################################################################'''
VfM2=[]
for i in range(tempi+1):
    VfM2.append( (np.absolute(np.fft.rfftn(vx[i])))**2 + (np.absolute(np.fft.rfftn(vy[i])))**2 + (np.absolute(np.fft.rfftn(vz[i])))**2 )
spettro_v , k = matti.radial_spectrum_2D(VfM2,tempi)
#%% visualizza VfM2   
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(VfM2[i][0,:,:],levels=np.logspace(-2,7,25),cmap=cm.Spectral_r,locator=mplt.ticker.LogLocator())
    plt.title('VfM2 (%s)' %a)
    plt.colorbar()
#%% visualizza spettro 2D
alpha=0
kolmogorov=np.zeros(len(spettro_v[0]),dtype=float)
for i in range(len(spettro_v[0])):
    kolmogorov[i]=5*10**13#*k[i]**(-5/3)
    
plt.figure(456)
plt.loglog(k[0:len(spettro_v[0])],kolmogorov)
for i in range(tempi+1):
    a = str(i)
    #plt.figure(i)
    plt.loglog(k[0:len(spettro_v[0])],spettro_v[i]*k[0:len(spettro_v[0])]**(alpha),label='%s'%a)
    plt.title('spettro 2D V, alpha=%f'%alpha)
plt.legend()
if (save==True): plt.savefig('/Users/Mattia/Desktop/risultati LDZ/4096^2 (init_turb. va=0.65 tmax=5)/immagini/spettroV.png', format='png',dpi=300,bbox_inches='tight')

#%% 
''' SPETTRO COMPRESSIBLE E SOLENOIDAL'''
spettro_v_com , spettro_v_sol , k = matti.radial_spectrum_2D_compressible_solenoidal(vx,vy,vz,tempi)
#%%
na=len(spettro_v_com[0])
alpha=2.2 ; beta=5/3
#alpha=0.0 ; beta=0.0
rettacom=np.zeros(na,dtype=float) ; rettasol=np.zeros(na,dtype=float)
for i in range(na):
    rettacom[i]=6*10**11
    rettasol[i]=8*10**10

wid=1.1

plt.figure(251)
plt.ylim(1*10**10,9*10**11)
#plt.ylim(1*10**9,2*10**11)
#plt.loglog(k[0:na],rettacom,'xkcd:black',linewidth=0.8)
#plt.loglog(k[0:na],rettasol,'xkcd:black',linewidth=0.8)
for i in range(tempi+1):
    a = str(i)
    #plt.figure(i)
    plt.loglog(k[0:na],spettro_v[i]*k[0:na]**(beta),'xkcd:green',linewidth=wid,label=r'$\mathcal{P}_v k^{5/3}$')
    plt.loglog(k[0:na],spettro_v_sol[i]*k[0:na]**(beta),'xkcd:brown',linewidth=wid,label=r'$\mathcal{P}_v^s k^{5/3}$')
    plt.loglog(k[0:na],spettro_v_com[i]*k[0:na]**(alpha),'xkcd:purple',linewidth=wid,label=r'$\mathcal{P}_v^c k^{2.2}$')
title='spettro 2D V, compr & sol, alpha=%1.2f, beta=%1.2f'%(alpha,beta)
plt.title('Compressible and solenoidal compensated spectra',fontsize='x-large')
plt.legend(loc='best',fontsize='x-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% distribuzione delle velocità
'''DISTRIBUZIONE V ##########'''
vtot_arr=[]
for i in range(tempi+1):
    if d*i+ti<10: c = str(d*i+ti) ; a = '0'+c
    else: a = str(d*i+ti)
    plt.figure(999)
    dis_v=np.histogram(vtot[i],100)
    plt.plot(dis_v[1][1:len(dis_v[1])],dis_v[0],label='%s'%a)
    plt.legend()
    vtot_arr.append(np.reshape(vtot[i],-1))

#%%    
mom, K, S = matti.momenti_distribuzione(vtot_arr,tempi)
plt.figure(777) ; plt.plot(t_out,K) ; plt.title('Kurtosis') ; plt.figure(888) ; plt.plot(t_out,S) ; plt.title('Skewness')

#%%
'''DENSITà DI CORRENTE ################################################################################################################'''
Jx , Jy , Jz = matti.Curl_migliorissimo(bx,by,bz,tempi)
#max_Jx = [] ; max_Jy =[] ; max_Jz =[]
#for t in range(tempi+1):
#    max_Jx.append(np.amax(Jx[t])) ; max_Jy.append(np.amax(Jy[t])) ; max_Jz.append(np.amax(Jz[t]))
#Jx=[] ;Jy=[]; Jz=[]
#
#for i in range(tempi+1):
#    Jx.append(np.zeros(1,dtype=np.float32)) ; Jy.append(np.zeros(1,dtype=np.float32)) ; Jz.append(np.zeros(1,dtype=np.float32))
#    Jx[i],Jy[i],Jz[i]=matti_C.Curl_migliorissimo(bx[i],by[i],bz[i])
    
#%% J perp
for i in range(tempi+1):
    a = str(i)
    plt.figure(i) 
    plt.contourf(np.sqrt(Jx[i][0,:,:]**2+Jy[i][0,:,:]**2),255,cmap=cm.Spectral)
    title='J_perpendicolare'
    plt.title(title)
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')
   
#%% J z
for i in [0]:
    vmax=np.amax(Jz[i]) ; vmin=np.amin(Jz[i]) ; vrif=max(vmax,np.abs(vmin))
#    norm = mplt.colors.Normalize(vmin,vmax)
#    colors= [[norm(vmin), "darkred"],[norm(vmin/2), "red"],[norm(0), "white"],[norm(vmax/2), "blue"],[norm(vmax), "darkblue"]]
#    cmap = mplt.colors.LinearSegmentedColormap.from_list("", colors)
    a = str(i)
    plt.figure(i)#,[10,3]) 
    plt.contourf(Jz[i][0,:,:],levels=np.linspace(-vrif,vrif,255),cmap=cm.seismic_r)
#    plt.xlim(0,600) ; plt.ylim(1650,1850) 
#    plt.xlim(50,650) ; plt.ylim(1050,1250)
    title='Jz'# particolare 2'
    plt.title(r'$(\nabla\times B)_z$')
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'500.png', format='png',dpi=500,bbox_inches='tight')
    

 #%%
plt.figure(100) ; plt.plot(t_out,np.sqrt(max_Jx)) ; plt.title('massimi Jx') 
plt.figure(101) ; plt.plot(t_out,np.sqrt(max_Jy)) ; plt.title('massimi Jy') 
plt.figure(102) ; plt.plot(t_out,np.sqrt(max_Jz)) ; plt.title('massimi Jz')    
    
#%%
'''RMS J #####'''
Var_Jx = matti.Varianza(Jx,tempi)
Var_Jy = matti.Varianza(Jy,tempi)
Var_Jz = matti.Varianza(Jz,tempi)
Var_Jtot = Var_Jx + Var_Jy + Var_Jz
#%%
plt.figure(200) ; plt.plot(t_out,np.sqrt(Var_Jx)) ; plt.title('Varianza Jx') 
plt.plot(t_out,np.sqrt(Var_Jy)) ; plt.title('Varianza Jy') 
plt.plot(t_out,np.sqrt(Var_Jz)) ; plt.title('Varianza Jz') 
plt.plot(t_out,np.sqrt(Var_Jtot)) ; plt.title('Varianza J totale')

#%%
'''CONFRONTO rotB dez/dt #####'''
confronto=Jz[0]/(Jz[0]+dezdt[0])
#vmax=np.amax(confronto) ; vmin=np.amin(confronto)
#tot=np.argmax(confronto)
#maxy=int(tot/2048)
#maxx=tot-maxy*2048
#tot=np.argmin(confronto)
#miny=int(tot/2048)
#minx=tot-miny*2048

plt.figure(136)
plt.contourf(confronto[0,:,:],levels=np.linspace(0,1.5,50),cmap=cm.gist_stern_r)
title='confronto RotB dezdt'
plt.title(title)
plt.colorbar()
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')

#%%
''' densità HELICITY J e B ################################################################################################################'''
Jx , Jy , Jz = matti.Curl_migliorissimo(bx,by,bz,tempi)
Ax , Ay , Az = matti.Potenziale_Vettore(Jx,Jy,Jz,tempi)

Helj=[] ; Helb=[]
for i in range(tempi+1):
    Helj.append( bx[i]*Jx[i] + by[i]*Jy[i] + bz[i]*Jz[i])
    Helb.append( bx[i]*Ax[i] + by[i]*Ay[i] + bz[i]*Az[i])

#%% Helicity J
for i in range(tempi+1):
    vmax=np.amax(Helj[i]) ; vmin=np.amin(Helj[i]) ; vrif=max(vmax,np.abs(vmin))
    norm = mplt.colors.Normalize(vmin,vmax)
    colors= [[norm(vmin), "darkgreen"],[norm(vmin/2), "green"],[norm(0.0), "white"],[norm(vmax/2), "purple"],[norm(vmax), "#4b006e"]]
#    colors= [[norm(vmin), "darkred"],[norm(vmin/2), "red"],[norm(0.0), "white"],[norm(vmax/2), "green"],[norm(vmax), "darkgreen"]]
    cmap = mplt.colors.LinearSegmentedColormap.from_list("", colors)
    a = str(i)
    plt.figure(i) 
    plt.contourf(Helj[i][0,:,:],255,cmap=cmap)
    title='Helicity-dens(J) giusto'
    plt.title('Current helicity density')
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'500.png', format='png',dpi=500,bbox_inches='tight')

#%% Helicity B
for i in range(tempi+1):
    vmax=np.amax(Helb[i]) ; vmin=np.amin(Helb[i]) ; vrif=max(vmax,np.abs(vmin))
    norm = mplt.colors.Normalize(vmin,vmax)
    colors= [[norm(vmin), "darkgreen"],[norm(vmin/2), "green"],[norm(0.0), "white"],[norm(vmax/2), "purple"],[norm(vmax), "#4b006e"]]
    cmap = mplt.colors.LinearSegmentedColormap.from_list("", colors)
    a = str(i)
    plt.figure(i) 
    plt.contourf(Helb[i][0,:,:],155,cmap=cmap)
    title='Helicity-dens(B) giusto'
    plt.title('Magnetic helicity density')
    plt.colorbar()
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')


#%% 
''' STIMA dE/dt ################################################################################################################'''
start=time.time()
dexdt=[] ; deydt=[] ; dezdt=[]
for i in range(tempi):
    dexdt.append((ex[i+1]-ex[i])/(t_out[i+1]-t_out[i]))
    deydt.append((ey[i+1]-ey[i])/(t_out[i+1]-t_out[i]))
    dezdt.append((ez[i+1]-ez[i])/(t_out[i+1]-t_out[i]))
    end=time.time()
matti.timing(start,end)
#%%
Var_dexdt = matti.Varianza(dexdt,tempi-1)
Var_deydt = matti.Varianza(deydt,tempi-1)
Var_dezdt = matti.Varianza(dezdt,tempi-1)
Var_detotdt = Var_dexdt + Var_deydt + Var_dezdt

#%%
plt.figure(200) 
plt.plot(t_out[0:tempi],np.sqrt(Var_dexdt),'xkcd:green',linewidth=0.8,label='rms_dexdt')
plt.plot(t_out[0:tempi],np.sqrt(Var_deydt),'xkcd:red',linewidth=0.8,label='rms_deydt')
plt.plot(t_out[0:tempi],np.sqrt(Var_dezdt),'xkcd:blue',linewidth=0.8,label='rms_dezdt')
plt.plot(t_out[0:tempi],np.sqrt(Var_detotdt),'xkcd:black',linewidth=0.8,label='rms_detotdt')
title='rms dEdt'
plt.legend()
plt.title(title)
if (save==True): plt.savefig(directory+'/immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%%
''' SPETTRI bz bperp Btot, vz vperp Vtot, pe rho ################################################################################################################'''
startino=time.time()
bzfM2=[] ; bperpfM2=[] ; vzfM2=[] ; vperpfM2=[] ; rhfM2=[] ; pefM2=[] ; BfM2=[] ; VfM2=[] ; b0=1
for i in range(tempi+1):
    bzfM2.append( (np.absolute(np.fft.rfftn(bz[i]-b0)))**2 )
    bperpfM2.append( (np.absolute(np.fft.rfftn(bx[i])))**2 + (np.absolute(np.fft.rfftn(by[i])))**2 )
    BfM2.append( (np.absolute(np.fft.rfftn(bx[i])))**2 + (np.absolute(np.fft.rfftn(by[i])))**2 + (np.absolute(np.fft.rfftn(bz[i]-b0)))**2 )
    vzfM2.append( (np.absolute(np.fft.rfftn(vz[i])))**2 )
    vperpfM2.append( (np.absolute(np.fft.rfftn(vx[i])))**2 + (np.absolute(np.fft.rfftn(vx[i])))**2 )
    VfM2.append( (np.absolute(np.fft.rfftn(vx[i])))**2 + (np.absolute(np.fft.rfftn(vy[i])))**2 + (np.absolute(np.fft.rfftn(vz[i])))**2 )
    rhfM2.append( (np.absolute(np.fft.rfftn(rh[i])))**2 )
    pefM2.append( (np.absolute(np.fft.rfftn(pe[i])))**2 )
    
spettro_bz , k = matti.radial_spectrum_2D(bzfM2,tempi)
spettro_bperp , k = matti.radial_spectrum_2D(bperpfM2,tempi)
spettro_B , k = matti.radial_spectrum_2D(BfM2,tempi)
spettro_vz , k = matti.radial_spectrum_2D(vzfM2,tempi)
spettro_vperp , k = matti.radial_spectrum_2D(vperpfM2,tempi)
spettro_V , k = matti.radial_spectrum_2D(VfM2,tempi)
spettro_rh , k = matti.radial_spectrum_2D(rhfM2,tempi)
spettro_pe , k = matti.radial_spectrum_2D(pefM2,tempi)
spettro_Etot=[] ; spettro_residuale=[]
for i in range(tempi+1):
    spettro_Etot.append(spettro_B[i] + spettro_V[i])
    spettro_residuale.append(spettro_B[i] - spettro_V[i])
endino=time.time()
matti.timing(startino,endino)

#%%
alpha=0
na=len(spettro_bz[0])
plt.figure(1001)
plt.ylim((1*10**5,8*10**11))
for i in range(tempi+1):
    plt.loglog(k[0:na],spettro_bz[i]*k[0:na]**(alpha),'xkcd:blue',linestyle='dashed',linewidth=0.8,label='bz - b0')
    plt.loglog(k[0:na],spettro_bperp[i]*k[0:na]**(alpha),'xkcd:blue',linestyle='dotted',linewidth=0.8,label='bperp')
    plt.loglog(k[0:na],spettro_B[i]*k[0:na]**(alpha),'xkcd:blue',linewidth=0.8,label='Btot - b0')
    plt.loglog(k[0:na],spettro_vz[i]*k[0:na]**(alpha),'xkcd:red',linestyle='dashed',linewidth=0.8,label='vz')
    plt.loglog(k[0:na],spettro_vperp[i]*k[0:na]**(alpha),'xkcd:red',linestyle='dotted',linewidth=0.8,label='vperp')
    plt.loglog(k[0:na],spettro_V[i]*k[0:na]**(alpha),'xkcd:red',linewidth=0.8,label='V tot')
    plt.loglog(k[0:na],spettro_rh[i]*k[0:na]**(alpha),'xkcd:green',linewidth=0.8,label='Densità')
    plt.loglog(k[0:na],spettro_pe[i]*k[0:na]**(alpha),'xkcd:purple',linewidth=0.8,label='Pressione')
    title='Spettri'
    plt.title(title)
plt.legend()
#if (save==True): plt.savefig(directory+'/immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%%
alpha=0
na=len(spettro_bz[0])
plt.figure(509)
plt.ylim((4*10**7,3*10**11))
for i in range(tempi+1):
    plt.loglog(k[0:na],spettro_Etot[i]*k[0:na]**(alpha),'xkcd:blue',linewidth=0.8,label='E tot')
    plt.loglog(k[0:na],spettro_residuale[i]*k[0:na]**(alpha),'xkcd:red',linewidth=0.8,label='E residuale')
    plt.loglog(k[0:na],spettro_B[i]*k[0:na]**(alpha),'xkcd:purple',linewidth=0.8,label='Btot - b0')
    title='Spettri Energia'
    plt.title(title)
plt.legend()
if (save==True): plt.savefig(directory+'/immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%%
alpha=1.3
kolmogorov=np.zeros(len(spettro_rh[0]),dtype=float)
for i in range(len(spettro_rh[0])):
    kolmogorov[i]=2*10**11#*k[i]**(-5/3)
    
plt.figure(712)
plt.loglog(k[0:len(spettro_bz[0])],kolmogorov,'xkcd:black',linewidth=0.8)
plt.ylim((1*10**10,1*10**12))
for i in range(tempi+1):
    plt.loglog(k[0:len(spettro_bz[0])],spettro_bz[i]*k[0:len(spettro_bz[0])]**(alpha),'xkcd:purple',linewidth=0.8,label='24')
    title='Spettro Radiale bz-b0, alpha=%1.2f'%alpha
    plt.title(title)
plt.legend()
if (save==True): plt.savefig(directory+'/immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

alpha=1.4
kolmogorov=np.zeros(len(spettro_rh[0]),dtype=float)
for i in range(len(spettro_rh[0])):
    kolmogorov[i]=2*10**11#*k[i]**(-5/3)

plt.figure(358)
plt.loglog(k[0:len(spettro_rh[0])],kolmogorov,'xkcd:black',linewidth=0.8)
plt.ylim((1*10**10,1*10**12))
for i in range(tempi+1):
    plt.loglog(k[0:len(spettro_rh[0])],spettro_rh[i]*k[0:len(spettro_rh[0])]**(alpha),'xkcd:green',linewidth=0.8,label='24')
    title='Spettro Radiale rho, alpha=%1.2f'%alpha
    plt.title(title)
plt.legend()
if (save==True): plt.savefig(directory+'/immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%%


        
#%%
for i in range(4,5):
#    vmax=np.amax(Helb[i]) ; vmin=np.amin(Helb[i]) ; vrif=max(vmax,np.abs(vmin))
#    norm = mplt.colors.Normalize(vmin,vmax)
#    colors= [[norm(vmin), "#4b006e"],[norm(vmin/2), "purple"],[norm(0.0), "white"],[norm(vmax/2), "green"],[norm(vmax), "darkgreen"]]
#    cmap = mplt.colors.LinearSegmentedColormap.from_list("", colors)
    a = str(i)
    plt.figure(5555) 
    plt.contourf(js[i][0,:,:],255,cmap=cm.Spectral)
    title='Emissività'
    plt.title(title)
    plt.colorbar()
#    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=1000,bbox_inches='tight')
    
#%%
plt.figure(952)
for t in range(tempi+1):
    plt.plot(I[t],linewidth=0.8,label='I')
    plt.plot(np.sqrt(U[t]**2+Q[t]**2),linewidth=0.8,label='Ipol')
    title='Intensità Sincrotrone Integrata'
    plt.title(title)
    plt.legend()
plt.figure(654)
for t in range(tempi+1):
    plt.plot(np.sqrt(U[t]**2+Q[t]**2)/I[t],linewidth=0.8,label='Ipol/I')
    title='Rapporto Polarizzazione'
    plt.title(title)
    plt.legend()
    
#%%
plt.figure(2345)
for t in range(tempi+1):
    a= str(t)
    plt.loglog(Itf[t],linewidth=0.8,label='%s'%a)
    title='Spettro Intensità'
    plt.title(title)
    plt.legend()  
    