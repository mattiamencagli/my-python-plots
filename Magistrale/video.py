# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 11:46:39 2019

@author: Mattia
"""

import numpy as np ; from scipy.io import FortranFile ; import pylab as plt ; import matplotlib.cm as cm #; import matplotlib.animation as an 
import funzioni_matti_CINECA as matti; import time ; import matplotlib as mplt
import matplotlib as mpl
mpl.use('Agg')

start=time.time()

#bx=[]; by=[]; bz=[]; bperp=[]; vx=[]; vy=[]; vz=[]; pe=[]; rh=[] 

#directory='/Users/Mattia/Desktop/risultati_LDZ/256^3 (init_turb. va=0.6 tmax=40)/'                   ; N=1024 ; n=256  ; nnz=256
#directory='/Users/Mattia/Desktop/risultati_LDZ/2048^2 (init_turb. va=0.5 tmax=36)/'                  ; N=256  ; n=2048 ; nnz=1
#directory='/Users/Mattia/Desktop/risultati_LDZ/2048^2 (init_turb. va=0.87 tmax=4.5)/'                ; N=16   ; n=2048 ; nnz=1
#directory='/Users/Mattia/Desktop/risultati_LDZ/2048^2 (init_turb. va=0.75 tmax=6)/'                  ; N=256  ; n=2048 ; nnz=1
directory='/Users/Mattia/Desktop/risultati_LDZ/2048^2 (init_turb. va=0.6 tmax=25)/'                  ; N=256  ; n=2048 ; nnz=1
#directory='/Users/Mattia/Desktop/risultati_LDZ/4096^2 (init_turb. va=0.65 tmax=5)/'                  ; N=256  ; n=4096 ; nnz=1
#directory='/Users/Mattia/Desktop/risultati_LDZ/1024^2 (init_turb. va=0.87 tmax=2)/'                  ; N=256  ; n=1024 ; nnz=1
#directory='/Users/Mattia/Desktop/risultati_LDZ/2048^2 (init_turb_rms. va=0.57 tmax=20)/'             ; N=1024 ; n=2048 ; nnz=1
#directory='/Users/Mattia/Documents/Tesi_Mag/Echo-LDZ/'                                               ; N=4    ; n=128  ; nnz=1



d=1        #ogni quanti tempi leggo?
ti=1        #tempo iniziale
tf=99      #tempo finale

save=True

tempi = int((tf-ti)/d) 
t_out=np.zeros(tempi+1,dtype=float)

bx=np.zeros(tempi+1,dtype=np.float32)
by=np.zeros(tempi+1,dtype=np.float32)
bz=np.zeros(tempi+1,dtype=np.float32)
Jx=np.zeros(tempi+1,dtype=np.float32) 
Jy=np.zeros(tempi+1,dtype=np.float32)
Jz=np.zeros(tempi+1,dtype=np.float32)

for i in range(tempi+1):
    array=np.zeros((nnz,n,n,12),dtype=np.float32) #ha bisogno di annullarlo ogni volta per motivi poco chiari.....
    
    if d*i+ti<1000: a = str(d*i+ti)
    if d*i+ti<100:  c = str(d*i+ti) ; a = '0'+c
    if d*i+ti<10:   c = str(d*i+ti) ; a = '00'+c
    
#    if d*i+ti<100:  a = str(d*i+ti) 
#    if d*i+ti<10:   c = str(d*i+ti) ; a = '0'+c
    
    print(' --- out%s.dat' %a)
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
    vx = ( array[:,:,:,1] )
    vy = ( array[:,:,:,2] )
    vz = ( array[:,:,:,3] )
#    pe.append( array[:,:,:,4] )
#    se.append( array[:,:,:,5] )
#    ex.append( array[:,:,:,6] )
#    ey.append( array[:,:,:,7] )
#    ez.append( array[:,:,:,8] )
#    bx = ( array[:,:,:,9] )
#    by = ( array[:,:,:,10] )
#    bz = ( array[:,:,:,11] )

    divV = matti.Divergence(vx,vy,vz)
    
    vmax=np.amax(divV) ; vmin=np.amin(divV) ; vrif=max(vmax,np.abs(vmin))
    
    a = str(d*i+ti)
    print('figure(%s)'%a)
    fig=plt.figure(i)
    ax=fig.add_axes([0.1,0.1,0.81,0.81])
    
    norm = mplt.colors.Normalize(vmin,vmax)
    colors= [[norm(vmin), "darkblue"],[norm(vmin/2), "blue"],[norm(0.0), "white"],[norm(vmax/2), "red"],[norm(vmax), "darkred"]]
    cmap = mplt.colors.LinearSegmentedColormap.from_list("", colors)
    
    plt.contourf(divV[0,:,:],255,cmap=cmap)
    title='Comprimibilità'
    plt.title(title+r'  $\nabla\cdot v$')
    plt.colorbar()#ticks=[-vrif,-vrif2,-vrif3,0,vrif3,vrif2,vrif])
#    cbar.set_ticklabels([mn,mn2,mn3,'   0',mx3,mx2,mx])
    if (save==True): plt.savefig(directory+'immagini/video/'+title+a+'.png', format='png',dpi=500)#,bbox_inches='tight')
    plt.close()
    

    end1=time.time()
    print('!!!!!!!!!!!!!! Fine ciclo numero  : - ',d*i+ti)
    matti.timing(start,end1)
    f.close()

#del(data0,data01,data02,data03,ix1,ix2,iy1,iy2,iz1,iz2,arr,ixmin,ixmax,array,n,p,nx,ny,nz,nv,nnz,N)

end=time.time()
matti.timing(start,end)

##%%
#ticks=np.zeros(9,dtype=int)
#space = int(len(bx[0][0,0,:])/8)
#for i in range(1,9):
#    ticks[i] = ticks[i-1]+space
#    
#fig = plt.figure()
#ax  = plt.axes(xticks=(ticks), yticks=(ticks))
#title='B_perp'
#plt.title(title)
#plt.ioff()
#vmax= np.amax(bperp[:])
#vmin= np.amin(bperp[:])
#
##colorbar = plt.colorbar(ticks=np.linspace(vmin,vmax,vmax-vmin+1))
#
#def animate(i):
#    imm = plt.contourf(bperp[i][0,:,:],levels=np.linspace(vmin,vmax,255),cmap=cm.YlGnBu_r)
#    return imm
#    
#anim=an.FuncAnimation(fig, animate,frames=tempi+1, interval=1000, repeat=False, blit=False)
#anim.save(directory+title+'.mp4', dpi=1000, writer= an.FFMpegWriter())



