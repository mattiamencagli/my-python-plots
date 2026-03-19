# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 10:57:46 2019

@author: Mattia
"""
import numpy as np ; import math ; from numba import jit ; import time ; from scipy import stats
'''###############################################################################'''
def timing(start,end):
    
    t=end-start
    ms = str(int(t*1000))
    h = str(int(t/3600)) ; t=t-int(h)*3600
    m = str(int(t/60)) ; t=t-int(m)*60
    s = str(int(t))  
    
    if int(h)<10: h='0'+h  
    else: h=h 
    
    if int(m)<10: m='0'+m
    else: m=m
    
    if int(s)<10: s='0'+s  
    else: s=s 
    
    if int(ms)<10: ms='00'+ms
    elif int(ms)<100: ms='0'+ms
    elif int(ms)<1000: ms=ms
    else : ms = ms[len(ms)-3]+ms[len(ms)-2]+ms[len(ms)-1]
    
    print('format :  h  m  s  ms')
    print('time  =  {}:{}:{}:{}'.format(h,m,s,ms))
    print()

'''###############################################################################'''
@jit(cache=True) 
def MaxMinMatrix(f_in,tempi):
    
    start=time.time()
    
    nz=len(f_in[0][:,0,0]) ; ny=len(f_in[0][0,:,0]) ; nx=len(f_in[0][0,0,:])
    
    maxx=[] ; minx=[] ; maxy=[] ; miny=[] ; maxz=[] ; minz=[]
    
    for i in range(tempi*1):
        tot=np.argmax(f_in[i])
        maxz.append(int(tot/(ny*nz)))
        maxy.append(int((tot-maxz[i]*nx*ny)/ny))
        maxx.append(tot-maxz[i]*nx*ny-maxy[i]*ny)
        
        tot=np.argmin(f_in[i])
        minz.append(int(tot/(ny*nz)))
        miny.append(int((tot-minz[i]*nx*ny)/ny))
        minx.append(tot-minz[i]*nx*ny-miny[i]*ny)

    end=time.time()
    print()
    print('END - MaxMinMatrix')
    timing(start,end)
        
    return maxz , maxy , maxx , minz , miny , minx

'''###############################################################################'''
@jit(cache=True)
def Media(f_in,tempi):
    
    start=time.time()
    nz=len(f_in[0][:,0,0]) ; ny=len(f_in[0][0,:,0]) ; nx=len(f_in[0][0,0,:])
    punti = nx*ny*nz
    
    ava_f=np.zeros(tempi+1,dtype=float)
    for t in range(tempi+1):
        for w in range(nz):
            for j in range(ny):
                for i in range(nx):
                    ava_f[t] += f_in[t][w,j,i]
        ava_f[t] = ava_f[t]/punti 
    
    end=time.time()
    print()
    print('END - Media')
    timing(start,end)
        
    return ava_f



'''###############################################################################'''
@jit(cache=True)
def radial_spectrum_primi_vicini_2D(f_in,tempi):
    start=time.time()
    
    ny = len( f_in[0][0,:,0] ) ; nx = len( f_in[0][0,0,:] ) 
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) 
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
    for i in range(nx,ny):
        ky[i]=-ny+i   
    
    spettro_f_finale = [] ; k_base=np.zeros(nx,dtype=float)
        
    k = np.zeros((ny,nx),dtype=float)
    for j in range(ny):
        for i in range(nx):
            k[j,i] = math.sqrt(kx[i]**2+ky[j]**2)
        
    for i in range(0,nx):
        k_base[i] =  k[0,i]
    
    spettro_f = []
    for i in range(tempi+1):
        spettro_f.append( np.zeros(nx,dtype=float) )
    
    for t in range(tempi+1):
        for i in range(0,nx):
            for j in range(0,ny):      
                l = np.argmin(np.abs((k_base-k[j,i])))
                spettro_f[t][l] += f_in[t][0,j,i]
                  
    for t in range(tempi+1):
        spettro_f_finale.append(np.delete(spettro_f[t],[nx-1]))      #l'ultimo bin lo tolgo perchè...
     
    end =time.time()
    print()
    print('END - radial_spectrum_primi_vicini_2D')
    timing(start,end)
        
    return spettro_f_finale , k_base    

'''###############################################################################'''
@jit(cache=True)
def radial_spectrum_2D(f_in,tempi):
    start=time.time()
    
    ny = len( f_in[0][0,:,0] ) ; nx = len( f_in[0][0,0,:] ) 
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) 
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
    for i in range(nx,ny):
        ky[i]=-ny+i   
    
    spettro_f_finale = [] ; k_base=np.zeros(nx,dtype=float)
    
    k = np.zeros((ny,nx),dtype=float)
    for j in range(ny):
        for i in range(nx):
            k[j,i] = math.sqrt(kx[i]**2+ky[j]**2)
        
    for i in range(0,nx):
        k_base[i] =  k[0,i]
    
    spettro_f = []
    for i in range(tempi+1):
        spettro_f.append( np.zeros(nx+1,dtype=float) )
    
    for t in range(tempi+1):
        for i in range(0,nx):
            for j in range(0,ny):      
                l = np.argmin(np.abs(k_base-k[j,i]))
                dist = k_base[l]-k[j,i]
                if dist>0 :
                    spettro_f[t][l] += f_in[t][0,j,i]*(1-dist)
                    spettro_f[t][l-1] += f_in[t][0,j,i]*(dist)
                elif dist<0 :
                    spettro_f[t][l] += f_in[t][0,j,i]*(1-np.abs(dist))
                    spettro_f[t][l+1] += f_in[t][0,j,i]*(np.abs(dist))
                else:# (k_base - k[j,i])=0 :
                    spettro_f[t][l] += f_in[t][0,j,i]
                    
    for t in range(tempi+1):
        spettro_f_finale.append(np.delete(spettro_f[t],[nx-1,nx]))      #l'ultimo bin lo tolgo perchè...
    
    end =time.time()
    print()
    print('END - radial_spectrum_2D')
    timing(start,end)
    
    return spettro_f_finale , k_base

'''###############################################################################'''
@jit(cache=True)
def radial_spectrum_2D_compressible_solenoidal(fx,fy,fz,tempi):
    start=time.time()
    
    fxt= [] ; fyt = [] ; fzt =[]
    for t in range(tempi+1):
        fxt.append(np.fft.rfftn(fx[t]))
        fyt.append(np.fft.rfftn(fy[t]))
        fzt.append(np.fft.rfftn(fz[t]))
    
    ny = len( fxt[0][0,:,0] ) ; nx = len( fyt[0][0,0,:] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int)
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
    for i in range(nx,ny):
        ky[i]=-ny+i   
    
    spettro_fc_finale = [] ; spettro_fs_finale = [] ; k_base=np.zeros(nx,dtype=float)
    
    k = np.zeros((ny,nx),dtype=float)
    for j in range(ny):
        for i in range(nx):
            k[j,i] = math.sqrt(kx[i]**2+ky[j]**2)
    
    compr=[] ; solen=[]
    k[0,0]=1.0 #; kx[0]=1.0 ; ky[0]=1.0 #cosi non scazza il punto 0,0
    for t in range(tempi+1):
        compr.append( np.abs( ( ( kx*fxt[t] + np.transpose((ky*np.transpose(fyt[t],(0,2,1))),(0,2,1)) )/k )**2) )
        solen.append( np.abs( ( (np.transpose((ky*np.transpose(fzt[t],(0,2,1))),(0,2,1)))**2 + (kx*fzt[t])**2 \
                               + ( kx*fyt[t] - np.transpose((ky*np.transpose(fxt[t],(0,2,1))),(0,2,1)) )**2 )/k**2 ) )
    k[0,0]=0.0 #; kx[0]=0.0 ; ky[0]=0.0  #cosi rimetto a posto il giusto k
    
    for i in range(0,nx):
        k_base[i] =  k[0,i]
    
    spettro_fc = [] ; spettro_fs = [] 
    for i in range(tempi+1):
        spettro_fc.append( np.zeros(nx+1,dtype=float) )
        spettro_fs.append( np.zeros(nx+1,dtype=float) )
    
    for t in range(tempi+1):
        for i in range(0,nx):
            for j in range(0,ny):      
                l = np.argmin(np.abs(k_base-k[j,i]))
                dist = k_base[l]-k[j,i]
                if dist>0 :
                    spettro_fc[t][l] += compr[t][0,j,i]*(1-dist)
                    spettro_fc[t][l-1] += compr[t][0,j,i]*(dist)
                    spettro_fs[t][l] += solen[t][0,j,i]*(1-dist)
                    spettro_fs[t][l-1] += solen[t][0,j,i]*(dist)
                elif dist<0 :
                    spettro_fc[t][l] += compr[t][0,j,i]*(1-np.abs(dist))
                    spettro_fc[t][l+1] += compr[t][0,j,i]*(np.abs(dist))
                    spettro_fs[t][l] += solen[t][0,j,i]*(1-np.abs(dist))
                    spettro_fs[t][l+1] += solen[t][0,j,i]*(np.abs(dist))
                else:# (k_base - k[j,i])=0 :
                    spettro_fc[t][l] += compr[t][0,j,i]
                    spettro_fs[t][l] += solen[t][0,j,i]

    for t in range(tempi+1):
        spettro_fc_finale.append(np.delete(spettro_fc[t],[nx-1,nx]))      #l'ultimo bin lo tolgo perchè...
        spettro_fs_finale.append(np.delete(spettro_fs[t],[nx-1,nx]))
    
    end =time.time()
    print()
    print('END - radial_spectrum_2D_compressible_solenoidal')
    timing(start,end)
    
    return spettro_fc_finale , spettro_fs_finale, k_base

'''###############################################################################'''
@jit(cache=True)
def linee_campo_B_2D(bx,by,tempi):
    start=time.time()
    
    bxf=[];byf=[] 
    for t in range(tempi+1):
        bxf.append(np.fft.rfftn(bx[t]))
        byf.append(np.fft.rfftn(by[t]))
        
    ny = len( bxf[0][0,:,0] ) ; nx = len( bxf[0][0,0,:] ) 
    inv_kx = np.zeros(nx, dtype=float) ; inv_ky = np.zeros(ny, dtype=float) #il modo 0 invertendolo mi da zero... cosi non divido per zero.
    
    for i in range(1,nx):
        inv_kx[i]=1/i ; inv_ky[i]=1/i 
    for i in range(nx,ny):
        inv_ky[i]=1/(-ny+i)  
        
    I=0+1j
    f_Az_x=[] ; f_Az_y=[] ; Az_x=[] ; Az_y=[] ; Az=[]
    for t in range(tempi+1):
        f_Az_x.append(-I*(np.transpose( (np.transpose(bxf[t],(0,2,1)))*inv_ky,(0,2,1) )))
        f_Az_y.append(I*byf[t]*inv_kx)
        Az_x.append(np.fft.irfftn(f_Az_x[t]))
        Az_y.append(np.fft.irfftn(f_Az_y[t]))
        Az.append( (Az_x[t] + Az_y[t])/2 )
    
    end =time.time()
    print()
    print('END - linee_campo_B_2D')
    timing(start,end)
    
    return Az_x , Az_y , Az

'''###############################################################################'''
@jit(cache=True)
def gyrotropic_spectrum_3D(f_in,tempi):
    start=time.time()
    
    ny = len( f_in[0][0,:,0] ) ; nx = len( f_in[0][0,0,:] ) ; nz = len( f_in[0][:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i; kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    for i in range(nx,nz):
        kz[i]=-ny+i
    
    spettro_f_finale = []
        
    k_perp = np.zeros((nz,ny),dtype=float)
    for j in range(nz):
        for i in range(ny):
            k_perp[j,i] = math.sqrt(kz[j]**2+ky[i]**2)
            
    k_base = np.zeros(nx,dtype=float)
    for i in range(0,nx):
        k_base[i] =  k_perp[0,i]
    
    spettro_f = []
    for i in range(tempi+1):
        spettro_f.append( np.zeros((nx,nx+1),dtype=float) )
    
    for t in range(tempi+1):
        for j in range(0,nz):
            for i in range(0,ny): 
                for w in range(0,nx):
                    l = np.argmin(np.abs((k_base-k_perp[j,i])))
                    dist = k_base[l]-k_perp[j,i]
                    if dist>0 :
                        spettro_f[t][w,l] += f_in[t][j,i,w]*(1-dist)
                        spettro_f[t][w,l-1] += f_in[t][j,i,w]*(dist)
                    elif dist<0 :
                        spettro_f[t][w,l] += f_in[t][j,i,w]*(1-np.abs(dist))
                        spettro_f[t][w,l+1] += f_in[t][j,i,w]*(np.abs(dist))
                    else:# (k_base - k[j,i])=0 :
                        spettro_f[t][w,l] += f_in[t][j,i,w]
    
    for t in range(tempi+1):
        spettro_f_finale.append( np.delete(spettro_f[t],nx,axis=1))     

    end =time.time()
    print()
    print('END - gyrotropic_spectrum_3D')
    timing(start,end)
    
    return spettro_f_finale , kx , k_base

'''###############################################################################'''   
@jit(cache=True)
def spettri_ridotti(f_in,spettro_in,tempi):
    start=time.time()
    
    ny = len( f_in[0][0,:,0] ) ; nx = len( f_in[0][0,0,:] ) ; nz = len( f_in[0][:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i; kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    for i in range(nx,nz):
        kz[i]=-ny+i
    
    spettro_f_finale_iso=[]; spettro_f_finale_perp=[]; spettro_f_finale_parall=[]
            
    k_mod = np.zeros((nx,nx),dtype=float) #raggio sfera
    for i in range(nx):
        for j in range(nx):
            k_mod[j,i] = math.sqrt(kx[j]**2+kx[i]**2)
      
    k_base = np.zeros(nx,dtype=float)          
    for i in range(0,nx):
        k_base[i] =  k_mod[0,i]
    
    spettro_f_iso = []
    for i in range(tempi+1):
        spettro_f_iso.append( np.zeros(nx+1,dtype=float) )
        
    for t in range(tempi+1):
        for i in range(0,nx):
            for j in range(0,nx):      
                l = np.argmin(np.abs(k_base-k_mod[j,i]))
                dist = k_base[l]-k_mod[j,i]
                if dist>0 :
                    spettro_f_iso[t][l] += spettro_in[t][j,i]*(1-dist)
                    spettro_f_iso[t][l-1] += spettro_in[t][j,i]*(dist)
                elif dist<0 :
                    spettro_f_iso[t][l] += spettro_in[t][j,i]*(1-np.abs(dist))
                    spettro_f_iso[t][l+1] += spettro_in[t][j,i]*(np.abs(dist))
                else:# (k_base - k[j,i])=0 :
                    spettro_f_iso[t][l] += spettro_in[t][j,i]
                                     
    for t in range(tempi+1):
        spettro_f_finale_iso.append(np.delete(spettro_f_iso[t],[nx-1,nx])) 

    spettro_f_perp_rid = [] ; spettro_f_parall_rid = []
    for i in range(tempi+1):
        spettro_f_perp_rid.append( np.zeros(nx,dtype=float) )
        spettro_f_parall_rid.append( np.zeros(nx,dtype=float) )
        
    for t in range(tempi+1):
        for j in range(0,nx):
            for i in range(0,nx):      
                spettro_f_perp_rid[t][i] += spettro_in[t][j,i] 
                spettro_f_parall_rid[t][j] += spettro_in[t][j,i] 
                
    for t in range(tempi+1):
        spettro_f_finale_perp.append(np.delete(spettro_f_perp_rid[t],[nx-1]))     
        spettro_f_finale_parall.append(np.delete(spettro_f_parall_rid[t],[nx-1]))
        
    end =time.time()
    print()
    print('END - spettri_ridotti')
    timing(start,end)

    return spettro_f_finale_iso, spettro_f_finale_perp, spettro_f_finale_parall

'''###############################################################################'''
#@jit(cache=True)
def Curl(fx,fy,fz,tempi):
    start=time.time()
    
    ny = len( fx[0][0,:,0] ) ; nx = len( fx[0][0,0,:] ) ; nz = len( fx[0][:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
        if nz>nx: kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    if nz>nx:
        for i in range(nx,nz):
            kz[i]=-ny+i
    
    Jt_x= [] ; Jt_y = [] ; Jt_z =[]
    for t in range(tempi+1):
        Jt_x.append(np.zeros((nz,ny,nx),dtype=complex))
        Jt_y.append(np.zeros((nz,ny,nx),dtype=complex))
        Jt_z.append(np.zeros((nz,ny,nx),dtype=complex))
    
    J_x= [] ; J_y = [] ; J_z =[]
    bxt= [] ; byt = [] ; bzt =[]
    
    I=0+1j
    
    for t in range(tempi+1):
        bxt.append(np.fft.rfftn(fx[t]))
        byt.append(np.fft.rfftn(fy[t]))
        bzt.append(np.fft.rfftn(fz[t]))
        for w in range(nz):
            for j in range(ny):
                for i in range(nx):
                    Jt_x[t][w,j,i] = ( ky[j]*bzt[t][w,j,i] - kz[w]*byt[t][w,j,i] )*I
                    Jt_y[t][w,j,i] = ( kz[w]*bxt[t][w,j,i] - kx[i]*bzt[t][w,j,i] )*I
                    Jt_z[t][w,j,i] = ( kx[i]*byt[t][w,j,i] - ky[j]*bxt[t][w,j,i] )*I
    for t in range(tempi+1):    
        J_x.append(np.fft.irfftn(Jt_x[t]))
        J_y.append(np.fft.irfftn(Jt_y[t]))
        J_z.append(np.fft.irfftn(Jt_z[t]))
    
    end =time.time()
    print()
    print('END - Curl')
    timing(start,end)
    
    return J_x , J_y , J_z

'''###############################################################################'''
#@jit(cache=True)
def Curl_migliore(fx,fy,fz,tempi):
    start=time.time()
    
    ny = len( fx[0][0,:,0] ) ; nx = len( fx[0][0,0,:] ) ; nz = len( fx[0][:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
        if nz>nx: kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    if nz>nx:
        for i in range(nx,nz):
            kz[i]=-ny+i
    
    Jt_x1= [] ; Jt_y1 = [] ; Jt_z1 =[] ; Jt_x2= [] ; Jt_y2 = [] ; Jt_z2 =[]
    for t in range(tempi+1):
        Jt_x1.append(np.zeros((nz,ny,nx),dtype=complex))
        Jt_y1.append(np.zeros((nz,ny,nx),dtype=complex))
        Jt_z1.append(np.zeros((nz,ny,nx),dtype=complex))
        Jt_x2.append(np.zeros((nz,ny,nx),dtype=complex))
        Jt_y2.append(np.zeros((nz,ny,nx),dtype=complex))
        Jt_z2.append(np.zeros((nz,ny,nx),dtype=complex))
        
    J_x= [] ; J_y = [] ; J_z =[] ; bxt= [] ; byt = [] ; bzt =[]
    I=0+1j
    
    for t in range(tempi+1):
        bxt.append(np.fft.rfftn(fx[t]))
        byt.append(np.fft.rfftn(fy[t]))
        bzt.append(np.fft.rfftn(fz[t]))
        for w in range(nz):
            Jt_x2[t][w,:,:] = ( kz[w]*byt[t][w,:,:] ) 
            Jt_y1[t][w,:,:] = ( kz[w]*bxt[t][w,:,:] )
        for j in range(ny):
            Jt_x1[t][:,j,:] = ( ky[j]*bzt[t][:,j,:] )    
            Jt_z2[t][:,j,:] = ( ky[j]*bxt[t][:,j,:] )
        for i in range(nx):  
            Jt_y2[t][:,:,i] = ( kx[i]*bzt[t][:,:,i] )
            Jt_z1[t][:,:,i] = ( kx[i]*byt[t][:,:,i] )   
        J_x.append(np.fft.irfftn( (Jt_x1[t] - Jt_x2[t])*I ))
        J_y.append(np.fft.irfftn( (Jt_y1[t] - Jt_y2[t])*I ))
        J_z.append(np.fft.irfftn( (Jt_z1[t] - Jt_z2[t])*I ))  
    
    end =time.time()
    print()
    print('END - Curl_migliore')
    timing(start,end)
    
    return J_x , J_y , J_z

'''###############################################################################'''
@jit(cache=True)
def Curl_migliorissimo(fx,fy,fz,tempi): #richiede troppa ram quindi sl mio pc non migliora quanto dovrebbe
    start=time.time()
    
    bxt= [] ; byt = [] ; bzt =[]
    for t in range(tempi+1):
        bxt.append(np.fft.rfftn(fx[t]))
        byt.append(np.fft.rfftn(fy[t]))
        bzt.append(np.fft.rfftn(fz[t]))
    
    ny = len( bxt[0][0,:,0] ) ; nx = len( bxt[0][0,0,:] ) ; nz = len( bxt[0][:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
        if nz>nx: kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    if nz>nx:
        for i in range(nx,nz):
            kz[i]=-ny+i
    
#    Jt_x1= [] ; Jt_y1 = [] ; Jt_z1 =[] ; Jt_x2= [] ; Jt_y2 = [] ; Jt_z2 =[]
#    for t in range(tempi+1):
#        Jt_x1.append(np.zeros((nz,ny,nx),dtype=complex))
#        Jt_y1.append(np.zeros((nz,ny,nx),dtype=complex))
#        Jt_z1.append(np.zeros((nz,ny,nx),dtype=complex))
#        Jt_x2.append(np.zeros((nz,ny,nx),dtype=complex))
#        Jt_y2.append(np.zeros((nz,ny,nx),dtype=complex))
#        Jt_z2.append(np.zeros((nz,ny,nx),dtype=complex))
#        
#    J_x= [] ; J_y = [] ; J_z =[] 
#    I=0+1j
#    
#    for t in range(tempi+1):
#        Jt_y2[t] = ( kx*bzt[t] )
#        Jt_z1[t] = ( kx*byt[t] )     
#        Jt_x2[t] = np.transpose( (kz*np.transpose(byt[t],(2,1,0))) , (2,1,0) ) 
#        Jt_y1[t] = np.transpose( (kz*np.transpose(bxt[t],(2,1,0))) , (2,1,0) )      
#        Jt_x1[t] = np.transpose( (ky*np.transpose(bzt[t],(0,2,1))) , (0,2,1) )
#        Jt_z2[t] = np.transpose( (ky*np.transpose(bxt[t],(0,2,1))) , (0,2,1) )        
#        J_x.append(np.fft.irfftn( (Jt_x1[t] - Jt_x2[t])*I ))
#        J_y.append(np.fft.irfftn( (Jt_y1[t] - Jt_y2[t])*I ))
#        J_z.append(np.fft.irfftn( (Jt_z1[t] - Jt_z2[t])*I ))  

    J_x= [] ; J_y = [] ; J_z =[] 
    I=0-1j  ;  ''' COL MENO É GIUSTO '''
    
    for t in range(tempi+1):   
        J_x.append(np.fft.irfftn( (np.transpose( (ky*np.transpose(bzt[t],(0,2,1))) , (0,2,1) ) - np.transpose( (kz*np.transpose(byt[t],(2,1,0))) , (2,1,0) ))*I ))
        J_y.append(np.fft.irfftn( (np.transpose( (kz*np.transpose(bxt[t],(2,1,0))) , (2,1,0) ) - kx*bzt[t])*I ))
        J_z.append(np.fft.irfftn( (kx*byt[t] - np.transpose( (ky*np.transpose(bxt[t],(0,2,1))) , (0,2,1) ))*I ))  
    
    end =time.time()
    print()
    print('END - Curl_migliorissimo')
    timing(start,end) 
    
    return J_x , J_y , J_z

'''###############################################################################'''
@jit(cache=True)
def Divergence(fx,fy,fz,tempi): #richiede troppa ram quindi sl mio pc non migliora quanto dovrebbe
    start=time.time()
    
    uxt= [] ; uyt = [] ; uzt =[]
    for t in range(tempi+1):
        uxt.append(np.fft.rfftn(fx[t]))
        uyt.append(np.fft.rfftn(fy[t]))
        uzt.append(np.fft.rfftn(fz[t]))
    
    ny = len( uxt[0][0,:,0] ) ; nx = len( uxt[0][0,0,:] ) ; nz = len( uxt[0][:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
        if nz>nx: kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    if nz>nx:
        for i in range(nx,nz):
            kz[i]=-ny+i

    div=[]
    I=0-1j     ; ''' COL MENO É GIUSTO '''
    
    for t in range(tempi+1):   
        div.append(np.fft.irfftn( (kx*uxt[t] + np.transpose((ky*np.transpose(uyt[t],(0,2,1))),(0,2,1)) + np.transpose((kz*np.transpose(uzt[t],(2,1,0))),(2,1,0)))*I ))
    
    end =time.time()
    print()
    print('END - Divergence')
    timing(start,end) 
    
    return div

'''###############################################################################'''
@jit(cache=True)
def Potenziale_Vettore(Jx,Jy,Jz,tempi): #richiede troppa ram quindi sl mio pc non migliora quanto dovrebbe
    start=time.time()
    
    jxt= [] ; jyt = [] ; jzt =[]
    for t in range(tempi+1):
        jxt.append(np.fft.rfftn(Jx[t]))
        jyt.append(np.fft.rfftn(Jy[t]))
        jzt.append(np.fft.rfftn(Jz[t]))
    
    ny = len( jxt[0][0,:,0] ) ; nx = len( jxt[0][0,0,:] ) ; nz = len( jxt[0][:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
        if nz>nx: kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    if nz>nx:
        for i in range(nx,nz):
            kz[i]=-ny+i
    
    invK2=np.ones((nz,ny,nx),dtype=float)
    for w in range(nz):
        for j in range(ny):
            for i in range(nx):
                if (kz[w]==0 and ky[j]==0 and kx[i]==0): 
                    print('salto kx=0, ky=0, kz=0')
                else:
                    invK2[w,j,i]=1/(kx[i]**2+ky[j]**2+kz[w]**2)
                
    A_x= [] ; A_y = [] ; A_z =[] 
    
    for t in range(tempi+1):   
        A_x.append(np.fft.irfftn( jxt[t]*invK2 ))
        A_y.append(np.fft.irfftn( jyt[t]*invK2 ))
        A_z.append(np.fft.irfftn( jzt[t]*invK2 )) 
    
    end =time.time()
    print()
    print('END - Potenziale_Vettore')
    timing(start,end) 
    
    return A_x , A_y , A_z

'''###############################################################################'''
@jit(cache=True)
def Potenziale_Vettore_daB(bx,by,bz,tempi): #richiede troppa ram quindi sl mio pc non migliora quanto dovrebbe
    start=time.time()
    
    bxt= [] ; byt = [] ; bzt =[]
    for t in range(tempi+1):
        bxt.append(np.fft.rfftn(bx[t]))
        byt.append(np.fft.rfftn(by[t]))
        bzt.append(np.fft.rfftn(bz[t]))
    
    ny = len( bxt[0][0,:,0] ) ; nx = len( bxt[0][0,0,:] ) ; nz = len( bxt[0][:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
        if nz>nx: kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    if nz>nx:
        for i in range(nx,nz):
            kz[i]=-ny+i
    
    invK2=np.ones((nz,ny,nx),dtype=float)
    for w in range(nz):
        for j in range(ny):
            for i in range(nx):
                if (w==0 and j==0 and i==0): 
                    print('salto kx=0, ky=0, kz=0')
                else:
                    invK2[w,j,i]=1/(kx[i]**2+ky[j]**2+kz[w]**2)
                
    A_x= [] ; A_y = [] ; A_z =[] 

    I=0+1j
    
    for t in range(tempi+1):   
        A_x.append(np.fft.irfftn( (np.transpose( (ky*np.transpose(bzt[t],(0,2,1))) , (0,2,1) ) - np.transpose( (kz*np.transpose(byt[t],(2,1,0))) , (2,1,0) ))*I*invK2 ))
        A_y.append(np.fft.irfftn( (np.transpose( (kz*np.transpose(bxt[t],(2,1,0))) , (2,1,0) ) - kx*bzt[t])*I*invK2 ))
        A_z.append(np.fft.irfftn( (kx*byt[t] - np.transpose( (ky*np.transpose(bxt[t],(0,2,1))) , (0,2,1) ))*I*invK2 )) 
    
    end =time.time()
    print()
    print('END - Potenziale_Vettore_daB')
    timing(start,end) 
    
    return A_x , A_y , A_z

'''##############################################################################'''
@jit(cache=True)
def Varianza(f_in,tempi):
    start=time.time()
    
    ny = len( f_in[0][0,:,0] ) ; nx = len( f_in[0][0,0,:] ) ; nz = len( f_in[0][:,0,0] )
    f_somma=np.zeros(tempi+1,dtype=float) ; f_quadratico_somma=np.zeros(tempi+1,dtype=float)
    
    for t in range(tempi+1):
        for w in range(nz):
            for j in range(ny):
                for i in range(nx):
                    f_somma[t] += f_in[t][w,j,i]
                    f_quadratico_somma[t] += (f_in[t][w,j,i])**2
    
    f_medio = f_somma/(nx*ny*nz)
    f_medio_quadro = (f_medio)**2
    f_quadratico_medio = f_quadratico_somma/(nx*ny*nz)
    Varianza_f = f_quadratico_medio - f_medio_quadro
                
    end =time.time()
    print()
    print('END - Varianza')
    timing(start,end)
    
    return Varianza_f

'''##############################################################################'''
@jit(cache=True) #questo funziona, ne son sicuro
def Pdf_2D(f_in,tempi):    
    start=time.time()

    ny = len( f_in[0][0,:,0] ) ; nx = len( f_in[0][0,0,:] )

    pdf_x = [] ; pdf_y = []
    for t in range(tempi+1):
        pdf_x.append([]) ; pdf_y.append([])
        for l in range(1,int(nx/2)+1):
            pdf_x[t].append([])
        for l in range(1,int(ny/2)+1):
            pdf_y[t].append([])
    
    for t in range(tempi+1):
        for l in range(1,int(nx/2)+1):
            for i in range(nx):
                if i+l < nx :
                    pdf_x[t][l-1].append( f_in[t][0,:,i+l] - f_in[t][0,:,i] )  
                elif i+l >= nx :
                    pdf_x[t][l-1].append( f_in[t][0,:,i+l-nx] - f_in[t][0,:,i] )  
            pdf_x[t][l-1] = np.asarray(pdf_x[t][l-1]).reshape(-1)
        for l in range(1,int(ny/2)+1):
            for i in range(ny):
                if i+l < ny :
                    pdf_y[t][l-1].append( f_in[t][0,i+1,:] - f_in[t][0,i,:] )  
                elif i+l >= nx :
                    pdf_y[t][l-1].append( f_in[t][0,i+l-ny,:] - f_in[t][0,i,:] )  
            pdf_y[t][l-1] = np.asarray(pdf_y[t][l-1]).reshape(-1)
        
    end=time.time()
    print()
    print('END - Pdf_2D')
    timing(start,end)
    
    return pdf_x , pdf_y

'''##############################################################################'''
@jit(cache=True) #questo funziona, ne son sicuro
def Pdf(f_in,tempi,flag):    
    start=time.time()

    ny = len( f_in[0][0,:,0] ) ; nx = len( f_in[0][0,0,:] ) ; nz = len( f_in[0][0,0,:] ) ; pdf=[]

    if flag=='x':
        n=nx
        for t in range(tempi+1):
            pdf.append([])
            for l in range(1,int(n/2)+1):
                pdf[t].append([])    
                
        for t in range(tempi+1):
            for l in range(1,int(n/2)+1):
                for i in range(n):
                    if i+l < n :
                        pdf[t][l-1].append( f_in[t][:,:,i+l] - f_in[t][:,:,i] )  
                    elif i+l >= n :
                        pdf[t][l-1].append( f_in[t][:,:,i+l-n] - f_in[t][:,:,i] )  
                pdf[t][l-1] = np.asarray(pdf[t][l-1]).reshape(-1)
        
    elif flag=='y': 
        n=ny
        for t in range(tempi+1):
            pdf.append([])
            for l in range(1,int(n/2)+1):
                pdf[t].append([])
        
        for t in range(tempi+1):
            for l in range(1,int(n/2)+1):
                for i in range(n):
                    if i+l < n :
                        pdf[t][l-1].append( f_in[t][:,i+l,:] - f_in[t][:,i,:] )  
                    elif i+l >= n :
                        pdf[t][l-1].append( f_in[t][:,i+l-n,:] - f_in[t][:,i,:] )  
                pdf[t][l-1] = np.asarray(pdf[t][l-1]).reshape(-1)
    
    elif flag=='z': 
        n=nz
        for t in range(tempi+1):
            pdf.append([])
            for l in range(1,int(n/2)+1):
                pdf[t].append([])
        
        for t in range(tempi+1):
            for l in range(1,int(n/2)+1):
                for i in range(n):
                    if i+l < n :
                        pdf[t][l-1].append( f_in[t][i+l,:,:] - f_in[t][i,:,:] )  
                    elif i+l >= n :
                        pdf[t][l-1].append( f_in[t][i+l-n,:,:] - f_in[t][i,:,:] )  
                pdf[t][l-1] = np.asarray(pdf[t][l-1]).reshape(-1)
                
    end=time.time()
    print()
    print('END - Pdf')
    timing(start,end)
    
    return pdf

'''##############################################################################'''
@jit(cache=True)
def mom_2D_miei(pdf_in,tempi,bins): #sono più lenti e meno precisi della funzione di python :(
    start=time.time()
    
    hist=[] ; base_hist=[]
    for t in range(tempi+1):
        hist.append(np.histogram(pdf_in[t],bins=bins)[0])
        base_hist.append(np.histogram(pdf_in[t],bins=bins)[1])

    mom=np.zeros((tempi+1,5),dtype=float)
    for t in range(tempi+1):
        for j in range(5):
            for i in range(len(hist[t])):
                mom[t,j] += hist[t][i]*(base_hist[t][i+1]-base_hist[t][i])*( ((base_hist[t][i+1]+base_hist[t][i])*0.5)**j ) 
        mom[t,:] = mom[t,:]/mom[t,0]
    
    
    end =time.time()
    print()
    print('END - mom_2D_miei')
    timing(start,end)
    
    return mom

'''##############################################################################'''
@jit(cache=True)
def momenti(pdf_in,tempi):
    start=time.time()
    nx = 2*len( pdf_in[0][:] )
            
    mom = np.zeros((tempi+1,int(nx/2),5),dtype=float)
    for t in range(tempi+1):
        for l in range(int(nx/2)):
            for j in range(5):
                mom[t,l,j] = stats.moment(pdf_in[t][l],moment=j)
            mom[t,l,:] = mom[t,l,:]/mom[t,l,0]
    
    K = mom[:,:,4]/((mom[:,:,2])**2) - 3 
    S = mom[:,:,3]/((mom[:,:,2])**1.5)
    
    end =time.time()
    print()
    print('END - momenti')
    timing(start,end)
    
    return mom , K , S 

'''##############################################################################'''
@jit(cache=True)
def momenti_distribuzione(arr,tempi):
    start=time.time()
            
    mom = np.zeros((tempi+1,5),dtype=float)
    for t in range(tempi+1):
        for j in range(5):
            mom[t,j] = stats.moment(arr[t],moment=j)
        mom[t,:] = mom[t,:]/mom[t,0]
    
    K = mom[:,4]/((mom[:,2])**2) - 3 
    S = mom[:,3]/((mom[:,2])**1.5)
    
    end =time.time()
    print()
    print('END - momenti_distribuzione')
    timing(start,end)
    
    return mom , K , S 
