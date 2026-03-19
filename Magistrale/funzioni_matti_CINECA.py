# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 10:57:46 2019

@author: Mattia

In questa versione le funzioni vogliono array e non liste di array.
Mi serve per output troppo grossi (es: 512^3)
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
    
    nz=len(f_in[:,0,0]) ; ny=len(f_in[0,:,0]) ; nx=len(f_in[0,0,:])
    
    tot=np.argmax(f_in)
    maxz = (int(tot/(ny*nz)))
    maxy = (int((tot-maxz*nx*ny)/ny))
    maxx = (tot-maxz*nx*ny-maxy*ny)
    
    tot=np.argmin(f_in)
    minz = (int(tot/(ny*nz)))
    miny = (int((tot-minz*nx*ny)/ny))
    minx = (tot-minz*nx*ny-miny*ny)

    end=time.time()
    print()
    print('END - MaxMinMatrix')
    timing(start,end)
        
    return maxz , maxy , maxx , minz , miny , minx

'''###############################################################################'''
@jit(cache=True)
def Media(f_in):
    
    start=time.time()
    nz=len(f_in[:,0,0]) ; ny=len(f_in[0,:,0]) ; nx=len(f_in[0,0,:])
    punti = nx*ny*nz
    
    ava_f=0.0
    for w in range(nz):
        for j in range(ny):
            for i in range(nx):
                ava_f += f_in[w,j,i]
    ava_f = ava_f/punti 
    
    end=time.time()
    print()
    print('END - Media')
    timing(start,end)
        
    return ava_f
    
'''###############################################################################'''
@jit(cache=True)
def radial_spectrum_primi_vicini_2D(f_in):
    start=time.time()
    
    
    ny = len( f_in[0,:,0] ) ; nx = len( f_in[0,0,:] ) 
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) 
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
    for i in range(nx,ny):
        ky[i]=-ny+i   
    
    k = np.zeros((ny,nx),dtype=float)
    for j in range(ny):
        for i in range(nx):
            k[j,i] = math.sqrt(kx[i]**2+ky[j]**2)
        
    k_base = k[0,:]
    
    spettro_f = np.zeros(nx,dtype=float) 
    
    for i in range(0,nx):
        for j in range(0,ny):      
            l = np.argmin(np.abs((k_base-k[j,i])))
            spettro_f[l] += f_in[0,j,i]
                  

    spettro_f_finale = np.delete(spettro_f,[nx-1])      #l'ultimo bin lo tolgo perchè...
     
    end =time.time()
    print()
    print('END - radial_spectrum_primi_vicini_2D')
    timing(start,end)
        
    return spettro_f_finale , k_base    

'''###############################################################################'''
@jit(cache=True)
def radial_spectrum_2D(f_in):
    start=time.time()
    
    ny = len( f_in[0,:,0] ) ; nx = len( f_in[0,0,:] ) 
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) 
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
    for i in range(nx,ny):
        ky[i]=-ny+i   
    
    k = np.zeros((ny,nx),dtype=float)
    for j in range(ny):
        for i in range(nx):
            k[j,i] = math.sqrt(kx[i]**2+ky[j]**2)
        
    k_base =  k[0,:]
    
    spettro_f = np.zeros(nx+1,dtype=float) 
    

    for i in range(0,nx):
        for j in range(0,ny):      
            l = np.argmin(np.abs(k_base-k[j,i]))
            dist = k_base[l]-k[j,i]
            if dist>0 :
                spettro_f[l] += f_in[0,j,i]*(1-dist)
                spettro_f[l-1] += f_in[0,j,i]*(dist)
            elif dist<0 :
                spettro_f[l] += f_in[0,j,i]*(1-np.abs(dist))
                spettro_f[l+1] += f_in[0,j,i]*(np.abs(dist))
            else:# (k_base - k[j,i])=0 :
                spettro_f[l] += f_in[0,j,i]
                    
    spettro_f_finale = np.delete(spettro_f,[nx-1,nx])      #l'ultimo bin lo tolgo perchè...
    
    end =time.time()
    print()
    print('END - radial_spectrum_2D')
    timing(start,end)
    
    return spettro_f_finale , k_base

'''###############################################################################'''
@jit(cache=True)
def radial_spectrum_2D_compressible_solenoidal(fx,fy,fz):
    start=time.time()
    
    fxt = np.fft.rfftn(fx)
    fyt = np.fft.rfftn(fy)
    fzt = np.fft.rfftn(fz)
    
    ny = len( fxt[0,:,0] ) ; nx = len( fyt[0,0,:] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int)
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
    for i in range(nx,ny):
        ky[i]=-ny+i   
    
    k_base=np.zeros(nx,dtype=float)
    
    k = np.zeros((ny,nx),dtype=float)
    for j in range(ny):
        for i in range(nx):
            k[j,i] = math.sqrt(kx[i]**2+ky[j]**2)
    
    k[0,0]=1.0 #; kx[0]=1.0 ; ky[0]=1.0 #cosi non scazza il punto 0,0
    compr = ( np.abs( ( ( kx*fxt + np.transpose((ky*np.transpose(fyt,(0,2,1))),(0,2,1)) )/k )**2) )
    solen = ( np.abs( ( (np.transpose((ky*np.transpose(fzt,(0,2,1))),(0,2,1)))**2 + (kx*fzt)**2 \
                           + ( kx*fyt - np.transpose((ky*np.transpose(fxt,(0,2,1))),(0,2,1)) )**2 )/k**2 ) )
    k[0,0]=0.0 #; kx[0]=0.0 ; ky[0]=0.0  #cosi rimetto a posto il giusto k
    
    for i in range(0,nx):
        k_base[i] =  k[0,i]
    
    spettro_fc = ( np.zeros(nx+1,dtype=float) )
    spettro_fs = ( np.zeros(nx+1,dtype=float) )
    
    for i in range(0,nx):
        for j in range(0,ny):      
            l = np.argmin(np.abs(k_base-k[j,i]))
            dist = k_base[l]-k[j,i]
            if dist>0 :
                spettro_fc[l] += compr[0,j,i]*(1-dist)
                spettro_fc[l-1] += compr[0,j,i]*(dist)
                spettro_fs[l] += solen[0,j,i]*(1-dist)
                spettro_fs[l-1] += solen[0,j,i]*(dist)
            elif dist<0 :
                spettro_fc[l] += compr[0,j,i]*(1-np.abs(dist))
                spettro_fc[l+1] += compr[0,j,i]*(np.abs(dist))
                spettro_fs[l] += solen[0,j,i]*(1-np.abs(dist))
                spettro_fs[l+1] += solen[0,j,i]*(np.abs(dist))
            else:# (k_base - k[j,i])=0 :
                spettro_fc[l] += compr[0,j,i]
                spettro_fs[l] += solen[0,j,i]

    spettro_fc_finale = (np.delete(spettro_fc,[nx-1,nx]))      #l'ultimo bin lo tolgo perchè...
    spettro_fs_finale = (np.delete(spettro_fs,[nx-1,nx]))
    
    end =time.time()
    print()
    print('END - radial_spectrum_2D_compressible_solenoidal')
    timing(start,end)
    
    return spettro_fc_finale , spettro_fs_finale, k_base

'''###############################################################################'''
@jit(cache=True)
def linee_campo_B_2D(bx,by):
    start=time.time()
    
    bxf = np.fft.rfftn(bx)
    byf = np.fft.rfftn(by)
        
    ny = len( bxf[0,:,0] ) ; nx = len( bxf[0,0,:] ) 
    inv_kx = np.zeros(nx, dtype=float) ; inv_ky = np.zeros(ny, dtype=float) #il modo 0 invertendolo mi da zero... cosi non divido per zero.
    
    for i in range(1,nx):
        inv_kx[i]=1/i ; inv_ky[i]=1/i 
    for i in range(nx,ny):
        inv_ky[i]=1/(-ny+i)  
        
    I=0+1j

    f_Az_x = -I*(np.transpose( (np.transpose(bxf,(0,2,1)))*inv_ky,(0,2,1) ))
    f_Az_y =  I*byf*inv_kx
    Az_x   = np.fft.irfftn(f_Az_x)
    Az_y   = np.fft.irfftn(f_Az_y)
    Az     =  (Az_x + Az_y)/2 
    
    end =time.time()
    print()
    print('END - linee_campo_B_2D')
    timing(start,end)
    
    return Az_x , Az_y , Az

'''###############################################################################'''
@jit(cache=True)
def gyrotropic_spectrum_3D(f_in):
    start=time.time()
    
    ny = len( f_in[0,:,0] ) ; nx = len( f_in[0,0,:] ) ; nz = len( f_in[:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i; kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    for i in range(nx,nz):
        kz[i]=-ny+i
        
    k_perp = np.zeros((nz,ny),dtype=float)
    for j in range(nz):
        for i in range(ny):
            k_perp[j,i] = math.sqrt(kz[j]**2+ky[i]**2)
            
    k_base =  k_perp[0,:]
    
    spettro_f =  np.zeros((nx,nx+1),dtype=float) 

    for j in range(0,nz):
        for i in range(0,ny): 
            for w in range(0,nx):
                l = np.argmin(np.abs((k_base-k_perp[j,i])))
                dist = k_base[l]-k_perp[j,i]
                if dist>0 :
                    spettro_f[w,l] += f_in[j,i,w]*(1-dist)
                    spettro_f[w,l-1] += f_in[j,i,w]*(dist)
                elif dist<0 :
                    spettro_f[w,l] += f_in[j,i,w]*(1-np.abs(dist))
                    spettro_f[w,l+1] += f_in[j,i,w]*(np.abs(dist))
                else:# (k_base - k[j,i])=0 :
                    spettro_f[w,l] += f_in[j,i,w]
    

    spettro_f_finale = np.delete(spettro_f,nx,axis=1)    

    end =time.time()
    print()
    print('END - gyrotropic_spectrum_3D')
    timing(start,end)
    
    return spettro_f_finale , kx , k_base

'''###############################################################################'''   
@jit(cache=True)
def spettri_ridotti(f_in,spettro_in):
    start=time.time()
    
    ny = len( f_in[0,:,0] ) ; nx = len( f_in[0,0,:] ) ; nz = len( f_in[:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i; kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    for i in range(nx,nz):
        kz[i]=-ny+i
            
    k_mod = np.zeros((nx,nx),dtype=float) #raggio sfera
    for i in range(nx):
        for j in range(nx):
            k_mod[j,i] = math.sqrt(kx[j]**2+kx[i]**2)
      
    k_base = k_mod[0,:]
    
    spettro_f_iso = np.zeros(nx+1,dtype=float) 
        
    for i in range(0,nx):
        for j in range(0,nx):      
            l = np.argmin(np.abs(k_base-k_mod[j,i]))
            dist = k_base[l]-k_mod[j,i]
            if dist>0 :
                spettro_f_iso[l] += spettro_in[j,i]*(1-dist)
                spettro_f_iso[l-1] += spettro_in[j,i]*(dist)
            elif dist<0 :
                spettro_f_iso[l] += spettro_in[j,i]*(1-np.abs(dist))
                spettro_f_iso[l+1] += spettro_in[j,i]*(np.abs(dist))
            else:# (k_base - k[j,i])=0 :
                spettro_f_iso[l] += spettro_in[j,i]
                                     

    spettro_f_finale_iso = np.delete(spettro_f_iso,[nx-1,nx])

    spettro_f_perp_rid =  np.zeros(nx,dtype=float) 
    spettro_f_parall_rid = np.zeros(nx,dtype=float) 
        
    for j in range(0,nx):
        for i in range(0,nx):      
            spettro_f_perp_rid[i] += spettro_in[j,i] 
            spettro_f_parall_rid[j] += spettro_in[j,i] 
                

    spettro_f_finale_perp = np.delete(spettro_f_perp_rid,[nx-1])
    spettro_f_finale_parall = np.delete(spettro_f_parall_rid,[nx-1])
        
    end =time.time()
    print()
    print('END - spettri_ridotti')
    timing(start,end)

    return spettro_f_finale_iso, spettro_f_finale_perp, spettro_f_finale_parall

'''###############################################################################'''
#@jit(cache=True)
def Curl(fx,fy,fz):
    start=time.time()
    
    ny = len( fx[0,:,0] ) ; nx = len( fx[0,0,:] ) ; nz = len( fx[:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
        if nz>nx: kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    if nz>nx:
        for i in range(nx,nz):
            kz[i]=-ny+i
    
    Jt_x = np.zeros((nz,ny,nx),dtype=complex)
    Jt_y = np.zeros((nz,ny,nx),dtype=complex)
    Jt_z = np.zeros((nz,ny,nx),dtype=complex)
    
    I=0+1j
    
    bxt = np.fft.rfftn(fx)
    byt = np.fft.rfftn(fy)
    bzt = np.fft.rfftn(fz)
    for w in range(nz):
        for j in range(ny):
            for i in range(nx):
                Jt_x[w,j,i] = ( ky[j]*bzt[w,j,i] - kz[w]*byt[w,j,i] )*I
                Jt_y[w,j,i] = ( kz[w]*bxt[w,j,i] - kx[i]*bzt[w,j,i] )*I
                Jt_z[w,j,i] = ( kx[i]*byt[w,j,i] - ky[j]*bxt[w,j,i] )*I
   
    J_x = np.fft.irfftn(Jt_x)
    J_y = np.fft.irfftn(Jt_y)
    J_z = np.fft.irfftn(Jt_z)
    
    end =time.time()
    print()
    print('END - Curl')
    timing(start,end)
    
    return J_x , J_y , J_z

'''###############################################################################'''
@jit(cache=True)
def Curl_migliore(fx,fy,fz): #secondo me questa funzione ha qualcosa di sbaglaito con nx,ny,nz
    start=time.time()
    
    ny = len( fx[0,:,0] ) ; nx = len( fx[0,0,:] ) ; nz = len( fx[:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
        if nz>nx: kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    if nz>nx:
        for i in range(nx,nz):
            kz[i]=-ny+i
    
    Jt_x1 = np.zeros((nz,ny,nx),dtype=complex)
    Jt_y1 = np.zeros((nz,ny,nx),dtype=complex)
    Jt_z1 = np.zeros((nz,ny,nx),dtype=complex)
    Jt_x2 = np.zeros((nz,ny,nx),dtype=complex)
    Jt_y2 = np.zeros((nz,ny,nx),dtype=complex)
    Jt_z2 = np.zeros((nz,ny,nx),dtype=complex)
        
    I=0+1j
    

    bxt = np.fft.rfftn(fx)
    byt = np.fft.rfftn(fy)
    bzt = np.fft.rfftn(fz)
    for w in range(nz):
        Jt_x2[w,:,:] = ( kz[w]*byt[w,:,:] ) 
        Jt_y1[w,:,:] = ( kz[w]*bxt[w,:,:] )
    for j in range(ny):
        Jt_x1[:,j,:] = ( ky[j]*bzt[:,j,:] )    
        Jt_z2[:,j,:] = ( ky[j]*bxt[:,j,:] )
    for i in range(nx):  
        Jt_y2[:,:,i] = ( kx[i]*bzt[:,:,i] )
        Jt_z1[:,:,i] = ( kx[i]*byt[:,:,i] )   
    J_x = np.fft.irfftn( (Jt_x1 - Jt_x2)*I )
    J_y = np.fft.irfftn( (Jt_y1 - Jt_y2)*I )
    J_z = np.fft.irfftn( (Jt_z1 - Jt_z2)*I )  
    
    end =time.time()
    print()
    print('END - Curl_migliore')
    timing(start,end)
    
    return J_x , J_y , J_z

'''###############################################################################'''
@jit(cache=True)
def Curl_migliorissimo(fx,fy,fz): #richiede troppa ram quindi sl mio pc non migliora quanto dovrebbe
    start=time.time()
    
    bxt = np.fft.rfftn(fx)
    byt = np.fft.rfftn(fy)
    bzt = np.fft.rfftn(fz)
    
    ny = len( bxt[0,:,0] ) ; nx = len( bxt[0,0,:] ) ; nz = len( bxt[:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
        if nz>nx: kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    if nz>nx:
        for i in range(nx,nz):
            kz[i]=-ny+i

    I=0-1j ; ''' COL MENO É GIUSTO '''
    
    J_x = np.fft.irfftn( (np.transpose( (ky*np.transpose(bzt,(0,2,1))) , (0,2,1) ) - np.transpose( (kz*np.transpose(byt,(2,1,0))) , (2,1,0) ))*I )
    J_y = np.fft.irfftn( (np.transpose( (kz*np.transpose(bxt,(2,1,0))) , (2,1,0) ) - kx*bzt)*I )
    J_z = np.fft.irfftn( (kx*byt - np.transpose( (ky*np.transpose(bxt,(0,2,1))) , (0,2,1) ))*I )
    
    end =time.time()
    print()
    print('END - Curl_migliorissimo')
    timing(start,end) 
    
    return J_x , J_y , J_z

'''###############################################################################'''
@jit(cache=True)
def Divergence(fx,fy,fz): #richiede troppa ram quindi sl mio pc non migliora quanto dovrebbe
    start=time.time()

    uxt = np.fft.rfftn(fx)
    uyt = np.fft.rfftn(fy)
    uzt = np.fft.rfftn(fz)
    
    ny = len( uxt[0,:,0] ) ; nx = len( uxt[0,0,:] ) ; nz = len( uxt[:,0,0] )
    kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)
    for i in range(0,nx):
        kx[i]=i ; ky[i]=i 
        if nz>nx: kz[i]=i
    for i in range(nx,ny):
        ky[i]=-ny+i   
    if nz>nx:
        for i in range(nx,nz):
            kz[i]=-ny+i
    
    I=0-1j  ; ''' COL MENO É GIUSTO '''
    
    div = np.fft.irfftn( (kx*uxt + np.transpose((ky*np.transpose(uyt,(0,2,1))),(0,2,1)) + np.transpose((kz*np.transpose(uzt,(2,1,0))),(2,1,0)))*I )
    
    end =time.time()
    print()
    print('END - Divergence')
    timing(start,end) 
    
    return div

'''###############################################################################'''
@jit(cache=True)
def Potenziale_Vettore(Jx,Jy,Jz): 
    start=time.time()
    
    jxt = (np.fft.rfftn(Jx))
    jyt = (np.fft.rfftn(Jy))
    jzt = (np.fft.rfftn(Jz))
    
    ny = len( jxt[0,:,0] ) ; nx = len( jxt[0,0,:] ) ; nz = len( jxt[:,0,0] )
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
    
    A_x = (np.fft.irfftn( jxt*invK2 ))
    A_y = (np.fft.irfftn( jyt*invK2 ))
    A_z = (np.fft.irfftn( jzt*invK2 )) 
    
    end =time.time()
    print()
    print('END - Potenziale_Vettore')
    timing(start,end) 
    
    return A_x , A_y , A_z

'''###############################################################################'''
@jit(cache=True)
def Potenziale_Vettore_daB(bx,by,bz): 
    start=time.time()

    bxt = (np.fft.rfftn(bx))
    byt = (np.fft.rfftn(by))
    bzt = (np.fft.rfftn(bz))
    
    ny = len( bxt[0,:,0] ) ; nx = len( bxt[0,0,:] ) ; nz = len( bxt[:,0,0] )
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

    I=0+1j
    
    A_x = (np.fft.irfftn( (np.transpose( (ky*np.transpose(bzt,(0,2,1))) , (0,2,1) ) - np.transpose( (kz*np.transpose(byt,(2,1,0))) , (2,1,0) ))*I*invK2 ))
    A_y = (np.fft.irfftn( (np.transpose( (kz*np.transpose(bxt,(2,1,0))) , (2,1,0) ) - kx*bzt)*I*invK2 ))
    A_z = (np.fft.irfftn( (kx*byt - np.transpose( (ky*np.transpose(bxt,(0,2,1))) , (0,2,1) ))*I*invK2 )) 
    
    end =time.time()
    print()
    print('END - Potenziale_Vettore_daB')
    timing(start,end) 
    
    return A_x , A_y , A_z

'''##############################################################################'''
@jit(cache=True)
def Varianza(f_in):
    start=time.time()
    
    ny = len( f_in[0,:,0] ) ; nx = len( f_in[0,0,:] ) ; nz = len( f_in[:,0,0] )
    f_somma=0. ; f_quadratico_somma=0.
    
    for w in range(nz):
        for j in range(ny):
            for i in range(nx):
                f_somma += f_in[w,j,i]
                f_quadratico_somma += (f_in[w,j,i])**2
    
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
def Pdf_2D(f_in):    
    start=time.time()

    ny = len( f_in[0,:,0] ) ; nx = len( f_in[0,0,:] )

    pdf_x = [] ; pdf_y = []
    for l in range(1,int(nx/2)+1):
        pdf_x = []
    for l in range(1,int(ny/2)+1):
        pdf_y = []

    for l in range(1,int(nx/2)+1):
        for i in range(nx):
            if i+l < nx :
                pdf_x[l-1] =  f_in[0,:,i+l] - f_in[0,:,i] 
            elif i+l >= nx :
                pdf_x[l-1] =  f_in[0,:,i+l-nx] - f_in[0,:,i]  
        pdf_x[l-1] = np.asarray(pdf_x[l-1]).reshape(-1)
    for l in range(1,int(ny/2)+1):
        for i in range(ny):
            if i+l < ny :
                pdf_y[l-1] =  f_in[0,i+1,:] - f_in[0,i,:]
            elif i+l >= nx :
                pdf_y[l-1] =  f_in[0,i+l-ny,:] - f_in[0,i,:]  
        pdf_y[l-1] = np.asarray(pdf_y[l-1]).reshape(-1)
        
    end=time.time()
    print()
    print('END - Pdf_2D')
    timing(start,end)
    
    return pdf_x , pdf_y

'''##############################################################################'''
@jit(cache=True) #questo funziona, ne son sicuro
def Pdf(f_in,flag):    
    start=time.time()

    ny = len( f_in[0,:,0] ) ; nx = len( f_in[0,0,:] ) ; nz = len( f_in[0,0,:] ) ; pdf=[]

    if flag=='x':
        n=nx
        pdf = []
        for l in range(1,int(n/2)+1):
            pdf = []    
                
        for l in range(1,int(n/2)+1):
            for i in range(n):
                if i+l < n :
                    pdf[l-1] =  f_in[:,:,i+l] - f_in[:,:,i] 
                elif i+l >= n :
                    pdf[l-1] =  f_in[:,:,i+l-n] - f_in[:,:,i]
            pdf[l-1] = np.asarray(pdf[l-1]).reshape(-1)
        
    elif flag=='y': 
        n=ny
        pdf = []
        for l in range(1,int(n/2)+1):
            pdf = []
        
        for l in range(1,int(n/2)+1):
            for i in range(n):
                if i+l < n :
                    pdf[l-1] =  f_in[:,i+l,:] - f_in[:,i,:]
                elif i+l >= n :
                    pdf[l-1] =  f_in[:,i+l-n,:] - f_in[:,i,:]
            pdf[l-1] = np.asarray(pdf[l-1]).reshape(-1)
    
    elif flag=='z': 
        n=nz
        pdf = []
        for l in range(1,int(n/2)+1):
            pdf = []
        
        for l in range(1,int(n/2)+1):
            for i in range(n):
                if i+l < n :
                    pdf[l-1] =  f_in[i+l,:,:] - f_in[i,:,:]
                elif i+l >= n :
                    pdf[l-1] =  f_in[i+l-n,:,:] - f_in[i,:,:]
            pdf[l-1] = np.asarray(pdf[l-1]).reshape(-1)
                
    end=time.time()
    print()
    print('END - Pdf')
    timing(start,end)
    
    return pdf

'''##############################################################################'''
@jit(cache=True)
def mom_2D_miei(pdf_in,bins): #sono più lenti e meno precisi della funzione di python 
    start=time.time()
    
    hist = np.histogram(pdf_in,bins=bins)[0]
    base_hist = np.histogram(pdf_in,bins=bins)[1]

    mom=np.zeros(5,dtype=float)
    for j in range(5):
        for i in range(len(hist)):
            mom[j] += hist[i]*(base_hist[i+1]-base_hist[i])*( ((base_hist[i+1]+base_hist[i])*0.5)**j ) 
    mom[:] = mom[:]/mom[0]
    
    end =time.time()
    print()
    print('END - mom_2D_miei')
    timing(start,end)
    
    return mom

'''##############################################################################'''
@jit(cache=True)
def momenti(pdf_in):
    start=time.time()
    nx = 2*len( pdf_in[:] )
            
    mom = np.zeros((int(nx/2),5),dtype=float)
    for l in range(int(nx/2)):
        for j in range(5):
            mom[l,j] = stats.moment(pdf_in[l],moment=j)
        mom[l,:] = mom[l,:]/mom[l,0]
    
    K = mom[:,4]/((mom[:,2])**2) - 3 
    S = mom[:,3]/((mom[:,2])**1.5)
    
    end =time.time()
    print()
    print('END - momenti')
    timing(start,end)
    
    return mom , K , S 

'''##############################################################################'''
@jit(cache=True)
def momenti_distribuzione(arr):
    start=time.time()
            
    mom = np.zeros(5,dtype=float)
    for j in range(5):
        mom[j] = stats.moment(arr,moment=j)
    mom[:] = mom[:]/mom[0]
    
    K = mom[4]/((mom[2])**2) - 3 
    S = mom[3]/((mom[2])**1.5)
    
    end =time.time()
    print()
    print('END - momenti_distribuzione')
    timing(start,end)
    
    return mom , K , S 
