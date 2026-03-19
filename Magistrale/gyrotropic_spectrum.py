# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 10:08:48 2019

@author: Mattia
"""
import funzioni_matti as matti
import numpy as np ; import h5py ;  import pylab as plt ; import matplotlib as mplt ; import matplotlib.cm as cm ; import time 

start = time.time()
B0=1.0  #solo lungo x 
tempi=40

LF=True

rh=[]; pe=[]; RHfM2=[]; PEfM2=[]; t_out=np.zeros(tempi+1,dtype=float)
bx=[]; by=[]; bz=[]; bzfM2=[]; byfM2=[]; bxfM2=[]; bxfM2_noB0=[]
vx=[]; vy=[]; vz=[]; vzfM2=[]; vyfM2=[]; vxfM2=[]
BfM2=[]; BfM2_noB0=[]; VfM2=[] ; const_matrix=[]       #MODULO QUADRO DEL B TOTALE (come somma dei moduli quadri di bx, by e bz)

directory='/Users/Mattia/Desktop/Risultati L.Franci/256^3(init_onda_2D. LF. 256tempi)(Lz=64Lx)/'
#directory='/Users/Mattia/Desktop/risultati SLandi/128x128x128(init_onda_2D-- 2 onde, una x-y l altra y-z[tmax=4.4; 45°(1,-1) e bo (2,3)])/'
#directory='/Users/Mattia/Desktop/risultati SLandi/256x256x256(init_turb-- 3D)/'

o=1
for i in range(tempi+1):
    if o*i<10: c = str(o*i) ; a = '0'+c
    else: a = str(o*i)
    print(a)
    h=h5py.File(directory+'out%s0.h5' %a)
    t_out[i] = ( h.get('time')[...] )
    bx.append( h.get('bx')[...] )
    by.append( h.get('by')[...] )
    bz.append( h.get('bz')[...] )
#    bxfM2.append( (np.absolute(np.fft.rfftn(bx[i])))**2 ) #così ho [x,y,z] mi serve poichè ho b0 lungo z e la mia routine considera b0 lungo x
#    byfM2.append( (np.absolute(np.fft.rfftn(by[i])))**2 )
#    bzfM2.append( (np.absolute(np.fft.rfftn(bz[i])))**2 )
    bxfM2.append( (np.absolute(np.fft.rfftn(np.transpose(bx[i],(2,1,0)))))**2 ) #così ho [x,y,z] mi serve poichè ho b0 lungo z e la mia routine considera b0 lungo x
    byfM2.append( (np.absolute(np.fft.rfftn(np.transpose(by[i],(2,1,0)))))**2 )
    bzfM2.append( (np.absolute(np.fft.rfftn(np.transpose(bz[i],(2,1,0)))))**2 )
    BfM2.append( bxfM2[i]+byfM2[i]+bzfM2[i] )
#    vx.append( h.get('vx')[...] )
#    vy.append( h.get('vy')[...] )
#    vz.append( h.get('vz')[...] )
#    vxfM2.append( (np.absolute(np.fft.rfftn(vx[i])))**2 )
#    vyfM2.append( (np.absolute(np.fft.rfftn(vy[i])))**2 )
#    vzfM2.append( (np.absolute(np.fft.rfftn(vz[i])))**2 )
#    VfM2.append( vxfM2[i]+vyfM2[i]+vzfM2[i] )    
#    bxfM2_noB0.append( (np.absolute(np.fft.rfft2(bx[i]-B0)))**2 )
#    BfM2_noB0.append( bxfM2_noB0[i]+byfM2[i]+bzfM2[i] )    
#    rh.append( h.get('rh')[...] )      ; pe.append( h.get('pe')[...] )
#    RHfM2.append( (np.absolute(np.fft.rfftn(rh[i])))**2 )        ; PEfM2.append( (np.absolute(np.fft.rfftn(pe[i])))**2 )
    
#for i in range(tempi+1):
#    const_matrix.append( np.zeros((nz,ny,nx))+1.)
    
ny = len( bzfM2[0][0,:,0] ) ; nx = len( bzfM2[0][0,0,:] ) ; nz = len( bzfM2[0][:,0,0] )
kx = np.zeros(nx, dtype=int) ; ky = np.zeros(ny, dtype=int) ; kz = np.zeros(nz, dtype=int)

for i in range(0,nx):
    kx[i]=i ; ky[i]=i 
    if nz>nx: kz[i]=i
for i in range(nx,ny):
    ky[i]=-ny+i   
if nz>nx:
    for i in range(nx,nz):
        kz[i]=-ny+i

end = time.time()

matti.timing(start,end)

if LF==True :
    kx = np.zeros(nx, dtype=float) 
    for i in range(0,nx):
        kx[i]=i/64
        
#%%

out_file = open("tempi_per_output.dat","w")
for i in range(tempi+1):
    if i<10: c = str(i) ; a = '0'+c
    else: a = str(i)
    out_file.write(t_out[i])
out_file.close()

#%%
#radial_spectrum_2D(const_matrix,spettro_const_matrix,k,tempi)
#y = math.pi*k   #deve tornare così dato che sto integrando "1" su una semicirconferenza
#plt.figure(0) #ne ploto uno solo.. tanto son tutti uguali
#plt.plot(k[0:nx-1],spettro_const_matrix[0])
#plt.plot(k[0:nx-1],y[0:nx-1])
#plt.title('spettro matrice costante (%s)' %a)
#
##%%
#gyrotropic_spectrum_3D(const_matrix,spettro_const_matrix,k,tempi)
##y = math.pi*k   #deve tornare così dato che sto integrando "1" su una semicirconferenza
#plt.figure(1) #ne ploto uno solo.. tanto son tutti uguali
#plt.plot(k[0:nx-1],spettro_const_matrix[0])
##plt.plot(k[0:nx-1],y[0:nx-1])
#plt.title('spettro girotropico 3D matrice costante con pesi')
#    
#%%
'''PROBABILITY DENSITY FUNCTION ########################################################'''
pdfx_bx = matti.Pdf(bx,tempi,'x')
pdfy_bx = matti.Pdf(bx,tempi,'y')
#%%
momx_bx , kurt_x , skew_x = matti.momenti(pdfx_bx,tempi)
momy_bx , kurt_y , skew_y = matti.momenti(pdfy_bx,tempi)
#bins=100000
#mom_bx_miei = mom_2D_miei(pdf_bx,tempi,bins)
#%% kurtosis e skewness a una scala
plt.figure(0) ; plt.plot(t_out,kurt_x[:,250]) ; plt.plot(t_out,np.zeros(len(t_out))) ; plt.title('Kurtosis') 
plt.figure(1) ; plt.plot(t_out,skew_x[:,250]) ; plt.plot(t_out,np.zeros(len(t_out))) ; plt.title('Skewness') 
#%%
plt.figure(0); plt.title('Kurtosis lungo x') 
for l in range(int(nx/16)): #così plotto una scala ogni 8
    plt.plot(t_out,kurt_x[:,int(l*8)])
plt.figure(1); plt.title('Kurtosis lungo y') 
for l in range(int(nx/16)): #così plotto una scala ogni 8
    plt.plot(t_out,kurt_y[:,int(l*8)])
#%%
plt.figure(2); plt.title('Skewness lungo x') 
for l in range(int(nx/16)): #così plotto una scala ogni 8
    plt.plot(t_out,skew_x[:,int(l*8)])
plt.figure(3); plt.title('Skewness lungo y') 
for l in range(int(nx/16)): #così plotto una scala ogni 8
    plt.plot(t_out,skew_y[:,int(l*8)])
#%%
for t in range(tempi+1):
    a=str(t)
    plt.figure(t)
    plt.hist(pdfx_bx[t][250],bins=100,log=True)
    plt.title('distribuzione delta bx(%s)'%a)
    
#%%
'''DENSITà DI CORRENTE ########################################################'''
J_x , J_y , J_z = matti.Curl_migliorissimo(bx,by,bz,tempi)
max_Jx = [] ; max_Jy =[] ; max_Jz =[]
for t in range(tempi+1):
    max_Jx.append(np.amax(J_x[t])) ; max_Jy.append(np.amax(J_y[t])) ; max_Jz.append(np.amax(J_z[t]))
#%% provo a stampare su file h5 le correnti:
for t in range(tempi+1):
    if t<10: c = str(t) ; a = '0'+c
    else: a = str(t)
    hf = h5py.File('../outJ0%s.h5' %a,'w')
    hf.create_dataset('jx', data=np.float32(J_x[t]))
    hf.create_dataset('jy', data=np.float32(J_y[t]))
    hf.create_dataset('jz', data=np.float32(J_z[t]))
    hf.close()    
    
hf = h5py.File('../outJ000.h5','r') ; d=0
for k in hf.keys():
    d=d+1
    print('dataset %s ='%d,k)
hf.close()

#%%
plt.figure(100) ; plt.plot(t_out,np.sqrt(max_Jx)) ; plt.title('massimi Jx') 
plt.figure(101) ; plt.plot(t_out,np.sqrt(max_Jy)) ; plt.title('massimi Jy') 
plt.figure(102) ; plt.plot(t_out,np.sqrt(max_Jz)) ; plt.title('massimi Jz') 

#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(J_x[i][0,:,:],levels=np.linspace(-2.7,2.7,255),cmap=cm.jet)
    plt.colorbar() ; plt.title('J(x) (%s)' %a)
    
#%%
'''RMS J #####'''
Var_Jx = matti.Varianza(J_x,tempi)
Var_Jy = matti.Varianza(J_y,tempi)
Var_Jz = matti.Varianza(J_z,tempi)
Var_Jtot = Var_Jx + Var_Jy + Var_Jz
#%%
plt.figure(0) ; plt.plot(t_out,np.sqrt(Var_Jx)) ; plt.title('Varianza Jx') 
plt.figure(1) ; plt.plot(t_out,np.sqrt(Var_Jy)) ; plt.title('Varianza Jy') 
plt.figure(2) ; plt.plot(t_out,np.sqrt(Var_Jz)) ; plt.title('Varianza Jz') 
plt.figure(3) ; plt.plot(t_out,np.sqrt(Var_Jtot)) ; plt.title('Varianza J totale')
#%%
'''CAMPO MAGNETICO (2D)########################################################'''
spettro_B , k = matti.radial_spectrum_2D(BfM2,tempi)
#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.loglog(k[0:nx-1],spettro_B[i]*k[0:nx-1]**2)
    plt.ylim((10**(-7),10**13))
    plt.title('spettro B(%s) modulo quadro con pesi' %a)
    
#%%
spettro_B_noB0 , k = matti.radial_spectrum_2D(BfM2_noB0,tempi)
#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.loglog(k[0:nx-1],spettro_B_noB0[i])
    plt.ylim((10**(-9),10**9))
    plt.title('spettro B(%s) [senza B0] modulo quadro con pesi' %a)
    
#%%  
'''CAMPO MAGNETICO ########################################################'''
spettro_B_gyro = matti.gyrotropic_spectrum_3D(BfM2,tempi)
#%% parallel = asse x . ma in questo caso, avendo trasposto all'inizio, in realta parallel = asse z
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(spettro_B_gyro[i],cmap=cm.jet,levels=np.logspace(-7, 14, 255),locator=mplt.ticker.LogLocator())
    mplt.pyplot.yscale('log') ; mplt.pyplot.xscale('log') ; plt.ylim((1,nx-1)) ; plt.xlim((1,nx-1))
    plt.ylabel('k parallel') ; plt.xlabel('k perpendicular') ; plt.colorbar() ; plt.title('Spettro B gyrotropico (%s)' %a)

#%%
spettro_B_iso, spettro_B_perp, spettro_B_parall = matti.spettri_ridotti(BfM2,spettro_B_gyro,tempi)

#%%
alpha=0.
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.loglog(ky[0:nx-1],spettro_B_perp[i]*(ky[0:nx-1])**(alpha),label='Spettro perpendicolare')
    plt.loglog(kx[0:nx-1],spettro_B_parall[i]*(kx[0:nx-1])**(alpha),label='Spettro parallelo')
#    plt.loglog(kx[0:nx-1],spettro_B_iso[i],label='Spettro isotropico')
    mplt.pyplot.legend()
    plt.ylim((10**-5,10**13))

#%%





























#%%
    '''VELOCITà ################################################################'''
spettro_V_gyro = matti.gyrotropic_spectrum_3D(VfM2,tempi)
#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(spettro_V_gyro[i],cmap=cm.jet,levels=np.logspace(-1, 13, 255),locator=mplt.ticker.LogLocator())
    mplt.pyplot.yscale('log') ; mplt.pyplot.xscale('log') ; plt.ylim((1,nx-1)) ; plt.xlim((1,nx-1))
    plt.ylabel('k parallel') ; plt.xlabel('k perpendicular') ; plt.colorbar() ; plt.title('Spettro V gyrotropico (%s)' %a)

#%%
spettro_V_iso, spettro_V_perp, spettro_V_parall = matti.spettri_ridotti(VfM2,spettro_V_gyro,tempi)
#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.loglog(k[0:nx-1],spettro_V_perp[i])  
    plt.loglog(k[0:nx-1],spettro_V_parall[i])         
    plt.loglog(k[0:nx-1],spettro_V_iso[i])
    plt.title("'V': blue(K_perp) - orange(K_parallel) - green(K_isotropic)")
    plt.ylim((10**-1,10**13))

#%%
    '''DENSITà ###############################################################'''
spettro_RH_gyro = matti.gyrotropic_spectrum_3D(RHfM2,tempi)
#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(spettro_RH_gyro[i],cmap=cm.jet,levels=np.logspace(-1, 13, 255),locator=mplt.ticker.LogLocator())
    mplt.pyplot.yscale('log') ; mplt.pyplot.xscale('log') ; plt.ylim((1,nx-1)) ; plt.xlim((1,nx-1))
    plt.ylabel('k parallel') ; plt.xlabel('k perpendicular') ; plt.colorbar() ; plt.title('Spettro RH gyrotropico (%s)' %a)

#%%
spettro_RH_iso, spettro_RH_perp, spettro_RH_parall = matti.spettri_ridotti(RHfM2,spettro_RH_gyro,tempi)
#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.loglog(k[0:nx-1],spettro_RH_perp[i])  
    plt.loglog(k[0:nx-1],spettro_RH_parall[i])         
    plt.loglog(k[0:nx-1],spettro_RH_iso[i])
    plt.title("'RH': blue(K_perp) - orange(K_parallel) - green(K_isotropic)")
    plt.ylim((10**-1,10**13))
#%%
    '''PRESSIONE #############################################################'''
spettro_PE_gyro = matti.gyrotropic_spectrum_3D(PEfM2,tempi)
#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.contourf(spettro_PE_gyro[i],cmap=cm.jet,levels=np.logspace(-1, 13, 255),locator=mplt.ticker.LogLocator())
    mplt.pyplot.yscale('log') ; mplt.pyplot.xscale('log') ; plt.ylim((1,nx-1)) ; plt.xlim((1,nx-1))
    plt.ylabel('k parallel') ; plt.xlabel('k perpendicular') ; plt.colorbar() ; plt.title('Spettro PE gyrotropico (%s)' %a)
#%%
spettro_PE_iso, spettro_PE_perp, spettro_PE_parall = matti.spettri_ridotti(PEfM2,spettro_PE_gyro,tempi)
#%%
for i in range(tempi+1):
    a = str(i)
    plt.figure(i)
    plt.loglog(k[0:nx-1],spettro_PE_perp[i])  
    plt.loglog(k[0:nx-1],spettro_PE_parall[i])         
    plt.loglog(k[0:nx-1],spettro_PE_iso[i])
    plt.title("'PE': blue(K_perp) - orange(K_parallel) - green(K_isotropic)")
    plt.ylim((10**-1,10**13))