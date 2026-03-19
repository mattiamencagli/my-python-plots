# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:09:16 2019

@author: Mattia
"""
import funzioni_matti as matti ; from scipy import stats
import numpy as np ; import h5py ;  import pylab as plt ; import matplotlib as mplt ; import matplotlib.cm as cm ; import time 

start=time.time()

cartella='/Users/Mattia/Desktop/risultati_LDZ/'

#directory=cartella+'512^3 (init_turb. va=0.8 tmax=16)/'             ;tf=40   ; threeD  =  True  ; post=False
directory=cartella+'512^3 (init_turb. va=0.6 tmax=40)/'             ;tf=159  ; threeD  =  True  ; post=False ; nome=r'$F_3$'
#directory=cartella+'512^3 (init_turb. va=0.6 tmax=40)/'             ;tf=159  ; threeD  =  True  ; post=True ; nome=r'$F_3$'  ; NNN=512 ; dd=3
#directory=cartella+'256^3 (init_turb. va=0.6 tmax=40)/'             ;tf=160  ; threeD  =  True  ; post=False ; nome=r'$E_3$'
#directory=cartella+'256^3 (init_turb. va=0.6 tmax=40)/'             ;tf=160  ; threeD  =  True  ; post=True ; nome=r'$E_3$'  ; NNN=256 ; dd=3
#directory=cartella+'2048^2 (init_turb. va=0.5 tmax=36)/'            ;tf=73   ; threeD  =  False ; post=False
#directory=cartella+'2048^2 (init_turb. va=0.87 tmax=4.5)/'          ;tf=16   ; threeD  =  False ; post=False ; nome=r'$B_2$' ; NNN=2048 ; dd=2
#directory=cartella+'2048^2 (init_turb. va=0.75 tmax=6)/'            ;tf=23   ; threeD  =  False ; post=False
#directory=cartella+'2048^2 (init_turb. va=0.6 tmax=25)/'            ;tf=98   ; threeD  =  False ; post=False ; nome=r'$A_2$' ; NNN=2048 ; dd=2
#directory=cartella+'2048^2 (init_turb. va=0.6 tmax=25)(fitta)/'     ;tf=43   ; threeD  =  False ; post=False
#directory=cartella+'4096^2 (init_turb. va=0.65 tmax=5)/'            ;tf=20   ; threeD  =  False ; post=False ; nome=r'$D_2$' ; NNN=4096 ; dd=2

#directory=cartella+'2048^2 (init_turb_rms. va=0.57 tmax=20)/'       ;tf=85   ; threeD  =  False ; post=False ; nome=r'$C_2$' ; NNN=2048 ; dd=2
#directory=cartella+'512^3 (init_turb_rms. va=0.18 tmax=40)/'        ;tf=160  ; threeD  =  True  ; post=False
  
del(cartella)

d=1         #ogni quanti tempi leggo?
ti=0        #tempo iniziale

save    =  True        #salva i plot?
wid     =  1.1         #larghezza linee dei plot

tempi = int((tf-ti)/d) 

rms_rh=np.zeros(tempi+1,dtype=np.float32) ; rms_pe=np.zeros(tempi+1,dtype=np.float32) ; max_rh=np.zeros(tempi+1,dtype=np.float32)
max_pe=np.zeros(tempi+1,dtype=np.float32) ; min_rh=np.zeros(tempi+1,dtype=np.float32) ; min_pe=np.zeros(tempi+1,dtype=np.float32)
max_Jx=np.zeros(tempi+1,dtype=np.float32) ; max_Jy=np.zeros(tempi+1,dtype=np.float32) ; max_Jz=np.zeros(tempi+1,dtype=np.float32) 
max_vx=np.zeros(tempi+1,dtype=np.float32) ; max_vy=np.zeros(tempi+1,dtype=np.float32) ; max_vz=np.zeros(tempi+1,dtype=np.float32)
max_bx=np.zeros(tempi+1,dtype=np.float32) ; max_by=np.zeros(tempi+1,dtype=np.float32) ; max_bz=np.zeros(tempi+1,dtype=np.float32)
rms_bx=np.zeros(tempi+1,dtype=np.float32) ; rms_by=np.zeros(tempi+1,dtype=np.float32) ; rms_bz=np.zeros(tempi+1,dtype=np.float32)
rms_vx=np.zeros(tempi+1,dtype=np.float32) ; rms_vy=np.zeros(tempi+1,dtype=np.float32) ; rms_vz=np.zeros(tempi+1,dtype=np.float32)
rms_Jx=np.zeros(tempi+1,dtype=np.float32) ; rms_Jy=np.zeros(tempi+1,dtype=np.float32) ; rms_Jz=np.zeros(tempi+1,dtype=np.float32)
rms_vorx=np.zeros(tempi+1,dtype=np.float32) ; rms_vory=np.zeros(tempi+1,dtype=np.float32) ; rms_vorz=np.zeros(tempi+1,dtype=np.float32)
rms_Jtot=np.zeros(tempi+1,dtype=np.float32) ; rms_btot=np.zeros(tempi+1,dtype=np.float32) ; rms_vtot=np.zeros(tempi+1,dtype=np.float32) 
rms_vortot=np.zeros(tempi+1,dtype=np.float32) ; rms_Vatot=np.zeros(tempi+1,dtype=np.float32)

helicityJ=np.zeros(tempi+1,dtype=np.float32) ; helicityB=np.zeros(tempi+1,dtype=np.float32) 
helicityJ_g=np.zeros(tempi+1,dtype=np.float32) ; helicityB_g=np.zeros(tempi+1,dtype=np.float32) 
chw=np.zeros(tempi+1,dtype=np.float32) ; ch=np.zeros(tempi+1,dtype=np.float32) 
    
#sbagliati
#ava_ch=np.zeros(tempi+1,dtype=np.float32) ; ava_chro=np.zeros(tempi+1,dtype=np.float32) 
#rms_ch=np.zeros(tempi+1,dtype=np.float32) ; rms_chro=np.zeros(tempi+1,dtype=np.float32)

''' ROBA ALQUANTO INUTILE... '''
#rms_zpx=np.zeros(tempi+1,dtype=np.float32) ; rms_zpy=np.zeros(tempi+1,dtype=np.float32) ; rms_zpz=np.zeros(tempi+1,dtype=np.float32)
#rms_zmx=np.zeros(tempi+1,dtype=np.float32) ; rms_zmy=np.zeros(tempi+1,dtype=np.float32) ; rms_zmz=np.zeros(tempi+1,dtype=np.float32)
#rms_zptot=np.zeros(tempi+1,dtype=np.float32) ; rms_zmtot=np.zeros(tempi+1,dtype=np.float32)
#Crohel_ava=np.zeros(tempi+1,dtype=np.float32) ; Crohel_ava2=np.zeros(tempi+1,dtype=np.float32) ; Crohel_rms=np.zeros(tempi+1,dtype=np.float32) 
#rms_mod_zp=np.zeros(tempi+1,dtype=np.float32) ; rms_mod_zm=np.zeros(tempi+1,dtype=np.float32)
#ava_zp=np.zeros(tempi+1,dtype=np.float32) ; ava_zm=np.zeros(tempi+1,dtype=np.float32)
#ava_zp2=np.zeros(tempi+1,dtype=np.float32) ; ava_zm2=np.zeros(tempi+1,dtype=np.float32)

t_out=np.zeros(tempi+1,dtype=np.float32) 

if threeD==True : 
    spettro_B_gyro=[] ; spettro_B_iso=[] ; spettro_B_perp=[] ; spettro_B_parall=[]
    spettro_V_gyro=[] ; spettro_V_iso=[] ; spettro_V_perp=[] ; spettro_V_parall=[]
    spettro_Etot_perp=[] ; spettro_Etot_parall=[] ; spettro_Etot_gyro=[]
    K_para=[] ; K_perp=[] 
else :   
    spettro_B_rad=[] ; spettro_V_rad=[] ; spettro_Etot_rad=[]
    
if post==False:
    for i in range(tempi+1):
        
        if d*i+ti<1000: a = str(d*i+ti)
        if d*i+ti<100:  c = str(d*i+ti) ; a = '0'+c
        if d*i+ti<10:   c = str(d*i+ti) ; a = '00'+c
        print('out_analysis%s.h5' %a)
        
        h=h5py.File(directory+'analisi/out_analysis%s.h5' %a)

        max_rh[i] = (h.get('max_rh')[...] )
        max_pe[i] = (h.get('max_pe')[...] )
        min_rh[i] = (h.get('min_rh')[...] )
        min_pe[i] = (h.get('min_pe')[...] )
        max_vx[i] = (h.get('max_vx')[...] )
        max_vy[i] = (h.get('max_vy')[...] )
        max_vz[i] = (h.get('max_vz')[...] )
        max_bx[i] = (h.get('max_bx')[...] )
        max_by[i] = (h.get('max_by')[...] )
        max_bz[i] = (h.get('max_bz')[...] )
        max_Jx[i] = (h.get('max_Jx')[...] )
        max_Jy[i] = (h.get('max_Jy')[...] )
        max_Jz[i] = (h.get('max_Jz')[...] )
        
        ''' IN REALTà QUESTE SONO LE VARIANZE ---> INFATTI NEI PLOT POI FAI LA RADICE'''
        ''' QUINDI SOMMARLI E BASTA VA BENE DATO CHE FACCIO LA RADICE ... NON DIMENTICARLO PER LA SETTORDICESIMA VOLTA PLEASE'''
        rms_rh[i] = (h.get('rms_rh')[...] )
        rms_pe[i] = (h.get('rms_pe')[...] )
        rms_bx[i] = (h.get('rms_bx')[...] )
        rms_by[i] = (h.get('rms_by')[...] )
        rms_bz[i] = (h.get('rms_bz')[...] )
        rms_btot[i] = (rms_bx[i] + rms_by[i] + rms_bz[i])
        rms_vx[i] = (h.get('rms_vx')[...] )
        rms_vy[i] = (h.get('rms_vy')[...] )
        rms_vz[i] = (h.get('rms_vz')[...] )
        rms_vtot[i] = (rms_vx[i] + rms_vy[i] + rms_vz[i])
        rms_Jx[i] = (h.get('rms_Jx')[...] )
        rms_Jy[i] = (h.get('rms_Jy')[...] )
        rms_Jz[i] = (h.get('rms_Jz')[...] )
        rms_Jtot[i] = (rms_Jx[i] + rms_Jy[i] + rms_Jz[i])
        rms_vorx[i] = (h.get('rms_vorx')[...] )
        rms_vory[i] = (h.get('rms_vory')[...] )
        rms_vorz[i] = (h.get('rms_vorz')[...] )
        rms_vortot[i] = (rms_vorx[i] + rms_vory[i] + rms_vorz[i])
        
        if  threeD==False :
            helicityJ[i] = (h.get('helicityJ')[...])
            helicityB[i] = (h.get('helicityB')[...])
            helicityJ_g[i] = (h.get('helicityJ_giusto')[...])
            helicityB_g[i] = (h.get('helicityB_giusto')[...])
            chw[i] = (h.get('chw_giusto')[...])
            ch[i] = (h.get('ch_giusto')[...])
            #sbaglaiti credo:
    #        ava_ch[i] = (h.get('ava_ch')[...])
    #        ava_chro[i] = (h.get('ava_chro')[...])
    #        rms_ch[i] = (h.get('rms_ch')[...])
    #        rms_chro[i] = (h.get('rms_chro')[...])
    
        
        ''' ROBA ALQUANTO INUTILE... '''
    #    rms_zpx[i] = (h.get('rms_zpx')[...] )
    #    rms_zpy[i] = (h.get('rms_zpy')[...] ) 
    #    rms_zpz[i] = (h.get('rms_zpz')[...] ) 
    #    rms_zmx[i] = (h.get('rms_zmx')[...] )
    #    rms_zmy[i] = (h.get('rms_zmy')[...] )
    #    rms_zmz[i] = (h.get('rms_zmz')[...] )
    #    rms_zptot[i] = rms_zpx[i] + rms_zpy[i] + rms_zpz[i]
    #    rms_zmtot[i] = rms_zmx[i] + rms_zmy[i] + rms_zmz[i] 
    #    rms_mod_zp[i] = (h.get('rms_mod_zp')[...] ) 
    #    rms_mod_zm[i] = (h.get('rms_mod_zm')[...] )
    #    ava_zp[i] = (h.get('ava_zp')[...] )
    #    ava_zm[i] = (h.get('ava_zm')[...] )
    #    ava_zp2[i] = (h.get('ava_zp')[...] )
    #    ava_zm2[i] = (h.get('ava_zm')[...] )
    #    Crohel_ava[i] = ( ava_zp[i]**2 - ava_zm[i]**2 )/( ava_zp[i]**2 + ava_zm[i]**2 )
    #    Crohel_rms[i] = ( rms_mod_zp[i] - rms_mod_zm[i] )/( rms_mod_zp[i] + rms_mod_zm[i] )
        
        t_out[i] = h.get('t_out')[...]
    
        
        if threeD==True :
            spettro_B_gyro.append( h.get('spettro_B_gyro')[...] )
            spettro_B_iso.append( h.get('spettro_B_iso')[...] )
            spettro_B_perp.append( h.get('spettro_B_perp')[...] )
            spettro_B_parall.append( h.get('spettro_B_parall')[...] )
            spettro_V_gyro.append( h.get('spettro_V_gyro')[...] )
            spettro_V_iso.append( h.get('spettro_V_iso')[...] )
            spettro_V_perp.append( h.get('spettro_V_perp')[...] )
            spettro_V_parall.append( h.get('spettro_V_parall')[...] )
            spettro_Etot_gyro.append( spettro_B_gyro[i] + spettro_B_gyro[i] )
            spettro_Etot_perp.append( spettro_B_perp[i] + spettro_V_perp[i] )
            spettro_Etot_parall.append( spettro_B_parall[i] + spettro_V_parall[i] )
            kperp  = h.get('kperp')[...]
            kparal = h.get('kparal')[...]
            K_para.append(spettro_Etot_gyro[i][:,0])
            K_perp.append(spettro_Etot_gyro[i][0,:])
            
        else :
            spettro_B_rad.append( h.get('spettro_B_rad')[...] )
            spettro_V_rad.append( h.get('spettro_V_rad')[...] )
            spettro_Etot_rad.append(spettro_B_rad[i] + spettro_V_rad[i])
            k = h.get('k')[...] 
    
        h.close()
    
else:
    for i in range(tempi+1):
        
        if d*i+ti<1000: a = str(d*i+ti)
        if d*i+ti<100:  c = str(d*i+ti) ; a = '0'+c
        if d*i+ti<10:   c = str(d*i+ti) ; a = '00'+c
        print('out_analysis%s.h5' %a)
        
        h=h5py.File(directory+'analisi_post/out_analysis%s.h5' %a)
        
        helicityJ_g[i] = (h.get('helicityJ')[...])
        helicityB_g[i] = (h.get('helicityB')[...])
        chw[i] = (h.get('chw_giusto')[...])
        ch[i] = (h.get('ch_giusto')[...])
        rms_Jx[i] = (h.get('rms_Jx')[...] )
        rms_Jy[i] = (h.get('rms_Jy')[...] )
        rms_Jz[i] = (h.get('rms_Jz')[...] )
        rms_Jtot[i] = np.sqrt(rms_Jx[i]*rms_Jx[i] + rms_Jy[i]*rms_Jy[i] + rms_Jz[i]*rms_Jz[i])
        rms_vorx[i] = (h.get('rms_vorx')[...] )
        rms_vory[i] = (h.get('rms_vory')[...] )
        rms_vorz[i] = (h.get('rms_vorz')[...] )
        rms_vortot[i] = np.sqrt(rms_vorx[i]*rms_vorx[i] + rms_vory[i]*rms_vory[i] + rms_vorz[i]*rms_vorz[i])
        max_Jx[i] = (h.get('max_Jx')[...] )
        max_Jy[i] = (h.get('max_Jy')[...] )
        max_Jz[i] = (h.get('max_Jz')[...] )
        t_out[i] = h.get('t_out')[...]
        h.close()
        
end=time.time()
matti.timing(start,end)


 #%% spettro_B_gyro
nx=len(spettro_B_gyro[0][0,:])
vmin=0 ; vmax=14
for i in [49]:
    a = str(i)
    plt.figure(i)
    title='Spettro_B_gyrotropico_max_turb(%s)' %a
    plt.contourf(spettro_B_gyro[i],cmap=cm.jet,levels=np.logspace(vmin,vmax, 255),locator=mplt.ticker.LogLocator())
    mplt.pyplot.yscale('log') ; mplt.pyplot.xscale('log') 
    plt.ylim((1,nx-1)) ; plt.xlim((1,nx-1)) 
    plt.ylabel('k parallel') ; plt.xlabel('k perpendicular') 
    plt.colorbar(ticks=np.logspace(vmin,vmax,vmax-vmin+1)) 
    plt.title(title) 
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')
        
#%% B ridotto
alpha=5/3
beta=3/2
gamma=7/2
#alpha=0 ; beta=0 ; gamma=0
nx=len(spettro_B_perp[0]) 
kolmogorov=np.zeros(nx,dtype=float)
for i in range(nx):
    kolmogorov[i]=1*10**15

plt.figure(25)
#plt.vlines(3, 10**9, 10**16,linestyles='dotted',colors='xkcd:purple',linewidth=wid)
#plt.vlines(nx -1 -3, 10**9, 10**16,linestyles='dotted',colors='xkcd:purple',linewidth=wid)
#plt.loglog(kperp[0:nx],kolmogorov,'xkcd:black',linewidth=wid)
#plt.loglog(kperp[0:nx],kolmogorov*6,'xkcd:black',linewidth=wid)
title='Spettri_ridotti_B_max_turb,alpha=%1.2f,beta=%1.2f,gamma=%1.2f'%(alpha,beta,gamma)
plt.ylim((10**14,7*10**15))
for i in [102]:
    a = str(i)
    plt.loglog(kperp[0:nx],spettro_B_perp[i]*(kperp[0:nx])**(alpha),'xkcd:green',linewidth=wid,label=r'$\mathcal{P}_{B}^{\perp} k^{5/3}_{\perp}$')
    plt.loglog(kperp[0:nx],spettro_B_perp[i]*(kperp[0:nx])**(beta),'xkcd:pink',linewidth=wid,label=r'$\mathcal{P}_{B}^\perp k^{3/2}_{\perp}$')
    plt.loglog(kparal[0:nx],spettro_B_parall[i]*(kparal[0:nx])**(gamma),'xkcd:orange',linewidth=wid,label=r'$\mathcal{P}_{B}^\parallel k^{3.5}_{\parallel}$')    
    #plt.loglog(kperp[0:nx],spettro_B_iso[i]*(kperp[0:nx])**(alpha),label='Spettro isotropico')
    plt.title('Compensated magnetic spectra',fontsize='x-large')
    plt.legend(fontsize='x-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%%
nx=len(spettro_B_perp[0]) ; l=3 ; punt=2*l+1 ; temp=102
slope_bperp=[] ; slope_bpara=[]
for i in range(l,nx-l):
    slope_bperp.append(stats.linregress(np.log(kperp[i-l:i+l]),np.log(spettro_B_perp[temp][i-l:i+l]))[0])
    slope_bpara.append(stats.linregress(np.log(kparal[i-l:i+l]),np.log(spettro_B_parall[temp][i-l:i+l]))[0])

slope_bperp = np.asarray(slope_bperp).reshape(-1)
slope_bpara = np.asarray(slope_bpara).reshape(-1)

kolmogorov=np.zeros(nx,dtype=float) ; kraiknan=np.zeros(nx,dtype=float) ; parall=np.zeros(nx,dtype=float) ; strano=np.zeros(nx,dtype=float) ; stranissimo=np.zeros(nx,dtype=float) 
for i in range(nx):
    kolmogorov[i]=-5/3
    kraiknan[i]=-3/2
    parall[i]=-2
    strano[i]=-3.5
    stranissimo[i]=-4.5

plt.figure(237)
title='Spectral Index B. fit a %i punti' %punt
plt.title('Spectral index (B) (%i points)'%punt,fontsize='x-large')
plt.semilogx(kperp[0:nx],kraiknan,'xkcd:pink',linewidth=wid,label='-3/2')
plt.semilogx(kperp[0:nx],kolmogorov,'xkcd:green',linewidth=wid,label='-5/3')
#plt.semilogx(kperp[0:nx],parall,'xkcd:purple',linewidth=wid,label='-2')
plt.semilogx(kperp[0:nx],strano,'xkcd:orange',linewidth=wid,label='-3.5')
#plt.semilogx(kperp[0:nx],stranissimo,'xkcd:teal',linewidth=wid,label='-4.5')
plt.semilogx(kperp[l:nx-l],slope_bperp,'xkcd:red',linewidth=wid,label=r'$\alpha_{\perp}$')
plt.semilogx(kparal[l:nx-l],slope_bpara,'xkcd:blue',linewidth=wid,label=r'$\alpha_{\parallel}$')
plt.xlim(l,len(slope_bperp)) ; plt.ylim(-6,0)#;  plt.ylim(-11,1)
plt.legend(fontsize='x-large', loc='upper right')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')


#%% B 2D
alpha=5/3
beta=3/2

nx=len(spettro_B_rad[0]) 
kolmogorov=np.zeros(nx,dtype=float)
for i in range(nx):
    kolmogorov[i]=2*10**14

#m = np.argmax(rms_Jtot) ; l=11
#media_spettri_B = sum(spettro_B_rad[m-l:m+l+1])/(2*l+1)

plt.figure(123)
plt.loglog(k[0:nx],kolmogorov,'xkcd:black',linewidth=wid)
title = 'Spettro_radiale_B,alpha=%1.2f,beta=%1.2f'%(alpha,beta)
plt.title('Compensated magnetic spectra', fontsize='x-large')
plt.vlines(15, 10**9, 10**16,linestyles='dotted',colors='xkcd:purple',linewidth=wid)
plt.vlines(nx -1 -15, 10**9, 10**16,linestyles='dotted',colors='xkcd:purple',linewidth=wid)
plt.ylim((3*10**11,9*10**12)) #A2=C2
for i in [85]:
    a = str(i)
    plt.loglog(k[0:nx],spettro_B_rad[i]*(k[0:nx])**(alpha),'xkcd:green',linewidth=wid,label=r'$\mathcal{P}_B \, k^{5/3}$')
    plt.loglog(k[0:nx],spettro_B_rad[i]*(k[0:nx])**(beta),'xkcd:pink',linewidth=wid,label=r'$\mathcal{P}_B \, k^{3/2}$')
#plt.loglog(k[0:nx],media_spettri_B*(k[0:nx])**(alpha),'xkcd:black',label='Sp medio')   
plt.legend(loc='upper left',fontsize='x-large')

#if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')
if (save==True): plt.savefig(directory+'immagini/'+title+'max_turb'+'.png', format='png',dpi=300,bbox_inches='tight')   
#if (save==True): plt.savefig(directory+'immagini/'+title+'_media.png', format='png',dpi=300,bbox_inches='tight')

#%%
nx=len(spettro_B_rad[0]) ; l=20 ; punt=2*l+1 ; temp=85
slope_brad=[]
for i in range(l,nx-l):
    slope_brad.append(stats.linregress(np.log(k[i-l:i+l]),np.log(spettro_B_rad[temp][i-l:i+l]))[0])

slope_brad = np.asarray(slope_brad).reshape(-1)

kolmogorov=np.zeros(nx,dtype=float) ; kraiknan=np.zeros(nx,dtype=float) 
for i in range(nx):
    kolmogorov[i]=-5/3
    kraiknan[i]=-3/2

plt.figure(237)
title='Spectral Index B. fit a %i punti'%punt
plt.title('Spectral index (B) (%i points)'%punt,fontsize='x-large')
plt.semilogx(k[0:nx],kolmogorov,'xkcd:green',linewidth=wid,label='-5/3')
plt.semilogx(k[0:nx],kraiknan,'xkcd:pink',linewidth=wid,label='-3/2')
plt.semilogx(k[l:nx-l],slope_brad,'xkcd:red',linewidth=wid,label=r'$\alpha$')
plt.xlim(l,len(slope_brad)) ; plt.ylim (-3,-0) #; plt.ylim(-11,1)
plt.legend(fontsize='x-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% spettro_V_gyro
nx=len(spettro_V_gyro[0][0,:])
vmin=0 ; vmax=14
for i in [102]:
    a = str(i)
    plt.figure(i)
    title='Spettro_V_gyrotropico_max_turb(%s)' %a
    plt.contourf(spettro_V_gyro[i],cmap=cm.jet,levels=np.logspace(vmin, vmax, 255),locator=mplt.ticker.LogLocator())
    mplt.pyplot.yscale('log') ; mplt.pyplot.xscale('log')
    plt.ylim((1,nx-1)) ; plt.xlim((1,nx-1))
    plt.ylabel('k parallel') ; plt.xlabel('k perpendicular') 
    plt.colorbar(ticks=np.logspace(vmin,vmax,vmax-vmin+1))
    plt.title(title)
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')


#%% V ridotto
alpha=5/3
beta=3/2
gamma=7/2
#alpha=0 ; beta=0 ; gamma=0
nx=len(spettro_V_perp[0]) 
kolmogorov=np.zeros(nx,dtype=float)
for i in range(nx):
    kolmogorov[i]=1.5*10**14
    
plt.figure(i)
plt.vlines(3, 10**9, 10**16,linestyles='dotted',colors='xkcd:purple',linewidth=wid)
plt.vlines(nx -1 -3, 10**9, 10**16,linestyles='dotted',colors='xkcd:purple',linewidth=wid)
#plt.loglog(kperp[0:nx],kolmogorov,'xkcd:black',linewidth=wid)
#plt.loglog(kperp[0:nx],kolmogorov*10.5,'xkcd:black',linewidth=wid)
title='Spettri_ridotti_V_max_turb,alpha=%1.2f,beta=%1.2f,gamma=%1.2f'%(alpha,beta,gamma)
plt.ylim((2*10**13,2*10**15))
for i in [102]:
    a = str(i)
    plt.loglog(kperp[0:nx],spettro_V_perp[i]*(kperp[0:nx])**(alpha),'xkcd:green',linewidth=wid,label=r'$\mathcal{P}_{v}^{\perp} k^{5/3}_{\perp}$')
    plt.loglog(kperp[0:nx],spettro_V_perp[i]*(kperp[0:nx])**(beta),'xkcd:pink',linewidth=wid,label=r'$\mathcal{P}_{v}^{\perp} k^{3/2}_{\perp}$')
    plt.loglog(kparal[0:nx],spettro_V_parall[i]*(kparal[0:nx])**(gamma),'xkcd:orange',linewidth=wid,label=r'$\mathcal{P}_{v}^{\parallel} k^{3.5}_{\parallel}$')
    #plt.loglog(kperp[0:nx],spettro_V_iso[i]*(kperp[0:nx])**(alpha),label='Spettro isotropico %s'%a)
    plt.title('Compensated velocity spectra',fontsize='x-large')
    plt.legend(fontsize='x-large',bbox_to_anchor=(0.31,0.44))
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%%
nx=len(spettro_V_perp[0]) ; l=3 ; punt=2*l+1 ; temp=102
slope_vperp=[] ; slope_vpara=[]
for i in range(l,nx-l):
    slope_vperp.append(stats.linregress(np.log(kperp[i-l:i+l]),np.log(spettro_V_perp[temp][i-l:i+l]))[0])
    slope_vpara.append(stats.linregress(np.log(kparal[i-l:i+l]),np.log(spettro_V_parall[temp][i-l:i+l]))[0])

slope_vperp = np.asarray(slope_vperp).reshape(-1)
slope_vpara = np.asarray(slope_vpara).reshape(-1)

kolmogorov=np.zeros(nx,dtype=float) ; kraiknan=np.zeros(nx,dtype=float) ; parall=np.zeros(nx,dtype=float) ; strano=np.zeros(nx,dtype=float) ; stranissimo=np.zeros(nx,dtype=float) 
for i in range(nx):
    kolmogorov[i]=-5/3
    kraiknan[i]=-3/2
    parall[i]=-2
    strano[i]=-3.5
    stranissimo[i]=-4.5

plt.figure(237)
title='Spectral Index V. fit a %i punti'%punt
plt.title('Spectral index (V) (%i points)'%punt,fontsize='x-large')
plt.semilogx(kperp[0:nx],kraiknan,'xkcd:pink',linewidth=wid,label='-3/2')
plt.semilogx(kperp[0:nx],kolmogorov,'xkcd:green',linewidth=wid,label='-5/3')
#plt.semilogx(kperp[0:nx],parall,'xkcd:purple',linewidth=wid,label='-2')
plt.semilogx(kperp[0:nx],strano,'xkcd:orange',linewidth=wid,label='-3.5')
#plt.semilogx(kperp[0:nx],stranissimo,'xkcd:teal',linewidth=wid,label='-4.5')
plt.semilogx(kperp[l:nx-l],slope_vperp,'xkcd:red',linewidth=wid,label=r'$\alpha_{\perp}$')
plt.semilogx(kparal[l:nx-l],slope_vpara,'xkcd:blue',linewidth=wid,label=r'$\alpha_{\parallel}$')
plt.xlim(l,len(slope_vperp)) ;  plt.ylim(-6,0)
plt.legend(fontsize='x-large', loc='upper right')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% V 2D
alpha=5/3
beta=3/2

nx=len(spettro_B_rad[0]) 
kolmogorov=np.zeros(nx,dtype=float)
for i in range(nx):
    kolmogorov[i]=1*10**13
    
#m = np.argmax(rms_Jtot) ; l=11
#media_spettri_V = sum(spettro_V_rad[m-l:m+l+1])/(2*l+1)

plt.figure(246)
#plt.loglog(k[0:nx],kolmogorov,'xkcd:black',linewidth=wid)
title = 'Spettro_radiale_V,alpha=%1.2f,beta=%1.2f'%(alpha,beta)
plt.title('Compensated velocity spectra', fontsize='x-large')
plt.vlines(15, 10**9, 10**16,linestyles='dotted',colors='xkcd:purple',linewidth=wid)
plt.vlines(nx -1 -15, 10**9, 10**16,linestyles='dotted',colors='xkcd:purple',linewidth=wid)
#plt.ylim((6*10**10,1*10**12)) #A2
plt.ylim((6*10**9,2.5*10**11)) #C2
for i in [51]:
    a = str(i)
    #plt.figure(i)
    plt.loglog(k[0:nx],spettro_V_rad[i]*(k[0:nx])**(alpha),'xkcd:green',linewidth=wid,label=r'$\mathcal{P}_v \, k^{5/3}$')
    plt.loglog(k[0:nx],spettro_V_rad[i]*(k[0:nx])**(beta),'xkcd:pink',linewidth=wid,label=r'$\mathcal{P}_v \, k^{3/2}$')
#plt.loglog(k[0:nx],media_spettri_V*(k[0:nx])**(alpha),'xkcd:black',label='Sp medio') 
plt.legend(loc='upper left',fontsize='x-large')

#if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')
if (save==True): plt.savefig(directory+'immagini/'+title+'max_turb'+'.png', format='png',dpi=300,bbox_inches='tight')
#if (save==True): plt.savefig(directory+'immagini/'+title+'_media.png', format='png',dpi=300,bbox_inches='tight')
    
#%%
nx=len(spettro_V_rad[0]) ; l=20 ; punt=2*l+1 ; temp=85
slope_vrad=[]
for i in range(l,nx-l):
    slope_vrad.append(stats.linregress(np.log(k[i-l:i+l]),np.log(spettro_V_rad[temp][i-l:i+l]))[0])

slope_vrad = np.asarray(slope_vrad).reshape(-1)

kolmogorov=np.zeros(nx,dtype=float) ; kraiknan=np.zeros(nx,dtype=float) 
for i in range(nx):
    kolmogorov[i]=-5/3
    kraiknan[i]=-3/2

plt.figure(237)
title='Spectral Index V. fit a %i punti'%punt
plt.title('Spectral index (V) (%i points)'%punt,fontsize='x-large')
plt.semilogx(k[0:nx],kolmogorov,'xkcd:green',linewidth=wid,label='-5/3')
plt.semilogx(k[0:nx],kraiknan,'xkcd:pink',linewidth=wid,label='-3/2')
plt.semilogx(k[l:nx-l],slope_vrad,'xkcd:red',linewidth=wid,label=r'$\alpha$')
plt.xlim(l,len(slope_vrad)) ;  plt.ylim(-3,0)
plt.legend(fontsize='x-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')


#%% cross helicity
title='cross helicity (giusta) new label'
plt.figure(2001) 
plt.plot(t_out,np.zeros((len(t_out))),'xkcd:white',label=nome)
plt.plot(t_out,np.zeros(len(chw)),'xkcd:purple',linewidth=wid,linestyle='dashed',label=r'$zero$ line')
plt.plot(t_out,chw,'xkcd:blue',linewidth=wid,label=r'$\mathcal{C\!H}_R$') 
plt.plot(t_out,ch,'xkcd:red',linewidth=wid,label=r'$\mathcal{C\!H}$')
plt.legend(fontsize='x-large', bbox_to_anchor=(1,0.58)) ; plt.title('Cross helicity',fontsize='xx-large')    #bbox_to_anchor=(1,0.65) C2
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

title='cross helicity (giusta) new'
plt.figure(2002) 
plt.plot(t_out,np.zeros((len(t_out))),'xkcd:white',label=nome)
plt.plot(t_out,np.zeros(len(chw)),'xkcd:purple',linewidth=wid,linestyle='dashed')
plt.plot(t_out,chw,'xkcd:blue',linewidth=wid)
plt.plot(t_out,ch,'xkcd:red',linewidth=wid)
plt.legend(fontsize='x-large', bbox_to_anchor=(1,0.97)) ; plt.title('Cross helicity',fontsize='xx-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%%

title='magnetic helicity new label'
plt.figure(2003) 
plt.plot(t_out,np.zeros((len(t_out))),'xkcd:white',label=nome)
plt.plot(t_out,np.zeros(len(helicityB)),'xkcd:purple',linewidth=wid,linestyle='dashed',label=r'$zero$ line')
plt.plot(t_out,helicityB_g*((2*np.pi)/NNN)**dd,'xkcd:blue',linewidth=wid,label=r'helicity $B$') 
plt.legend(fontsize='x-large') ; plt.title('Magnetic helicity',fontsize='xx-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

title='magnetic helicity new'
plt.figure(2004) 
plt.plot(t_out,np.zeros((len(t_out))),'xkcd:white',label=nome)
plt.plot(t_out,np.zeros(len(helicityB)),'xkcd:purple',linewidth=wid,linestyle='dashed')
plt.plot(t_out,helicityB_g*((2*np.pi)/NNN)**dd,'xkcd:blue',linewidth=wid) 
plt.legend(fontsize='x-large') ; plt.title('Magnetic helicity',fontsize='xx-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

title='current helicity new label'
plt.figure(2005) 
plt.plot(t_out,np.zeros((len(t_out))),'xkcd:white',label=nome)
plt.plot(t_out,np.zeros(len(helicityJ)),'xkcd:purple',linewidth=wid,linestyle='dashed',label=r'$zero$ line')
plt.plot(t_out,helicityJ_g*((2*np.pi)/NNN)**dd,'xkcd:red',linewidth=wid,label=r'helicity $J$') 
plt.legend(fontsize='x-large') ; plt.title('Current helicity',fontsize='xx-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

title='current helicity new'
plt.figure(2006) 
plt.plot(t_out,np.zeros((len(t_out))),'xkcd:white',label=nome)
plt.plot(t_out,np.zeros(len(helicityJ)),'xkcd:purple',linewidth=wid,linestyle='dashed')
plt.plot(t_out,helicityJ_g*((2*np.pi)/NNN)**dd,'xkcd:red',linewidth=wid) 
plt.legend(loc='upper right',fontsize='x-large') ; plt.title('Current helicity',fontsize='xx-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% rms J e vor
title='rms J new'
if post==True : title= 'rms J post'
plt.figure(200) 
tmax=t_out[np.argmax(np.sqrt(rms_Jtot))]
plt.plot(t_out,np.zeros((len(t_out))),'xkcd:white',label=nome)
plt.vlines(tmax, 0, np.amax(np.sqrt(rms_Jtot)),linestyles='dashed',label='t=%1.2f'%tmax,colors='xkcd:purple',linewidth=wid)
plt.plot(t_out,np.sqrt(rms_Jx),'xkcd:green',linewidth=wid,label=r'$\sigma_N$(J$x$)') 
plt.plot(t_out,np.sqrt(rms_Jy),'xkcd:red',linewidth=wid,label=r'$\sigma_N$(J$y$)')
plt.plot(t_out,np.sqrt(rms_Jz),'xkcd:blue',linewidth=wid,label=r'$\sigma_N$(J$z$)') 
plt.plot(t_out,np.sqrt(rms_Jtot),'xkcd:black',linewidth=wid,label=r'$\sigma_N$(J)') 
plt.legend(fontsize='x-large', loc='upper left') ; plt.title('RMS density current',fontsize='xx-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

title='rms J new no lab'
if post==True : title= 'rms J post'
plt.figure(201) 
tmax=t_out[np.argmax(np.sqrt(rms_Jtot))]
plt.plot(t_out,np.zeros((len(t_out))),'xkcd:white',label=nome)
plt.vlines(tmax, 0, np.amax(np.sqrt(rms_Jtot)),linestyles='dashed',label='t=%1.2f'%tmax,colors='xkcd:purple',linewidth=wid)
plt.plot(t_out,np.sqrt(rms_Jx),'xkcd:green',linewidth=wid) 
plt.plot(t_out,np.sqrt(rms_Jy),'xkcd:red',linewidth=wid)
plt.plot(t_out,np.sqrt(rms_Jz),'xkcd:blue',linewidth=wid) 
plt.plot(t_out,np.sqrt(rms_Jtot),'xkcd:black',linewidth=wid) 
plt.legend(fontsize='x-large', loc='upper left') ; plt.title('RMS density current',fontsize='xx-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

titlesave='rms vor new'
if post==True : titlesave= 'rms vor post' ; title='rms $\omega$ post'
plt.figure(250) 
tmax=t_out[np.argmax(np.sqrt(rms_vortot))]
plt.plot(t_out,np.zeros((len(t_out))),'xkcd:white',label=nome)
plt.vlines(tmax, 0, np.amax(np.sqrt(rms_vortot)),linestyles='dashed',label='t=%1.2f'%tmax,colors='xkcd:purple',linewidth=wid)
plt.plot(t_out,np.sqrt(rms_vorx),'xkcd:green',linewidth=wid,label=r'$\sigma_N$($\omega x$)') 
plt.plot(t_out,np.sqrt(rms_vory),'xkcd:red',linewidth=wid,label=r'$\sigma_N$($\omega y$)')
plt.plot(t_out,np.sqrt(rms_vorz),'xkcd:blue',linewidth=wid,label=r'$\sigma_N$($\omega z$)') 
plt.plot(t_out,np.sqrt(rms_vortot),'xkcd:black',linewidth=wid,label=r'$\sigma_N$($\omega$)') 
plt.legend(fontsize='x-large', loc='upper left') ; plt.title('RMS vorticity',fontsize='xx-large')
if (save==True): plt.savefig(directory+'immagini/'+titlesave+'.png', format='png',dpi=300,bbox_inches='tight')

titlesave='rms vor new no lab'
if post==True : titlesave= 'rms vor post' ; title='rms $\omega$ post'
plt.figure(251) 
tmax=t_out[np.argmax(np.sqrt(rms_vortot))]
plt.plot(t_out,np.zeros((len(t_out))),'xkcd:white',label=nome)
plt.vlines(tmax, 0, np.amax(np.sqrt(rms_vortot)),linestyles='dashed',label='t=%1.2f'%tmax,colors='xkcd:purple',linewidth=wid)
plt.plot(t_out,np.sqrt(rms_vorx),'xkcd:green',linewidth=wid) 
plt.plot(t_out,np.sqrt(rms_vory),'xkcd:red',linewidth=wid)
plt.plot(t_out,np.sqrt(rms_vorz),'xkcd:blue',linewidth=wid) 
plt.plot(t_out,np.sqrt(rms_vortot),'xkcd:black',linewidth=wid) 
plt.legend(fontsize='x-large', loc='upper left') ; plt.title('RMS vorticity',fontsize='xx-large')
if (save==True): plt.savefig(directory+'immagini/'+titlesave+'.png', format='png',dpi=300,bbox_inches='tight')


#%% rms b
title='rms_b'
plt.figure(300) 
plt.plot(t_out,np.sqrt(rms_bx),'xkcd:green',linewidth=wid,label=r'$\sigma_N(Bx)$') 
plt.plot(t_out,np.sqrt(rms_by),'xkcd:red',linewidth=wid,label=r'$\sigma_N(By)$')
plt.plot(t_out,np.sqrt(rms_bz),'xkcd:blue',linewidth=wid,label=r'$\sigma_N(Bz)$') 
plt.plot(t_out,np.sqrt(rms_btot),'xkcd:black',linewidth=wid,label=r'$\sigma_N(B)$') 
plt.legend(loc='upper right',fontsize='x-large') ; plt.title('RMS magnetic field',fontsize='x-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% rms v
title='rms_v'
plt.figure(400) 
plt.plot(t_out,np.sqrt(rms_vx),'xkcd:green',linewidth=wid,label=r'$\sigma_N(vx)$') 
plt.plot(t_out,np.sqrt(rms_vy),'xkcd:red',linewidth=wid,label=r'$\sigma_N(vy)$')
plt.plot(t_out,np.sqrt(rms_vz),'xkcd:blue',linewidth=wid,label=r'$\sigma_N(vz)$') 
plt.plot(t_out,np.sqrt(rms_vtot),'xkcd:purple',linewidth=wid,label=r'$\sigma_N(v)$') 
plt.legend(loc='upper right',fontsize='x-large') ; plt.title('RMS velocity field',fontsize='x-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')
 
#%% contronto rms b e v tot
title='confronto rms b v'
plt.figure(472,[9,3]) 
plt.plot(t_out,np.sqrt(rms_vtot)+0.15,'xkcd:purple',linewidth=wid,label=r'$\sigma_N(v)$ [translated]') 
plt.plot(t_out,np.sqrt(rms_btot),'xkcd:black',linewidth=wid,label=r'$\sigma_N(B)$') 
plt.legend(loc='upper right',fontsize='x-large') ; plt.title('Comparison total RMS of magnetic and velocity fields',fontsize='x-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% rms rh e pe
title='rms_rh_pe'
plt.figure(404) 
plt.plot(t_out,np.sqrt(rms_rh),'xkcd:red',linewidth=wid,label='rms rh') 
plt.plot(t_out,np.sqrt(rms_pe),'xkcd:blue',linewidth=wid,label='rms pe')
plt.legend() ; plt.title(title)
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

##%% rms cross helicity
#title='Cross Helicity'
#plt.figure(404) 
#plt.errorbar(t_out,ava_ch,yerr=rms_ch,capsize=1.7,color='xkcd:blue',linewidth=wid,label='CH ava')
#plt.errorbar(t_out,ava_chro,yerr=rms_chro,capsize=1.7,color='xkcd:red',linewidth=wid,label='CH(rho) ava') 
#plt.plot(t_out,np.zeros(len(ava_ch)),'xkcd:purple',linewidth=wid,linestyle='dashed',label='zero line')
#plt.title(title) ; plt.legend()
#if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%%
##%% rms Va
#title='rms_Va'
#plt.figure(451) 
#plt.plot(t_out,np.sqrt(rms_Vax),'xkcd:green',linewidth=wid,label='rms Vax') 
#plt.plot(t_out,np.sqrt(rms_Vay),'xkcd:red',linewidth=wid,label='rms Vay')
#plt.plot(t_out,np.sqrt(rms_Vaz),'xkcd:blue',linewidth=wid,label='rms Vaz') 
#plt.plot(t_out,np.sqrt(rms_Vatot),'xkcd:black',linewidth=wid,label='rms Vatot') 
#plt.legend() ; plt.title(title)
#if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')
#
##%% rms zp zm
#title='rms_zp_zm^2'
#plt.figure(496) 
#plt.plot(t_out,np.sqrt(rms_zpx),'xkcd:azure',linewidth=wid,label='rms zpx') 
#plt.plot(t_out,np.sqrt(rms_zpy),'xkcd:lightgreen',linewidth=wid,label='rms zpy')
#plt.plot(t_out,np.sqrt(rms_zpz),'xkcd:teal',linewidth=wid,label='rms zpz') 
#plt.plot(t_out,np.sqrt(rms_zptot),'xkcd:blue',linewidth=wid,label='rms zptot') 
#plt.plot(t_out,np.sqrt(rms_mod_zp),'xkcd:violet',linewidth=wid,label='rms mod zp') 
##plt.plot(t_out,np.sqrt(ava_zp),label='avarage zp')
#plt.plot(t_out,np.sqrt(rms_zmx),'xkcd:salmon',linewidth=wid,label='rms zmx') 
#plt.plot(t_out,np.sqrt(rms_zmy),'xkcd:tan',linewidth=wid,label='rms zmy')
#plt.plot(t_out,np.sqrt(rms_zmz),'xkcd:orange',linewidth=wid,label='rms zmz') 
#plt.plot(t_out,np.sqrt(rms_zmtot),'xkcd:red',linewidth=wid,label='rms zmtot')
#plt.plot(t_out,np.sqrt(rms_mod_zm),'xkcd:magenta',linewidth=wid,label='rms mod zm') 
##plt.plot(t_out,np.sqrt(ava_zm),label='avarage zm')
#plt.legend() ; plt.title(title)
#if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

##%% rms cross helicity
#title='rms_CrossHelicity'
#plt.figure(404) 
#plt.plot(t_out,Crohel_ava,'xkcd:blue',linewidth=wid,label='Crohel ava') 
#plt.plot(t_out,Crohel_rms,'xkcd:red',linewidth=wid,label='Crohel rms') 
#plt.plot(t_out,np.zeros(len(Crohel_ava)),'xkcd:green',linewidth=wid,label='zero line')
#plt.title(title) ; plt.legend()
#if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% massimi J
title='massimi_J'
plt.figure(100) 
plt.plot(t_out,np.sqrt(max_Jx),'xkcd:green',linewidth=wid,label='max Jx') 
plt.plot(t_out,np.sqrt(max_Jy),'xkcd:red',linewidth=wid,label='max Jy') 
plt.plot(t_out,np.sqrt(max_Jz),'xkcd:blue',linewidth=wid,label='max Jz') 
plt.legend() ; plt.title(title)
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% massimi b
title='massimi_B'
plt.figure(500) 
plt.plot(t_out,np.sqrt(max_bx),'xkcd:green',linewidth=wid,label='max bx') 
plt.plot(t_out,np.sqrt(max_by),'xkcd:red',linewidth=wid,label='max by') 
plt.plot(t_out,np.sqrt(max_bz),'xkcd:blue',linewidth=wid,label='max bz') 
plt.legend() ; plt.title(title)
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% massimi V
title='massimi_V new'
plt.figure(600) 
plt.plot(t_out,max_vx,'xkcd:green',linewidth=wid,label=r'max $vx$') 
plt.plot(t_out,max_vy,'xkcd:red',linewidth=wid,label=r'max $vy$') 
plt.plot(t_out,max_vz,'xkcd:blue',linewidth=wid,label=r'max $vz$') 
#plt.plot(t_out,np.sqrt(max_vz*max_vz+max_vy*max_vy+max_vx*max_vx),'xkcd:black',linewidth=wid,label=r'max $v$') #questa lasciala fare che è meglio
plt.legend(fontsize='x-large') ; plt.title('Maximum velocity', fontsize='xx-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% massimi e minimi pe e rh
title='max_e_min_pe_rh new'
plt.figure(700) 
plt.plot(t_out,max_pe,'xkcd:purple',linewidth=wid,label=r'max $P$') 
plt.plot(t_out,min_pe,'xkcd:pink',linewidth=wid,label=r'min $P$') 
plt.plot(t_out,max_rh,'xkcd:green',linewidth=wid,label=r'max $\rho$') 
plt.plot(t_out,min_rh,'xkcd:bright green',linewidth=wid,label=r'min $\rho$') 
plt.legend(fontsize='x-large', bbox_to_anchor=(0.31,0.40)) ; plt.title('Pressure and density', fontsize='xx-large')
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% Etot  2D
alpha=3/2
nx=len(spettro_Etot_rad[0]) 
kolmogorov=np.zeros(nx,dtype=float)
for i in range(nx):
    kolmogorov[i]=1.2*10**14

#m = np.argmax(rms_Jtot) ; l=11
#media_spettri_B = sum(spettro_B_rad[m-l:m+l+1])/(2*l+1)

plt.figure(285)
plt.loglog(k[0:nx],kolmogorov,'xkcd:black',linewidth=wid)
title = 'Spettro_radiale_Etot,alpha=%1.2f'%alpha
plt.title(title)
plt.ylim((1*10**12,5*10**14))
for i in [20]:
    a = str(i)
    plt.loglog(k[0:nx],spettro_Etot_rad[i]*(k[0:nx])**(alpha),linewidth=wid,label='Sp. %s'%a)
#plt.loglog(k[0:nx],media_spettri_B*(k[0:nx])**(alpha),'xkcd:black',label='Sp medio')   
plt.legend()

#if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')
if (save==True): plt.savefig(directory+'immagini/'+title+'max_turb'+'.png', format='png',dpi=300,bbox_inches='tight')   
#if (save==True): plt.savefig(directory+'immagini/'+title+'_media.png', format='png',dpi=300,bbox_inches='tight')

 #%% spettro_Etot_gyro
nx=len(spettro_Etot_gyro[0][0,:])
vmin=0 ; vmax=14
for i in [102]:
    a = str(i)
    #plt.figure(i,figsize=(12.8,9.6))
    plt.figure(i)
    title='Spettro_Etot_gyrotropico_max_turb(%s)' %a
    plt.contourf(spettro_Etot_gyro[i],cmap=cm.jet,levels=np.logspace(vmin,vmax, 255),locator=mplt.ticker.LogLocator())
    mplt.pyplot.yscale('log') ; mplt.pyplot.xscale('log') 
    plt.ylim((1,nx-1)) ; plt.xlim((1,nx-1)) 
    plt.ylabel('k parallel') ; plt.xlabel('k perpendicular') 
    plt.colorbar(ticks=np.logspace(vmin,vmax,vmax-vmin+1)) 
    plt.title(title) 
    if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%% K_para e K_perp
start=time.time()
npar = len(K_para[0])
nper = len(K_perp[0])
confronto=np.zeros((npar,nper),dtype=float)
for t in range(tempi+1):
    for j in range(npar):
        for i in range(nper):
            confronto[j,i]=np.abs(K_perp[t][i]-K_para[t][j])

relaz=[]
for t in range(tempi+1):
    relaz.append(np.zeros(npar,dtype=float))

for t in range(tempi+1):
    for i in range(npar):
        relaz[t][i]=np.argmin(confronto[:,i])

#m = np.argmax(rms_Jtot) ; l=10
#media_bal_K_paraB = sum(K_paraB[m-l:m+l+1])/(2*l+1)
#media_bal_K_perpB = sum(K_perpB[m-l:m+l+1])/(2*l+1)
        
end=time.time()
print('time : %i'%(end-start),'sec')

#%%
alpha=2/3
vmin=0 ; vmax=14
plt.figure(941)   
title='Balance Kpar e Kperp alpha=%1.2f'%alpha
for i in [102]:
    a = str(i)
    plt.contourf(spettro_Etot_gyro[i],cmap=cm.jet,levels=np.logspace(vmin,vmax, 29),locator=mplt.ticker.LogLocator())
    plt.loglog(kperp[0:nper],relaz[i],'xkcd:black',linewidth=1.,label='Balance %s'%a)
    plt.loglog(kperp[1:nper],(kperp[1:nper]*(2))**alpha,'xkcd:purple',linewidth=1.2,linestyle='dashed',label='Critical Balance')
#plt.loglog(media_bal_K_perpB,media_bal_K_paraB*media_bal_K_perpB**alpha,'xkcd:black',label='Avarage Balance')
plt.ylim((0.7,npar)) ; plt.xlim((0.7,nper))
plt.title(title)
plt.colorbar(ticks=np.logspace(vmin,vmax,vmax-vmin+1))
plt.ylabel('k parallel') ; plt.xlabel('k perpendicular') 

#piatto=np.zeros(len(K_perp[0]),dtype=float)
#piatto[:]=np.amin(media_bal_K_paraB*media_bal_K_perpB**alpha)
#plt.loglog(media_bal_K_perpB,piatto)

if (save==True): plt.savefig(directory+'immagini/'+title+'1.png', format='png',dpi=300,bbox_inches='tight')
        
#%% Etot ridotto
alpha=3/2
beta=4.5
nx=len(spettro_Etot_perp[0]) 
kolmogorov=np.zeros(nx,dtype=float)
for i in range(nx):
    kolmogorov[i]=9.5*10**14

plt.figure(490)
plt.loglog(kperp[0:nx],kolmogorov,'xkcd:black',linewidth=wid)
plt.loglog(kperp[0:nx],kolmogorov*6.3,'xkcd:black',linewidth=wid)
title='Spettri_ridotti_Etot_max_turb,alpha=%1.2f,beta=%1.2f'%(alpha,beta)
plt.ylim((9*10**13,7*10**15))
for i in [102]:
    a = str(i)
    plt.loglog(kperp[0:nx],spettro_B_perp[i]*(kperp[0:nx])**(alpha),linewidth=wid,label='Sp perp %s (alpha)'%a)
    plt.loglog(kparal[0:nx],spettro_B_parall[i]*(kparal[0:nx])**(beta)/15,linewidth=wid,label='Sp par %s (beta)'%a)
    #plt.loglog(kperp[0:nx],spettro_B_iso[i]*(kperp[0:nx])**(alpha),label='Spettro isotropico')
    plt.legend()
    plt.title(title)
if (save==True): plt.savefig(directory+'immagini/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

