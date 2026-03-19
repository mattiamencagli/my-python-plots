# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 10:08:39 2019

@author: Mattia
"""

import numpy as np
import pylab 
import codecs

directory='/Users/Mattia/Desktop/Risultati L.Franci/256^3(init_onda_2D. LF. 256tempi)/'

with codecs.open(directory+'rms.dat',encoding='utf-8-sig') as f:
    rms=np.loadtxt(f)

m=181
j = np.zeros(m,dtype=int) ; 
for i in range(0,m):   
    j[i] = np.array([i])

time = rms[0,0:m]

Med_rho = rms[1,0:m]
Med_vx=rms[2,0:m] ; Med_vy=rms[3,0:m] ; Med_vz=rms[4,0:m] ; Med_V=Med_vy+Med_vz+Med_vx
Med_pe = rms[5,0:m]
Med_bx=rms[6,0:m] ; Med_by=rms[7,0:m] ; Med_bz=rms[8,0:m] ; Med_B=Med_by+Med_bz+Med_bx
Med_T=rms[9,0:m] ; Med_Etot=rms[10,0:m] ; Med_Ekin=rms[11,0:m]

Var_rho = rms[12,0:m]
Var_vx=rms[13,0:m] ; Var_vy=rms[14,0:m] ; Var_vz=rms[15,0:m] ; Var_V=Var_vy+Var_vz+Var_vx
Var_pe = rms[16,0:m]
Var_bx=rms[17,0:m] ; Var_by=rms[18,0:m] ; Var_bz=rms[19,0:m] ; Var_B=Var_by+Var_bz+Var_bx
Var_T=rms[20,0:m] ; Var_Etot=rms[21,0:m] ; Var_Ekin=rms[22,0:m]

#%%
pylab.figure(0)
pylab.plot(j,time)
pylab.title('tempo')

#%%
pylab.figure(12)
pylab.plot(j,Var_rho)
pylab.title('varianza densità')

#%%
pylab.figure(13)
pylab.plot(j,Var_vx)
pylab.title('varianza vx')
pylab.figure(14)
pylab.plot(j,Var_vy)
pylab.title('varianza vy')
pylab.figure(15)
pylab.plot(j,Var_vz)
pylab.title('varianza vz')
pylab.figure(31)
pylab.plot(j,Var_V)
pylab.title('varianza V TOTALE')

#%%
pylab.figure(17)
pylab.plot(j,Var_bx)
pylab.title('varianza bx')
pylab.figure(18)
pylab.plot(j,Var_by)
pylab.title('varianza by')
pylab.figure(19)
pylab.plot(j,Var_bz)
pylab.title('varianza bz')
pylab.figure(30)
pylab.plot(j,Var_B)
pylab.title('varianza B TOTALE')

#%%
pylab.figure(20)
pylab.plot(j,Var_T)
pylab.title('varianza Temperatura')

pylab.figure(22)
pylab.plot(j,Var_Ekin)
pylab.title('varianza Energia cinetica')
