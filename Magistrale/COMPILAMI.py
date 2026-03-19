# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:13:02 2019

@author: Mattia
"""

print('VIA')

import numpy as np ; import pylab as plt 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm 

#import matplotlib as mpl
#mpl.use('Agg')
#
#x=np.linspace(0,2*np.pi,100)
#y=np.sin(x)
#
#plt.figure(668)
#
#plt.plot(x,y,label='hehe')
#plt.legend()
#
#plt.savefig('provaproca.png', format='png',dpi=70,bbox_inches='tight')


#x=np.linspace(0,100,101)
#
#a=np.concatenate((x[90:100],x[0:20]),axis=None)
#
#z=np.concatenate((np.linspace(429,511,511-429+1),np.linspace(0,40,41)),axis=None)


#%%
'''gamma lorentz grafico'''
#v=np.linspace(0,1,1000)
#gam=(1-v**2)**(-0.5)
#
#plt.figure(0)
#plt.ylim(0,10)
##plt.plot(v,gam,'xkcd:white',label=r'$E_3$')
#plt.plot(v,gam,'xkcd:white',label=r'$F_3$')
#plt.legend(fontsize='x-large')
#plt.savefig('/Users/Mattia/Desktop/Flabel.png', format='png',dpi=300,bbox_inches='tight')

#%% 
'''coordinate massimo marice'''

#tot=np.argmax(matrix)
#
#toty=int(tot/ny)
#totx=tot-toty*ny
#%%
'''friedrik'''
#
#modva=0
#cs=2
#s=modva/cs
#
#title='Friedrich diagram va.div.cs=%1.2f'%s
#
#theta=np.linspace(0,2*np.pi,3000)
#
#va=np.abs(modva*np.cos(theta))
#
#csf=np.sqrt(0.5*( cs**2 + modva**2 + np.sqrt( (cs**2 + modva**2)**2 -4*(cs*modva*np.cos(theta))**2 ) ))
#css=np.abs(np.sqrt(0.5*( cs**2 + modva**2 - np.sqrt( (cs**2 + modva**2)**2 -4*(cs*modva*np.cos(theta))**2 ) )))
#
#norm=max(np.amax(va),np.amax(csf),np.amax(css))
#
#plt.figure(9999)#,[6.4*2, 4.8*2])
#plt.polar(theta,css/norm,'xkcd:blue',label='css')
#plt.polar(theta,csf/norm,'xkcd:green',label='csf')
#plt.polar(theta,va/norm,'xkcd:red',label='va')
#plt.ylim(0,1.16)
#plt.savefig('/Users/Mattia/Documents/Tesi_Mag/Tesi/Friedrich/'+title+'.png', format='png',dpi=300,bbox_inches='tight')

#%%
'''grafico va '''

#B0=np.linspace(0,5,1000)
#beta=np.linspace(0,5,1000)
#
##va_cont=np.zeros((len(B0),len(beta)),dtype=float)
##for i in range(len(B0)):
##    for j in range(len(beta)):
##        va_cont[i,j]=np.sqrt(B0[i]**2/(1+(1 + 2*beta[j])*B0[i]**2))
#        
##plt.figure()  
##plt.contourf(va_cont,255,cmap=cm.coolwarm)    #levels=np.linspace(-1.0,1.0,255)  
##plt.colorbar()
#        
#B0,beta=np.meshgrid(B0,beta)
#va=np.sqrt(B0**2/(1+(1 + 2*beta)*B0**2))
#
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#surf =ax.plot_surface(B0, beta, va, cmap=cm.coolwarm,rstride=10, cstride=10, alpha=None,antialiased=False)
#ax.set_xlabel(r'$B_0$')
#ax.set_ylabel(r'$\beta$')
#ax.set_zlabel(r'$V_a$')
#ax.view_init(23, 125)
#fig.colorbar(surf, shrink=0.5, aspect=14)
#title='Va3D'
#plt.savefig('/Users/Mattia/Documents/Tesi_Mag/Tesi/Va3D/'+title+'.png', format='png',dpi=400)

