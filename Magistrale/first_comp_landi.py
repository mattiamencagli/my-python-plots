# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 13:50:34 2019

@author: Mattia
"""

import h5py ; import time 

start=time.time()

bx = [] ; by = [] ; bz = [] ; vx = [] ; vy = [] ; vz = [] ; pe = [] ; rh = []

tempi=1

directory='/Users/Mattia/Desktop/Risultati L.Franci/256^3(init_onda_2D. LF. 256tempi)(Lz=64Lx)/'

print()
d=0
hf = h5py.File(directory + 'out000.h5', 'r') #per visualizzare che dataset ci sono nel file H5
#for key in hf.keys():
#    d=d+1
#    print('dataset %s = '%d,key)
#hf.close()
#print()
#
#for i in range(tempi+1):
#    if i<10: c = str(i) ; a = '0'+c
#    else: a = str(i)
#    print(a)
#    h=h5py.File(directory + 'out0%s.h5' %a)
#    bx.append( h.get('bx')[...] )
#    by.append( h.get('by')[...] )
#    bz.append( h.get('bz')[...] )
#    vx.append( h.get('vx')[...] )
#    vy.append( h.get('vy')[...] )
#    vz.append( h.get('vz')[...] )
#    rh.append( h.get('rh')[...] )
#    pe.append( h.get('pe')[...] )

end=time.time()
t=end-start
print()

h = str(int(t/3600)) ; m = str(int(t/60)) ; s =str(int(t)) ; ms=str(int(t*1000))
if int(h)<10: h='0'+h  
else: h=h 
if int(m)<10: m='0'+m
else: m=m
if int(s)<10: s='0'+s  
else: h=h 
if int(ms)<10: ms='00'+ms
elif int(ms)<100: ms='0'+ms
elif int(ms)<1000: ms=ms
else : ms = ms[len(ms)-3]+ms[len(ms)-2]+ms[len(ms)-1]

print('        h  m  s  ms')
print('time = {}:{}:{}:{}'.format(h,m,s,ms))
