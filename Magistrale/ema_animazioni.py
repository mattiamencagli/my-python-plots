# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

from data_echo import *
import matplotlib.animation as animation
from subprocess import call       #to call shell commands
from matplotlib import rc

plt.rc('text',usetex=False)
plt.rc('font',family='serif')

#filename='test_0905_noct_nostaticfloor_reflex_orto'
#filename='test_0903_Athenafloor_u2v3D'
filename='test'

#If true, it bypasses the derived quantities
fast=True

#type of output file
hdf5=True
#hdf5=False

#image=True
image=False


if hdf5:
  ext='.h5'
else:
  ext='.dat'

#switch={'rh':    [[0,0],[0],True,[1e-7,1]],          #plot coord,var,log,[Qmin,Qmax]
#	'pg':    [[0,1],[4],True,[1e-9,1e-2]],
#        'se':    [[0,2],[0,4],True,[None,None]],
#	'T':     [[1,0],[0,4],True,[None,None]],
#        'glf-1': [[1,1],[1,2,3],True,[1e-3,49]],
#        'b2':    [[1,2],[6,7,8,9,10,11],True,[None,None]]}

#if image:
switch={'rh':    [(0),[0],True,[1e-7,1],(0,0)],          #plot coord,var,log,[Qmin,Qmax]
        'pg':    [(1),[4],True,[1e-9,1e-2],(0,1)],
        'se':    [(2),[0,4],True,[1e-6,1e8],(0,2)],
        'T':     [(3),[0,4],True,[1e-6,1e4],(1,0)],
        'glf-1': [(4),[1,2,3],True,[1e-3,4],(1,1)],
        'b2':    [(5),[6,7,8,9,10,11],True,[1e-7,1e-2],(1,2)]}
#else:
#    switch={'rh':    [(0,0),[0],True,[1e-7,1]],          #plot coord,var,log,[Qmin,Qmax]
#            'pg':    [(0,1),[4],True,[1e-9,1e-2]],
#            'se':    [(0,2),[0,4],True,[1e-6,1e8]],
#    	    'T':     [(1,0),[0,4],True,[1e-6,1e4]],
#            'glf-1': [(1,2),[1,2,3],True,[1e-3,4]],
#            'b2':    [(1,3),[6,7,8,9,10,11],True,[1e-7,1e-2]]}

#switch={'vx':    [0,[1],False,[1e-2,1]],          #plot coord,var,log,[Qmin,Qmax]
#	'vy':    [1,[2],False,[-1,-1]],
#        'vz':    [2,[3],True,[None,None]],
#	'bx':    [3,[9],False,[1e-8,None]],
#        'by':    [4,[10],False,[1e-8,None]],
#        'bz':    [5,[11],False,[1e-8,None]]}

#var=switch[key][1]


xmax=50

direc = (raw_input('Select directory: '))
if (direc==''):
   direc='./'
else:
   direc = direc+'/' 

nout = int(raw_input('Select number of outputs to process: '))

#fig, ax = plt.subplots(nrows=2, ncols=3)
if image:
    fig, ax = plt.subplots(nrows=1, ncols=6)
else:
    fig, ax = plt.subplots(nrows=2, ncols=3)
fig.set_size_inches(18.5, 10.5,forward=True)
fig.tight_layout

data00 = echo()
data00 = read_data(direc,'out000'+ext,data00,fast=False,hdf5=hdf5,verb=False)[0] 
if data00.case==8:
    case='disk'
else:
    case='tearing'

# animation function.  This is called sequentially
def animate(n):
    global Qmin,Qmax

    data = echo()
    inputfile = 'out%03i'%(n)+ext    
    print inputfile
    var=[]
    for i in switch.iterkeys():
	var.extend(switch[i][1])
    var=unique(var)
    print var
#    var=[0,1,2,3,4,9,10,11]
    data = read_data(direc,inputfile,data,fast,hdf5=hdf5,variables=var,verb=False)[0]      

    imm=[None] * 10
    for key in switch.iterkeys():
    	log=switch[key][2]
    	Qmin=switch[key][3][0]
    	Qmax=switch[key][3][1]

    	if key=='rh':
    		Q=abs(data.rh)
    	elif key=='pg':
    	    Q=abs(data.pg)
    	elif key=='T':
    	    Q=abs(data.pg/data.rh)
    	elif key=='se':
    	    Q=abs(data.pg/data.rh**data.gamma)
    	elif key=='glf-1':
    	    v2=(data.vx**2*data00.g_11+
    	       data.vy**2*data00.g_22+
    	       data.vz**2*data00.g_33+
    	       data.vx*data.vz*data00.g_13*2.)
    	    Q=(1.-v2)**(-0.5)-1
    	elif key=='b2':
    	    Q=(data.bx**2*data00.g_11+
    	       data.by**2*data00.g_22+
    	       data.bz**2*data00.g_33+
    	       data.bx*data.bz*data00.g_13*2.-
    	       data.ex**2*data00.g_11-
    	       data.ey**2*data00.g_22-
    	       data.ez**2*data00.g_33-
    	       data.ex*data.ez*data00.g_13*2.)
    	elif key=='vx':
    		Q=abs(data.vx)
    	elif key=='vy':
    		Q=abs(data.vy)
    	elif key=='vz':
    		Q=abs(data.vz)
    	elif key=='bx':
    		Q=abs(data.bx)
    	elif key=='by':
    		Q=abs(data.by)
    	elif key=='bz':
    		Q=abs(data.bz)
    	    
#   	 Q=data.vz
#   	 Q=data.bz**2*data00.g_33+       sqrt(data.v2)
#   	 Q[data.bx*data.bz*data00.g_13*2.isnan(Q)]=0.
#   	 Q[isinf(Q)]=0.
#   	 Q=data.bx*sqrt(data00.g_11)
#   	 Q=data.by*sqrt(data00.g_22)
#   	 Q=data.bz*sqrt(data00.g_33)
    	    
    	qmin=Qmin
    	qmax=Qmax
    	if (Qmin is None):
    	    qmin=Q.min()
    	if (Qmax is None):
    	    qmax=Q.max()
    	if (log and qmin==0.):
    	    qmin=qmax*1e-10

    	Q[where(Q<qmin)]=qmin
    	Q[where(Q>qmax)]=qmax

    	if case=='disk':
    	    if data.nz>1:
    	        equator=True
    	    else:
    	        equator=False
    	    equator=False
#    	    imm.append(data.image(Q,disk=False,r_hor=data00.r_h,xmax=xmax,log=log,anim=True,alpha=data00.alpha,beta=data00.betaz,P_c=data00.P_c,equator=equator,time=data.time,cb_range=[qmin*0.9999,qmax],fig=fig,ax=ax[switch[key][0][0],switch[key][0][1]],title=key))

        if image:
#### IMAGE ####
    	    imm[switch[key][0]],ax[switch[key][0]]=data.image(Q,disk=False,r_hor=data00.r_h,xmax=xmax,log=log,anim=True,alpha=data00.alpha,beta=data00.betaz,P_c=data00.P_c,equator=equator,time=data.time,cb_range=[qmin*0.9999,qmax],fig=fig,ax=ax[switch[key][0]],title=key+r" $\in$ [%2.0g,%2.0g], Time = %i"%(qmin,qmax,data.time),cbar=False)
        else:
#### CONTOUR ####
    	    imm[switch[key][0]],ax[switch[key][4]]=data.contour(Q,disk=False,log=log,anim=True,cb_range=[qmin*0.9999,qmax],ax=ax[switch[key][4]],title=key+r" $\in$ [%2.0g,%2.0g], Time = %i"%(qmin,qmax,data.time),cbar=False)
    return imm
    
# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

# call the animator.  blit=True means only re-draw the parts that have changed.
anim=animation.FuncAnimation(fig, animate,frames=nout, interval=1000,repeat=False,blit=False)
anim.save(filename=filename+'.mp4',writer=writer)
#anim.save(filename=filename+'.mp4',fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()

