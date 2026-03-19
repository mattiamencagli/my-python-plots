import numpy as np
import powerspectrum as pws
import multiprocessing as mp
import ep_tools as ep
import os


dnucor_max=10.

os.chdir('/scratch/seismo/papini/KEPLER/KIC005184732_kitty')

class chunks_kitty(object):
    """
      it's just an object defining the chunks containing spot's signals, of which to compute 
      Lomb-Scargle Periodograms to use for correlations

      Time=0 of the chunks is assumed to be the beginning of Q7.1

    """
    def __init__(self):
        self.KIC=5184732
        self.n=10  #number of chunks
        self.range=[] #ranges of chunks, manually defined.
        self.range.append([13.,30.])
        self.range.append([37.,50.])
        self.range.append([72.,87.])
        self.range.append([280.,295.])
        self.range.append([380.,400.])
        self.range.append([400.,420.])
        self.range.append([452.,467.])
        self.range.append([555.,595.])
        self.range.append([677.,692.])
        self.range.append([695.,712.])
        self.range.append([727.,742.])
        self.range.append([918.,933.])
        self.range.append([933.,955.])

        self.n=np.shape(self.range)[0]

        self.nl=0 #number of long_chunks
        self.ra_l=[] #ranges of long_chunks, manually defined. FIXED TO 30days
        self.ra_l.append([60.,90.])
        self.ra_l.append([134.,164.])
        self.ra_l.append([372.,402.])
        self.ra_l.append([563.,593.])
        self.ra_l.append([676.1,706.1])
        self.ra_l.append([790.,820.])
        self.ra_l.append([918.,948.])

        self.nl=np.shape(self.ra_l)[0]

        self.pmode_envelope=[1500.,2800.]
        self.pmode_to_corr=[1730.,2300.]
        

#interface for multiprocessing
def mp_spectrum (res,*args,**kwargs):
    res.put(pws.spectrum(*args,**kwargs))
    

#define some parameters

is_concatenated = True


filel=['kplr005184732_kasoc-ts_slc_v1.fits']
filel_llc=['kplr005184732_kasoc-ts_llc_v1.fits']

fitskw='PDCSAP_FLUX'
fitskw='FLUX'

#reading the lightcurve
if is_concatenated :
    import astropy.io.fits as pyfits
    fitsfile=pyfits.open(filel[0])
    time, light= np.array(fitsfile[1].data.field('TIME')), np.array(fitsfile[1].data.field(fitskw))
    fitsfile.close()
else :
    time, light= pws.makelightcurve('PDCSAP_FLUX',filelist=filel)

time =time -time[0]

#removing nans et similia    
light=light[np.isfinite(time )]
time = time[np.isfinite(time )]
time = time[np.isfinite(light)]
light=light[np.isfinite(light)]

#getting the ranges for the chunked lightcurve
chu=chunks_kitty()



kchu=np.where((time >= chu.ra_l[0][0]) & (time <= chu.ra_l[0][1]))

#setting workers for parallel processing
nwork=chu.nl-1
nwork_sub=4

kchu=[]
for i in chu.ra_l : kchu.append(np.where((time >= i[0]) & (time <=i[1]))[0])

ps=pws.spectrum(time[kchu[0]],light[kchu[0]],psd=True,n=nwork_sub)

freq_ref=ps[0].copy()


#setting machinery for parallel processing
res = mp.Queue() #will contain the result
#
workers = [mp.Process(target=mp_spectrum, args=(res,time[kchu[i+1]],light[kchu[i+1]]),
                      kwargs={'freq': freq_ref,'psd': True,'n':nwork_sub}) 
                      for i in range(nwork)]


#calculate spectra of chunks
for each in workers: each.start()
#

out=[ps]
for i in range(nwork):
    out.append(res.get())

for each in workers: each.join()
#
#

for each in workers: each.terminate()


#create COMPLEX intensity spectrum from LS-periodogram
#
#  I(t)=A(om)cos(om*t) + B(om)sin(om*t)
# I(om)=A(om)+iB(om)
out=np.asarray(out) #converting to np.ndarray
freq=out[0,0,:]
I_compl=out[:,2,:] +1j*out[:,3,:]

chu.freq=freq.copy()
chu.pws=out[:,1,:].copy()
chu.I_compl=I_compl.copy()
#correlate spectra
cmc=ep.self_correlate(I_compl,x=freq_ref*1e6,xmax=dnucor_max)

#average the correlations over different spot timeseries powerspectra
cmc_avg=cmc.average()



#load martin fitted frequencies
dummy=np.loadtxt('/scratch/seismo/papini/KEPLER/NIELSEN_table2.dat',dtype='str')
k=np.where(np.float64(dummy[0,1:]) == chu.KIC)
#dictionary of Martin's dataset
#dataset={}
#for i in range(dummy.shape[0]):
#    dataset[dummy[i,0]]=dummy[i,k[0]+1]
#chu.dataset=dataset.copy()

#list of Martin's dataset fit
fit_data=[dummy[:,0],dummy[:,1+k[0][0]]]

#identifying modesets
kfreq_l0=np.where( ['frequency_l0_50th' in s for s in fit_data[0][:].flatten()])[0]
khght_l0=np.where( ['height_l0_50th'    in s for s in fit_data[0][:].flatten()])[0]
kfreq_l1=np.where( ['frequency_l1_50th' in s for s in fit_data[0][:].flatten()])[0]
khght_l1=np.where( ['height_l1_50th'    in s for s in fit_data[0][:].flatten()])[0]
kfreq_l2=np.where( ['frequency_l2_50th' in s for s in fit_data[0][:].flatten()])[0]
khght_l2=np.where( ['height_l2_50th'    in s for s in fit_data[0][:].flatten()])[0]
kwdth   =np.where( ['width_50th'        in s for s in fit_data[0][:].flatten()])[0]
ksplit  =np.where( ['splitting_50th'    in s for s in fit_data[0][:].flatten()])[0]

dnu_split= np.mean(np.float64(fit_data[1][ksplit[:-1]]))

#mode object
modes=type('Mode',(),{'freq':[],'width':[],'height':[],'l':[],'modeset':[]})()

modes.dnu_split=dnu_split
#cycling through modeset
for i in range(int(fit_data[0][-1][2])) : #assumes that the modeset number is 
#                                         #in the third element of the string 
#                                         #of the last line of table2.dat
    #l=0
    modes.freq.append(  np.float(fit_data[1][kfreq_l0[i]]))
    modes.height.append(np.float(fit_data[1][khght_l0[i]]))
    modes.width.append( np.float(fit_data[1][kwdth[i]] ))
    modes.l.append(0)
    modes.modeset.append(i)
    #cmc.knui_mode.append(ep.find_nearest(cmc.nui,modes.freq[-1],index=True))
    #l=1
    modes.freq.append(  np.float(fit_data[1][kfreq_l1[i]]))
    modes.height.append(np.float(fit_data[1][khght_l1[i]]))
    modes.width.append( np.float(fit_data[1][kwdth[i]] ))
    modes.l.append(1)
    modes.modeset.append(i)
    #l=2
    modes.freq.append(  np.float(fit_data[1][kfreq_l2[i]]))
    modes.height.append(np.float(fit_data[1][khght_l2[i]]))
    modes.width.append( np.float(fit_data[1][kwdth[i]] ))
    modes.l.append(2)
    modes.modeset.append(i)    
    #cmc.knui_mode.append(ep.find_nearest(cmc.nui,modes.freq[-1],index=True))
modes.shape=np.shape(modes.freq)    
#correlating l=1 and l=2 multiplets
#selecting multiplet ranges
dnucorr=10.#dnucor_max
modes.ranges=np.ndarray((np.size(modes.freq),2) ,dtype=np.int32)

#finding the modes indexes
for i,jfreq in enumerate(modes.freq):
    if jfreq > 0.:
        modes.ranges[i,:]=np.where(
                          np.logical_and(cmc.nui >= jfreq -dnucorr 
                                        ,cmc.nui<= jfreq +dnucorr))[0][[0,-1]]
    else:
        modes.ranges[i,:]=[-1,-1]

# forcing ranges to same size
modes.ranges[:,1]=modes.ranges[:,0] + np.median(modes.ranges[:,1]-modes.ranges[:,0])

n_inrange=modes.ranges[0,1]-modes.ranges[0,0]

l1_mask=np.where([i == 1 for i in modes.l])[0]
l2_mask=np.where([i == 2 for i in modes.l])[0]
l1_mask=l1_mask[np.where(modes.ranges[l1_mask,0] > -1)]
l2_mask=l2_mask[np.where(modes.ranges[l1_mask,0] > -1)]


cmc.cmc_chunks=np.zeros(modes.shape+(n_inrange,cmc_avg.shape[1],),dtype=np.complex128)

cmc.cmc_over_chunks=np.zeros(modes.shape+(cmc_avg.shape[1],))

cmc.cmc_over_nl=np.zeros((2,modes.ranges[0,1]-modes.ranges[0,0],cmc_avg.shape[1]),dtype=np.complex128)


lmax=max(modes.l)

#getting chunks of correlations
cmc.get_chunks(Range=modes.ranges)
dnui=cmc.nui[modes.ranges.flatten()[0]:modes.ranges.flatten()[1]]
dnui-=dnui[n_inrange/2]
####NOW I HAVE
#
#   cmc.cmc_chunks[nchunks,nmodes,nui_in_range,nucor]
#
#   these can be averaged in different ways to see if there is a signal
#
# WHAT FOLLOWS COULD BE DONE MUCH MORE EFFICIENTLY (TIME CONSUMING AND MEMORY CONSUMING)
# BUT SINCE IT IS FAST ANYWAY I DO IT LIKE THIS TO MAKE LESS MISTAKES


#    --1st--    avg=abs(sum_nui, sum_chunks  ...)^2       [nmodes,nucor]
#               avg[i,:]=avg[i,:]/max(avg[i,:]) for i=1,N [nmodes,nucor]
#               plot(np.sum(avg,0)/N)                     [nucor]
#
#    --2nd--
#


#--FIRST--#
#autocorrelation of each average multiplet (i.e. autocorrelation of the 
# chunk averaged correlations <I(nu)*I(nu+dnu)>_chunks
# over the frequency interval encompassing the multiplet)
#
#  cmc.cmc_over_chunks[nl,nu+dnu]=
#       sum_[nu_nl] |<I(nu_nl)*I(nu_nl+dnucor)>_chunks|**2
#  normalized to one

cmc_chu_avg=np.sum(cmc.cmc_chunks,axis=0) #<I(nu_nl)*I(nu_nl+dnucor)>_chunks

cmc_chu_nui_avg=np.square(np.abs(cmc_chu_avg)) #|<I(nu_nl)*I(nu_nl+dnucor)>_chunks|^2
cmc_chu_nui_avg=cmc_chu_nui_avg.sum(axis=1) #\int_[nl] dnu |<I(nu_nl)*I(nu_nl+dnucor)>_chunks|^2

modes.cmc_norms=np.ndarray(modes.shape[0])
for i, j in enumerate(filter(lambda x: x>=0,modes.ranges[:,0])) :
    modes.cmc_norms[i]=np.max(cmc_chu_nui_avg[i,:])     
    cmc_chu_nui_avg[i,:]/=np.max(cmc_chu_nui_avg[i,:]) 

import matplotlib.pylab as plt

plt.figure(1)
    
plt.plot(cmc.dnucor,cmc_chu_nui_avg[l1_mask,:].sum(axis=0)/l1_mask.shape[0],linewidth=1,label='$\ell=1$')
plt.plot(cmc.dnucor,cmc_chu_nui_avg[l2_mask,:].sum(axis=0)/l2_mask.shape[0],linewidth=1,label='$\ell=2$')
plt.xlim((-5,5))
plt.xlabel(r'$\Delta\nu /\mu$Hz')
plt.ylabel(r'$\sum_{n=1}^N \frac{1}{N \int_[nl]<P(\nu)>d\nu}\int_{[n\ell]}|<I(\nu)*I(\nu+\Delta\nu)>|^2$d$\nu$')
[plt.axvline(x=i*dnu_split, color='red',alpha=0.5,linestyle='--') for i in np.arange(9)-4]
plt.legend()       


#--SECOND--#
#autocorrelation of each average multiplet (i.e. autocorrelation of the 
# chunk averaged correlations <I(nu)*I(nu+dnu)>_chunks
# over the frequency interval encompassing the multiplet)
#
#  cmc.cmc_over_chunks[nl,nu+dnu]=
#       sum_[nu_nl] |<I(nu_nl)*I(nu_nl+dnucor)>_chunks|**2
#  normalized to one
cmc_chu_modes_avg=np.square(np.abs(cmc_chu_avg)) #|<I(nu_nl)*I(nu_nl+dnucor)>_chunks|^2
cmc_chu_modes_avg=cmc_chu_nui_avg.sum(axis=1) #\int_[nl] dnu |<I(nu_nl)*I(nu_nl+dnucor)>_chunks|^2

#\sum_modes<I(nu_nl)*I(nu_nl+dnucor)>_chunks*1/C_[pnl](nu)
cmc_chu_l1_avg=np.sum(cmc_chu_avg[l1_mask,:,:].transpose()/np.sqrt(modes.cmc_norms[l1_mask]),axis=2) 
cmc_chu_l2_avg=np.sum(cmc_chu_avg[l2_mask,:,:].transpose()/np.sqrt(modes.cmc_norms[l2_mask]),axis=2) 

#plotting int_[nl]dnu |<I(nu)*I(nu+dnu)>_chunks_modes|^2
plt.figure(2)
plt.plot(cmc.dnucor,np.sum(np.square(np.abs(cmc_chu_l1_avg)),axis=1),label='l1')
plt.plot(cmc.dnucor,np.sum(np.square(np.abs(cmc_chu_l2_avg)),axis=1),label='l2')
plt.legend()



#PLOTTING SLICES AT \Delta m=+-m\omega_beta to seek for some correlation
#using rot splitting from martin

f, axarr = plt.subplots(5,lmax)

for j in range(lmax*4+1) :
    idnucor=np.where(cmc.dnucor >= dnu_split*(j-2*lmax))[0][0]
    [axarr.flatten()[j].plot(dnui,cmc_chu_avg[i,:,idnucor],label="$\ell=1, MS=$"+str(modes.modeset[i])) for i in l1_mask]
    axarr.flatten()[j].plot(dnui,np.abs(cmc_chu_avg[l1_mask,:,idnucor].sum(axis=0)),label="$\ell=1$, abs(Total)",linewidth=2)
    axarr.flatten()[j].set_title('$\Delta m = $'+str(j-2*lmax))
    #axarr.flatten()[j].legend()
axarr.flatten()[1].legend(bbox_to_anchor=[1.25,1])

f, axarr = plt.subplots(5,lmax)

for j in range(lmax*4+1) :
    idnucor=np.where(cmc.dnucor >= dnu_split*(j-2*lmax))[0][0]
    [axarr.flatten()[j].plot(dnui,cmc_chu_avg[i,:,idnucor],label="$\ell=2, MS=$"+str(modes.modeset[i])) for i in l2_mask]
    axarr.flatten()[j].plot(dnui,np.abs(cmc_chu_avg[l2_mask,:,idnucor].sum(axis=0)),label="$\ell=2$, abs(Total)",linewidth=2)
    axarr.flatten()[j].set_title('$\Delta m = $'+str(j-2*lmax))
    #axarr.flatten()[j].legend()
axarr.flatten()[1].legend(bbox_to_anchor=[1.25,1])

"""
selezionare le parti di cmc da correlare, fare il quadrato e sommare
test differenti modi di fare la media, scriverli tutti qui
"""

#OLD /ALTERNATIVE WAY OF DOING THE FIRST
if False :
    for i,j in enumerate(modes.ranges):
        if j[0] >= 0 :
            #autocorrelation of each average multiplet (i.e. autocorrelation of the 
            # chunk averaged correlations <I(nu)*I(nu+dnu)>_chunks
            # over the frequency interval encompassing the multiplet)
            #
            #  cmc.cmc_over_chunks[nl,nu+dnu]=
            #       sum_[nu_nl] |<I(nu_nl)*I(nu_nl+dnucor)>_chunks|**2
            #  normalized to one
    
            cmc.cmc_over_chunks[i,:]=np.sum(np.square(np.abs(cmc_avg[j[0]:j[1],:])),axis=0)
            cmc.cmc_over_chunks[i,:]/=np.max(cmc.cmc_over_chunks[i,:])
    
            if modes.l[i] == 1 : cmc.cmc_over_nl[0,:,:]+=cmc_avg[j[0]:j[1],:]
            if modes.l[i] == 2 : cmc.cmc_over_nl[1,:,:]+=cmc_avg[j[0]:j[1],:]
    #calculating <I(nu)*I(nu+dnu)>_chunks_nl and then the autocorrelation of those
    
    
    
    
    import matplotlib.pylab as plt
    #plotting average correlation over all the multiplets (l=1 and l=2)
    plt.figure(1)
    
    plt.plot(cmc.dnucor,cmc.cmc_over_chunks[l1_mask,:].sum(axis=0)/l1_mask.shape[0],linewidth=1,label='$\ell=1$')
    plt.plot(cmc.dnucor,cmc.cmc_over_chunks[l2_mask,:].sum(axis=0)/l2_mask.shape[0],linewidth=1,label='$\ell=2$')
    plt.xlim((-5,5))
    plt.xlabel(r'$\Delta\nu /\mu$Hz')
    plt.ylabel(r'$\sum_{n=1}^N \frac{1}{N \int_[nl]<P(\nu)>d\nu}\int_{[n\ell]}|<I(\nu)*I(\nu+\Delta\nu)>|^2$d$\nu$')
    [plt.axvline(x=i*dnu_split, color='red',alpha=0.5,linestyle='--') for i in np.arange(9)-4]
    plt.legend()       
