import multiprocessing as mp ; import sys

nwork = int(sys.argv[1])

ti = int(sys.argv[2])         #tempo iniziale
tf = int(sys.argv[3])         #tempo finale
d  = int(sys.argv[4])         #ogni quanti tempi leggo?

tempi = int((tf-ti)/d) + 1

if tempi != nwork : 
    print('ATTENZIONE:numero di output da analizzare diverso dal numero di processori dati') ; print('STOP')
    exit()

a=[]
for i in range(nwork):
    a.append('')

for i in range(nwork):
    if d*i+ti<1000: a[i] = (str(d*i+ti))
    if d*i+ti<100:  c = str(d*i+ti) ; a[i] = ('0'+c)
    if d*i+ti<10:   c = str(d*i+ti) ; a[i] = ('00'+c) 
    print(a[i])

def analysis_cineca(a):

    print('START')
    
    import numpy as np ; from scipy.io import FortranFile #import matplotlib as mplt 
    import time ; import h5py ; import funzioni_matti_CINECA as matti 
    
    start=time.time()
    
    #dir_lavoro='/Users/Mattia/Desktop/risultati LDZ/256^3 (init_turb. va=0.6 tmax=40)/'
    #dir_save='/'
    
    dir_lavoro='/marconi_scratch/userexternal/mmencagl/LDZ3D/'
    dir_save='/marconi/home/userexternal/mmencagl/risultati_analisi/512^3(va=0.8,tmax=16)/'
    dir_save='./'
    
    print('dir work : ',dir_lavoro)
    print('dir save : ',dir_save)
    
    N=1024        #quanti processori ho usato
    
    n=512
    nnz=512
    
    array=np.zeros((nnz,n,n,12),dtype=np.float32) #ha bisogno di annullarlo ogni volta per motivi poco chiari.....
    
    print('---LOADING--- out%s.dat' %a)
    f=FortranFile(dir_lavoro+'out%s.dat' %a)
    
    data0=f.read_reals(dtype=np.float32)
    t_out=data0[0]
     
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
       
    #rh.append( array[:,:,:,0] )
    vx = array[:,:,:,1] ; print('---LOAD--- vx')
    vy = array[:,:,:,2] ; print('---LOAD--- vy')
    vz = array[:,:,:,3] ; print('---LOAD--- vz')
    #pe.append( array[:,:,:,4] )
    #5,6,7,8 ? tre campi elettrici e l'entropia?
    bx = array[:,:,:,9] ;  print('---LOAD--- bx')
    by = array[:,:,:,10] ; print('---LOAD--- by')
    bz = array[:,:,:,11] ; print('---LOAD--- bz')
    #btot.append( np.sqrt(bx[i]**2 + by[i]**2 + (bz[i]-1.)**2) )
    #vtot.append( np.sqrt(vx[i]**2 + vy[i]**2 + (vz[i]-0.)**2) )
    #gamma.append( (1.0-vtot[i])**(-1/2) )
    
    f.close()
    
    del(data0,data01,data02,data03,ix1,ix2,iy1,iy2,iz1,iz2,arr,ixmin,ixmax,array,n,p,nx,ny,nz,nv,nnz,N)
    
    ''' B e V ############################################## '''
    print('building : BfM2 and VfM2')
    
    BfM2 = (np.absolute(np.fft.rfftn(np.transpose(bx,(2,1,0)))))**2 + (np.absolute(np.fft.rfftn(np.transpose(by,(2,1,0)))))**2 + (np.absolute(np.fft.rfftn(np.transpose(bz,(2,1,0)))))**2 
    VfM2 = (np.absolute(np.fft.rfftn(np.transpose(vx,(2,1,0)))))**2 + (np.absolute(np.fft.rfftn(np.transpose(vy,(2,1,0)))))**2 + (np.absolute(np.fft.rfftn(np.transpose(vz,(2,1,0)))))**2 
    
    spettro_B_gyro , kparal, kperp = matti.gyrotropic_spectrum_3D(BfM2)
    spettro_B_iso, spettro_B_perp, spettro_B_parall = matti.spettri_ridotti(BfM2,spettro_B_gyro)
    
    spettro_V_gyro , kparal, kperp = matti.gyrotropic_spectrum_3D(VfM2)
    spettro_V_iso, spettro_V_perp, spettro_V_parall = matti.spettri_ridotti(VfM2,spettro_V_gyro)
    
    rms_bx = matti.Varianza(bx) ; rms_by = matti.Varianza(by) ; rms_bz = matti.Varianza(bz)
    rms_vx = matti.Varianza(vx) ; rms_vy = matti.Varianza(vy) ; rms_vz = matti.Varianza(vz)
    	
    ''' J ############################################## '''
    
    Jx , Jy , Jz = matti.Curl_migliorissimo(bx,by,bz)
    
    max_Jx=(np.amax(Jx)) ; max_Jy=(np.amax(Jy)) ; max_Jz=(np.amax(Jz))
    
    rms_Jx = matti.Varianza(Jx) ; rms_Jy = matti.Varianza(Jy) ; rms_Jz = matti.Varianza(Jz)
    
    ''' stampo Hdf5 ############################################## '''
    hf = h5py.File(dir_save+'out_analysis%s.h5' %a,'w')
    hf.create_dataset('spettro_B_gyro', data=np.float32(spettro_B_gyro))
    hf.create_dataset('spettro_B_iso', data=np.float32(spettro_B_iso))
    hf.create_dataset('spettro_B_perp', data=np.float32(spettro_B_perp))
    hf.create_dataset('spettro_B_parall', data=np.float32(spettro_B_parall))
    hf.create_dataset('spettro_V_gyro', data=np.float32(spettro_V_gyro))
    hf.create_dataset('spettro_V_iso', data=np.float32(spettro_V_iso))
    hf.create_dataset('spettro_V_perp', data=np.float32(spettro_V_perp))
    hf.create_dataset('spettro_V_parall', data=np.float32(spettro_V_parall))
    hf.create_dataset('kparal', data=np.float32(kparal))
    hf.create_dataset('kperp', data=np.float32(kperp))
    hf.create_dataset('max_Jx', data=np.float32(max_Jx))
    hf.create_dataset('max_Jy', data=np.float32(max_Jy))
    hf.create_dataset('max_Jz', data=np.float32(max_Jz))
    hf.create_dataset('rms_bx', data=np.float32(rms_bx))
    hf.create_dataset('rms_by', data=np.float32(rms_by))
    hf.create_dataset('rms_bz', data=np.float32(rms_bz))
    hf.create_dataset('rms_vx', data=np.float32(rms_vx))
    hf.create_dataset('rms_vy', data=np.float32(rms_vy))
    hf.create_dataset('rms_vz', data=np.float32(rms_vz))
    hf.create_dataset('rms_Jx', data=np.float32(rms_Jx))
    hf.create_dataset('rms_Jy', data=np.float32(rms_Jy))
    hf.create_dataset('rms_Jz', data=np.float32(rms_Jz))
    hf.create_dataset('t_out', data=np.float32(t_out))
    hf.close()   
    print('#################### out_analysis%s.h5 ####################' %a)
    
    end=time.time()
    
    print()
    print('  -------------------------')
    print('  |   END  OF  ANALYSIS   |')
    print('  -------------------------')
    matti.timing(start,end)


print('prima di workers')
# mp.Queue() e mp.join() servono a passare le informazioni tra i processori. (CREDO)
workers = [mp.Process(target=analysis_cineca, args=(a[i],)) for i in range(nwork)]
print('dopo di workers')
for each in workers: each.start()
#print('tra start e join')
for each in workers: each.join()
#print('tra join e terminate')
#for each in workers: each.terminate()
print('FINEFINEFINEFINEFINEFINEFINEFINEFINEFINEFINEFINEFINEFINEFINE')

