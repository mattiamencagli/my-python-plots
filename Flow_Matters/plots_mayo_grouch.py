import numpy as np
import h5py as h5
import matplotlib
import matplotlib.pyplot as plt


LOG_DISTRIB = 1;
# options for filter out large droplets (R>RF)
FILTER = 1; 
RF = 80;

save_dir = '/home/matti/Pictures/image_test/'

cases={'256'       : {'path':'/home/matti/DATA/Mayonese/staggered_run_halfZ/', 'Nz':256, 'col':'#26547c', 'linestyle':'-', 'linewidth':2, 'marker':'s', 'marker2':'P', 'markersize':60},
       '512'       : {'path':'/home/matti/DATA/Mayonese/staggered_run/'      , 'Nz':512, 'col':'#ef476f', 'linestyle':'-', 'linewidth':2, 'marker':'o', 'marker2':'+', 'markersize':60},
       '768'       : {'path':'/home/matti/DATA/Mayonese/staggered_run_1.5Z/' , 'Nz':768, 'col':'#ffd166', 'linestyle':'-', 'linewidth':2, 'marker':'v', 'marker2':'x', 'markersize':60},
       }
iEND = len(cases)-1
int_keys = np.array([int(key) for key in cases.keys()])

#######################################################################################################################
# primo e secondo elemento: step e number_droplets
# terzo e quarto elemento: avarage volume e total volume [non per il file posvel]
# il resto sono vol e surf delle droplets
def read_nonunif_col(f,bool_twotwo,bool_vol):
    first_two_columns = []
    second_two_columns = []
    matrix_data = []
    with open(f, 'r') as file:
        for line in file:
            elements = line.split()
            if bool_twotwo:
                first_two_columns.append([int(elements[0]), int(elements[1])])
                second_two_columns.append([float(elements[2]), float(elements[3])]) 
                data = list(map(float, elements[4:]))
            else:
                first_two_columns.append([int(elements[0]), int(elements[1])])
                data = list(map(float, elements[2:]))
            matrix_data.append(data)
    max_length = max(len(row) for row in matrix_data)
    matrix_data_padded = np.array([row + [0] * (max_length - len(row)) for row in matrix_data])
    if bool_twotwo and bool_vol:
        #matrix_data_padded = ((3.0/(4.0*np.pi))*matrix_data_padded)**(1.0/3.0)
        second_two_columns = ((3.0/(4.0*np.pi))*np.array(second_two_columns))**(1.0/3.0) # from volume to radius
    return np.array(first_two_columns).T, matrix_data_padded, np.array(second_two_columns).T
#######################################################################################################################   
def comp_rad(V):
    return np.cbrt(3.0*V/(4.0*np.pi))
#######################################################################################################################   
def compute_Rava_D32_D33(steps,V,F,_RF):
    VF = 1e15
    if F: VF = 4*np.pi*(_RF**3)/3
    Vava = np.array([np.mean(V[i], where=((V[i]>0)&(V[i]<VF))) for i in range(steps)])
    D = 2*comp_rad(V)
    if F: D[D>_RF*2] = 0 # rimuovo colonne
    D32num   = np.array([np.sum(D[i]**3) for i in range(steps)])
    D32denum = np.array([np.sum(D[i]**2) for i in range(steps)])
    D[D==0] = 1 # cambio gli elementi nulli per non far sbroccare il logaritmo
    D33num   = np.array([np.sum(np.log(D[i])*D[i]**3) for i in range(steps)])
    D33denum = D32num
    # ritorno il Raggio con media aritmetica e il raggio (quindi *0.5) calcolato con D32 e D33
    return comp_rad(Vava), (0.5*D32num/D32denum), (0.5*np.exp(D33num/D33denum)) 
#######################################################################################################################
def read_field_fromfile(name, field, Nz, Ny, Nx):
    file = h5.File(name, "r")
    Nx = int(file.get("LBE3D").attrs["lbe_sx"][0])
    Ny = int(file.get("LBE3D").attrs["lbe_sy"][0])
    Nz = int(file.get("LBE3D").attrs["lbe_sz"][0])
    name_field="/LBE3D/"+field
    f = np.zeros((Nx,Ny,Nz))
    f[:,:,:] = file.get(str(name_field))[...]
    file.close()
    return f
#######################################################################################################################
def read_SST_antidiagonal(direc, Nz, Ny, Nx):
    SST = np.zeros((3,Nx,Ny,Nz))
    SST[0,:,:,:] = read_field_fromfile(direc+"SST_xy.600001.h5","SST_xy", Nz, Ny, Nx)
    SST[1,:,:,:] = read_field_fromfile(direc+"SST_xz.600001.h5","SST_xz", Nz, Ny, Nx) # <---------------------- !!!!!
    SST[2,:,:,:] = read_field_fromfile(direc+"SST_yz.600001.h5","SST_yz", Nz, Ny, Nx)
    return SST
#######################################################################################################################
def SST_average_alongZ(SST, Nz, Ny, Nx):
    SST_ave = np.zeros((3,Nz))
    SST_ave = np.sum(SST, axis=(1,2))/(Nx*Ny)
    return SST_ave
#######################################################################################################################
def read_posvel(f):
    step_e_num, d_posvel, _ = read_nonunif_col(f, 0, 0)
    d_pos = np.zeros((len(d_posvel[:,0]),int(len(d_posvel[0,:])/2)))
    d_vel = np.zeros_like(d_pos)
    for j in range(len(d_posvel[:,0])):
        d_pos[j,:] = np.array([d_posvel[j,i] for i in range(len(d_posvel[j,:])) if i % 6 <  3])
        d_vel[j,:] = np.array([d_posvel[j,i] for i in range(len(d_posvel[j,:])) if i % 6 >= 3])
    return d_pos, d_vel
#######################################################################################################################
def comp_distr_radii(direc, Nz, RF):
    name = "droplets_volume_rho1.dat"
    step_e_num, d_vol, avg_e_tot_vol = read_nonunif_col(direc+name, 1, 1)
    bool_direc = direc!="/home/matti/DATA/Mayonese/staggered_run/"
    if bool_direc: 
        name = "droplets_posvel_rho1.dat"
        DVOL = d_vol
    else: 
        name = "droplets_posvel_rho1_ALT.dat"
        name1 = "droplets_volume_rho1_ALT.dat"
        step_e_num_ALT, DVOL, avg_e_tot_vol_ALT = read_nonunif_col(direc+name1, 1, 1)
    d_pos, d_vel = read_posvel(direc+name)
    idx_sz = 0; idx_ez = Nz; stepz = 32
    RR_RdistrGAP = range(idx_sz,idx_ez,stepz)
    distr_Rava = np.zeros((len(d_pos[:,0]),int(Nz/stepz)))
    distr_RD32 = np.zeros_like(distr_Rava)
    distr_RD33 = np.zeros_like(distr_Rava)
    z = d_pos[:,2::3]
    j = 0
    if FILTER: 
        VF = 4*np.pi*(RF**3)/3
    else: 
        VF = 1e15
        RF = 1e15
    for i in RR_RdistrGAP:
        maskz = (z>i) & (z<=i+stepz)
        masked_d_vol = np.where(maskz, DVOL, 0)
        masked_d_diam = 2*comp_rad(masked_d_vol) # use diameters for D32 and D33
        tot_drop = np.sum(maskz, axis=1, where=(masked_d_vol<VF))
        tot_drop[tot_drop==0] = 1 # remove zeros for division
        distr_Rava[:,j] = np.sum(masked_d_vol, axis=1, where=(masked_d_vol<VF))/tot_drop
        D32num = np.sum(masked_d_diam**3, axis=1, where=(masked_d_diam<RF*2))
        D32denum = np.sum(masked_d_diam**2, axis=1, where=(masked_d_diam<RF*2))
        masked_d_diam[masked_d_diam==0] = 1  # remove zeros for logarithm
        D33num = np.sum(np.log(masked_d_diam)*masked_d_diam**3, axis=1, where=(masked_d_diam<RF*2))
        D33denum = D32num
        D32denum[D32denum==0] = 1 # remove zeros for division
        D33denum[D33denum==0] = 1 # remove zeros for division
        distr_RD32[:,j] = D32num/D32denum
        distr_RD33[:,j] = D33num/D33denum
        j += 1
    distr_Rava = ((3.0/(4.0*np.pi))*np.array(distr_Rava))**(1.0/3.0)
    distr_RD32 = 0.5*distr_RD32
    distr_RD33 = 0.5*np.exp(distr_RD33)
    distr_Rava_ret = np.repeat(distr_Rava[-1,:],stepz)
    distr_RD32_ret = np.repeat(distr_RD32[-1,:],stepz)
    distr_RD33_ret = np.repeat(distr_RD33[-1,:],stepz)
    return distr_Rava_ret, distr_RD32_ret, distr_RD33_ret




for i,case in enumerate(cases):

    ls   = cases[case]['linestyle']
    lw   = cases[case]['linewidth']
    col  = cases[case]['col']
    mk   = cases[case]['marker']
    mk2  = cases[case]['marker2']
    ms   = cases[case]['markersize']
    Nx   = 512
    Ny   = 512
    Nz   = cases[case]['Nz']
    path = cases[case]['path']

    #read Vfrac and Ndroplets
    name  = "output.csv"
    #head  = np.genfromtxt(path+name, delimiter=',\t', dtype=str,   max_rows=1); print(f"header = {head}")
    Vfrac = np.genfromtxt(path+name, delimiter=',\t', dtype=float, skip_header=True)
    
    # volume fraction overtime
    _ = plt.figure(1)
    _ = plt.plot(Vfrac[:,0], Vfrac[:,2], color=col, linestyle=ls, linewidth=lw, label="Gap = %s; final $V_{frac}$=%1.1f%%"%(case,Vfrac[-1,2]))
    if i==iEND:
        _ = plt.rc('axes', axisbelow=True)
        _ = plt.grid(alpha=0.2)
        _ = plt.xlabel("Steps")
        _ = plt.ylabel("$V_{droplet}$ / $V_{tot}$")
        equation = r"$\rho_{oil} = \rho_{oil}\cdot(1+f_{push})$" + "\n" + r"$\rho_{H_2O} = \rho_{H_2O}\cdot(1-f_{push}\cdot(\rho_{oil}/\rho_{H_2O}))$"
        box = dict(facecolor='white', edgecolor='blue', boxstyle='round,pad=0.3')
        #_ = plt.text(0.05, 0.80, equation, verticalalignment='bottom', horizontalalignment='left', transform=plt.gca().transAxes, fontsize=14, color='black', bbox=box)
        _ = plt.text(0.95, 0.05, equation, verticalalignment='bottom', horizontalalignment='right', transform=plt.gca().transAxes, fontsize=12, color='black', bbox=box)
        _ = plt.legend()

    # number of droplets overtime
    _ = plt.figure(2)
    _ = plt.plot(Vfrac[:,0], Vfrac[:,1], color=col, linestyle=ls, linewidth=lw, label="Gap = %s; final $N_{droplets}$=%d"%(case,Vfrac[-1,1]))
    if i==iEND:
        _ = plt.rc('axes', axisbelow=True)
        _ = plt.grid(alpha=0.2)
        _ = plt.xlabel("Steps")
        _ = plt.ylabel("$N_{droplets}$")
        _ = plt.legend()

    # compute D32
    name = "droplets_volume_rho1.dat"
    step_e_num, d_vol, avg_e_tot_vol = read_nonunif_col(path+name, 1, 0)
    Rava, D32, D33 = compute_Rava_D32_D33(len(step_e_num[0,:]), d_vol, FILTER, RF)
    
    # D32 overtime
    if FILTER: labF = " [filter R>%d]"%RF
    else: labF= ""
    _ = plt.figure(3)
    _ = plt.plot(step_e_num[0,:], D32, color=col, linestyle=ls, linewidth=lw, label=r"Gap = %s; final $0.5*D_{3,2}$=%1.1f"%(case,D32[-1]))
    if i==iEND:
        _ = plt.rc('axes', axisbelow=True)
        _ = plt.grid(alpha=0.2)
        _ = plt.xlabel('Steps')
        _ = plt.ylabel('Radius')
        _ = plt.title(r'$0.5*D_{3,2}$'+labF)
        _ = plt.legend()

    # final D32 as function of gap size
    _ = plt.figure(4)
    if i==0:
        _ = plt.plot(int_keys, np.array([D32[-1], D32[-1]*2, D32[-1]*3]), color='black', linewidth=1, linestyle='--')
    _ = plt.scatter(Nz, D32[-1], color=col, marker=mk, s=ms, label="Gap = %s"%case)
    if i==iEND:
        _ = plt.rc('axes', axisbelow=True)
        _ = plt.grid(alpha=0.2)
        _ = plt.ylabel('Radius')
        _ = plt.xlabel('Gap size')
        _ = plt.xticks(int_keys)
        _ = plt.legend()

    #read SST and average in XY plane
    SST = read_SST_antidiagonal(path, Nz, Ny, Nx)
    SSTave = SST_average_alongZ(SST, Nz, Ny, Nx)

    # shearstress along the gap (figures 5,6,7)
    if i==iEND: Nticks=7
    else: Nticks=5
    XT = np.linspace(0,Nz,Nticks)
    _ = plt.figure(5+i)
    _ = plt.plot(SSTave[1,:], color='red'  , linestyle=ls, linewidth=lw, label=r"$P_{xz}$") #important one
    _ = plt.plot(SSTave[0,:], color='blue' , linestyle=ls, linewidth=lw, label=r"$P_{xy}$")
    _ = plt.plot(SSTave[2,:], color='green', linestyle=ls, linewidth=lw, label=r"$P_{yz}$")
    aveP = np.average(SST[1,:])
    _ = plt.hlines(aveP, 0, Nz, color='black', linestyle = '--', linewidth=1, label=r"$\overline{P_{xz}}$=%1.5f"%aveP)
    # stdP = np.sqrt(np.var(f[1,10:-10]))
    # _ = plt.hlines([aveP+stdP,aveP-stdP], 0, N, color='black', linewidth=0.5, linestyle = '--')
    # _ = plt.fill_between(XT, aveP-stdP, aveP+stdP, color='yellow', alpha=0.2, label='±1 Std Dev Region')
    _ = plt.title(r"Gap = %s"%case)
    _ = plt.rc('axes', axisbelow=True)
    _ = plt.grid(alpha=0.2)
    _ = plt.ylabel('Shear stress')
    _ = plt.xlabel('Gap [Z axis]')
    _ = plt.xticks(XT)
    _ = plt.legend()
    _ = plt.savefig(save_dir+"SST_"+case+".pdf", bbox_inches='tight', transparent=True, format='pdf')

    #compute distribution D32
    distr_Rava, distr_D32, distr_D33 = comp_distr_radii(path, Nz, RF)
    #compute Ca from distribution
    sigma = 0.0235
    Ca_new = SSTave[1,:]*distr_D32/sigma
    aveCa_new = np.average(Ca_new)
    #read dimensionless number computed throught vrms
    name = "dimensionless_numbers"+case+".csv"
    dimless = np.loadtxt(path+name, delimiter=',', skiprows=1).T
    aveCa_old = dimless[3,-1]

    #comparison at the final steps of Ca (with vrms versus with SST)
    _ = plt.figure(8)
    _ = plt.scatter(Nz, aveCa_new, color=col, marker=mk,  s=ms, label="Gap = %s (with $SST$)"%case)
    _ = plt.scatter(Nz, aveCa_old, color=col, marker=mk2, s=ms, label="Gap = %s (with $u_{rms}$)"%case)
    if i==iEND:
        _ = plt.rc('axes', axisbelow=True)
        _ = plt.grid(alpha=0.2)
        _ = plt.title('Ca comparison along the gap at the end')
        _ = plt.ylabel(r'Capillarity number $Ca$')
        _ = plt.xlabel('Gap size')
        _ = plt.xticks(int_keys)
        _ = plt.legend()

    # Ca overtime
    _ = plt.figure(9)
    _ = plt.plot(dimless[0,:], dimless[3,:], color=col, linestyle=ls, linewidth=lw, label=r'Gap = %s; final $Ca$=%1.2f'%(case,dimless[3,-1]))
    if i==iEND:
        _ = plt.rc('axes', axisbelow=True)
        _ = plt.grid(alpha=0.2)
        _ = plt.xlabel('Steps')
        _ = plt.ylabel(r'Capillarity number $Ca$')
        _ = plt.title(r"$Ca = \rho_{tot} \cdot u_{rms}\cdot \nu / \sigma$")
        _ = plt.legend()

    # Re overtime
    _ = plt.figure(10)
    _ = plt.plot(dimless[0,:], dimless[1,:], color=col, linestyle=ls, linewidth=lw, label=r'Gap = %s; final $Re$=%1.2f'%(case,dimless[1,-1]))
    if i==iEND:
        _ = plt.rc('axes', axisbelow=True)
        _ = plt.grid(alpha=0.2)
        _ = plt.xlabel('Steps')
        _ = plt.ylabel(r'Reynolds number $Re$')
        _ = plt.title(r"$Re = u_{rms}\cdot N_z / \nu$")
        _ = plt.legend()

    print(f"case {case} done!")


_ = plt.figure(1)
_ = plt.savefig(save_dir+"V_frac.pdf", bbox_inches='tight', transparent=True, format='pdf')
_ = plt.figure(2)
_ = plt.savefig(save_dir+"N_drop.pdf", bbox_inches='tight', transparent=True, format='pdf')
_ = plt.figure(3)
_ = plt.savefig(save_dir+"D32.pdf", bbox_inches='tight', transparent=True, format='pdf')
_ = plt.figure(4)
_ = plt.savefig(save_dir+"D32_final_scale.pdf", bbox_inches='tight', transparent=True, format='pdf')
_ = plt.figure(8)
_ = plt.savefig(save_dir+"Ca_old_new_comparison.pdf", bbox_inches='tight', transparent=True, format='pdf')
_ = plt.figure(9)
_ = plt.savefig(save_dir+"Ca.pdf", bbox_inches='tight', transparent=True, format='pdf')
_ = plt.figure(10)
_ = plt.savefig(save_dir+"Re.pdf", bbox_inches='tight', transparent=True, format='pdf')
    