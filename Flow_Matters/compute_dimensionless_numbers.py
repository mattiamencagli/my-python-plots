import h5py as h5
import numpy as np
import os
import re
import sys

# primo e secondo elemento: step e number_droplets
# terzo e quarto elemento: avarage volume e total volume [non per il file posvel]
# il resto sono vol e surf delle droplets
def read_Rava(f):
    C = 4.0/(3.0*np.pi)
    c13 = 1.0/3.0
    els = []
    with open(f, 'r') as file:
        for line in file:
            elements = line.split()
            els.append([int(elements[0]), (C*float(elements[2]))**c13])
    return np.array(els)

def read_dim_fromfile(name):
    file = h5.File(name, "r")
    Nx = int(file.get("LBE3D").attrs["lbe_sx"][0])
    Ny = int(file.get("LBE3D").attrs["lbe_sy"][0])
    Nz = int(file.get("LBE3D").attrs["lbe_sz"][0])
    dt = int(file.get("LBE3D").attrs["lbe_diag_nsteps"][0])
    print("Nx Ny Nz = %d %d %d"%(Nx,Ny,Nz), file=sys.stderr)
    print("dt : %d"%(dt), file=sys.stderr)
    return Nx,Ny,Nz,dt
    file.close()

def find_outs(folder):
    files = os.listdir(folder)
    numbers = []
    for file in files:
        match = re.search(r"density_t.(\d+)\.h5", file)
        if match:
            numbers.append(int(match.group(1)))
    return np.sort(np.array(numbers))

def read_field_fromfile(name, field):
    file = h5.File(name, "r")
    f = np.zeros((Nx,Ny,Nz))
    name_field="/LBE3D/"+field
    f[:,:,:] = file.get(str(name_field))[...]
    file.close()
    return f

def swap_XZ(Nx,Nz):
    return Nx,Nz;

# constants of this simulation
viscosity_kin = 1/6
sigma = 0.0235
rho_tot = 1.36

num_run = str(sys.argv[1])
direc = "/data/storage19/mattia/DATA/stag_"+num_run+"_run_mayo/"

name=direc+"RUN/density_t.0.h5"
Nx,Ny,Nz,dt = read_dim_fromfile(name)
numbs = find_outs(direc+"RUN/")
numbs = numbs[1:]  # remove step 0
#numbs = numbs[-1:]  # only last step

#Nz,Nx = swap_XZ(Nx,Nz);

R_ava = read_Rava(direc+"droplets_volume_rho1.dat")

rho1 = np.zeros((Nx,Ny,Nz)); rho1[:,:,:] = np.nan;
#rho2 = np.zeros((Nx,Ny,Nz)); rho2[:,:,:] = np.nan;
vx   = np.zeros((Nx,Ny,Nz)); vx[:,:,:] = np.nan;
vy   = np.zeros((Nx,Ny,Nz)); vy[:,:,:] = np.nan;
vz   = np.zeros((Nx,Ny,Nz)); vz[:,:,:] = np.nan;
#v2   = np.zeros((Nx,Ny,Nz)); v2[:,:,:] = np.nan;
rhoA = np.nan;
vrms = np.nan;
Re_N = np.nan;
We_N = np.nan;
Ca_N = np.nan;

file = open("dimensionless_numbers"+num_run+".csv", "w")
file.write("STEP, Reynolds, Weber, Capillarity, Vrms, Rava\n")
file.flush()

for t in numbs:
    print("t = %d"%t,file=sys.stderr)
    #name=direc+"density_t.%d.h5"%t
    #rho1[:,:,:] = read_field_fromfile(name,"rho1")
    #rho2[:,:,:] = read_field_fromfile(name,"rho2")
    name=direc+"RUN/velocity_t.%d.h5"%t
    vx[:,:,:] = read_field_fromfile(name,"vx")
    vy[:,:,:] = read_field_fromfile(name,"vy")
    vz[:,:,:] = read_field_fromfile(name,"vz")
    print("   v read\n",file=sys.stderr)
    ###
    vrms = np.sqrt(np.sum(vx[:,:,:]*vx[:,:,:]+vy[:,:,:]*vy[:,:,:]+vz[:,:,:]*vz[:,:,:])/(Nx*Ny*Nz))
    #rhoA = np.sum(rho1[:,:,:])/(Nx*Ny*Nz)
    ###
    R_idx = int(t/dt)-1
    Re_N = vrms*Nz/viscosity_kin
    We_N = rho_tot*vrms*vrms*R_ava[R_idx,1]/sigma
    Ca_N = viscosity_kin*vrms*rho_tot/sigma
    ###
    file.write(f"{t:d}, {Re_N:1.20g}, {We_N:1.20g}, {Ca_N:1.20g}, {vrms:1.20g}, {R_ava[R_idx,1]:1.20g}\n")
    file.flush()

file.close()
