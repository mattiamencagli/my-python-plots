import numpy as np

''' constant parameters: '''
min_i = 1e-8
rad_to_deg = 180.0 / np.pi


################################################################################################################################
################################################################################################################################
def open_reg(direc, name, jump=1, checkstop=False):
    ''' open the files in which I write the data of the refularized systems '''
    f = open(direc + name)
    file=f.readlines()
    f.close()
    
    #masses of the bodies
    m=[]
    while(len(file[0])<20):
        m.append(float(file[0]))
        file.pop(0)
    
    #check how many time the reg sys has reached the synch with its center of mass
    stop=[]
    if(checkstop):
        i=0
        print("stops : ", end = '')
        while i<len(file):
            if(file[i]=="stop\n"):
                stop.append(i)
                file.pop(i)
                print("%d, "%i, end = '')
            i += 1
        if(len(stop)!=0):
            stop.pop()
        print()
    
    #the data
    o = np.loadtxt(file)
    print(o.shape)
    
    #remove some points to decrease the workload
    if(jump>1):
        keep = []
        for i in range(len(o)):
            if(i%jump == 0):
                keep.append(o[i])
        o=keep
    o = np.array(o)
    m = np.array(m)
    print(o.shape)
    
    return o, stop, m, len(m)


################################################################################################################################
################################################################################################################################
def change_units(o, m, m_unit, t_unit, v_unit, conv_pc_UA):
    ''' change the unit of measure '''
    o[:,0:3*len(m)] *= conv_pc_UA
    o[:,3*len(m):6*len(m)] *= v_unit
    o[:,6*len(m)] *= t_unit
    m *= m_unit
    print("units: M_sol, year, km/s, UA")
    return


################################################################################################################################
################################################################################################################################
def kepl_units(a,T,inc,ome,Ome,nu,t_unit,conv_pc_UA):
    ''' change the unit of measure '''
    a *= conv_pc_UA
    T *= t_unit
    inc *= rad_to_deg
    ome *= rad_to_deg
    Ome *= rad_to_deg
    nu *= rad_to_deg
    print("units: UA, years, deg")
    return


################################################################################################################################
################################################################################################################################
def acos2(num, den, posflag):
    ''' do the arccos**2 as in tsunami (FORSE NP.ARCCOS BASTA? CONTROLLA....) '''
    cose = num/den
    arccos = np.zeros_like(cose)
    
    mask1 = (cose > -1.0) & (cose < 1.0)
    arccos[mask1] = np.arccos(cose[mask1])
    
    mask2 = mask1 & (posflag < 0.0)
    arccos[mask2] *= -1.0
    
    mask3 = cose <= -1.0
    arccos[mask3] = np.pi
    
    #quindi per cose>1 si ha arccos = 0
    return arccos


def cart_to_kepl(m1, pos1, vel1, m2, pos2, vel2):
    ''' pass from cartesian coordinates to keplerian coordinates, it needs unit of meaure with G=1 (so the isteddas units are ok) '''
    
    dr = pos1 - pos2
    dv = vel1 - vel2

    r = np.sqrt( dr[:,0]*dr[:,0] + dr[:,1]*dr[:,1] + dr[:,2]*dr[:,2] )
    v2 = dv[:,0]*dv[:,0] + dv[:,1]*dv[:,1] + dv[:,2]*dv[:,2]

    mu = m1 + m2
    vcirc2 = mu/r

    #semimajor axis
    a = -mu / ( v2 - 2.0*vcirc2 )
    
    #period
    T = 2 * np.pi * np.sqrt(np.abs(a*a*a)/(m1+m2))

    #angular momentum vector
    lx = dr[:,1]*dv[:,2] - dr[:,2]*dv[:,1]
    ly = dr[:,2]*dv[:,0] - dr[:,0]*dv[:,2]
    lz = dr[:,0]*dv[:,1] - dr[:,1]*dv[:,0]
    l = np.sqrt( lx*lx + ly*ly + lz*lz )

    vdiff2 = v2 - vcirc2
    rvr = (dr[:,0]*dv[:,0] + dr[:,1]*dv[:,1] + dr[:,2]*dv[:,2])
    vr = rvr/r
    muinv = 1.0/mu

    #eccentricity
    ex = muinv * (vdiff2*dr[:,0] - rvr*dv[:,0])
    ey = muinv * (vdiff2*dr[:,1] - rvr*dv[:,1])
    ez = muinv * (vdiff2*dr[:,2] - rvr*dv[:,2])
    e = np.sqrt( ex*ex + ey*ey + ez*ez )
#     e = np.sqrt(1 + (v2-2.0*vcirc2)*l*l/(mu*mu))

    #inclination
    inc = acos2(lz, l, np.ones_like(l))

    #line of the nodes
    nx = -ly
    ny =  lx
    n = np.sqrt( nx*nx + ny*ny )

    #Longitude of ascending node
    Ome = acos2(nx, n, ny)

    true_long = np.zeros_like(l)
    peri_long = np.zeros_like(l)
    ome = np.zeros_like(l)
    nu = np.zeros_like(l)
    ome_nu = np.zeros_like(l)
    
    mask1 = (inc < min_i) | (inc > np.pi - min_i)
    no_mask1 = ~mask1
    mask2 = (inc < np.pi * 0.5)
    m_comb1 = mask1 & mask2
    m_comb2 = mask1 & (~mask2)    
#     m_comb3 = no_mask1 & mask2
#     m_comb4 = no_mask1 & (~mask2)
    
    #planar case
    true_long[mask1] = acos2(dr[mask1,0], r[mask1], dr[mask1,1]) #true longitude (planar)
    peri_long[mask1] = acos2(ex[mask1], e[mask1], ey[mask1]) #longitude of pericenter (planar)
    
    ome[m_comb1] = peri_long[m_comb1] - Ome[m_comb1] #argument of pericenter (planar, prograde)
    nu[m_comb1] = true_long[m_comb1] - peri_long[m_comb1] #true anomaly (planar, prograde)
    
    ome[m_comb2] = Ome[m_comb2] - peri_long[m_comb2] #argument of pericenter (planar, retrograde)
    nu[m_comb2] = peri_long[m_comb2] - true_long[m_comb2] #true anomaly (planar, retrograde)
            
    #non-planar case
    ome_nu[no_mask1] = acos2(nx[no_mask1]*dr[no_mask1,0] + ny[no_mask1]*dr[no_mask1,1], n[no_mask1]*r[no_mask1], dr[no_mask1,2])
    ome[no_mask1] = acos2(nx[no_mask1]*ex[no_mask1] + ny[no_mask1]*ey[no_mask1], n[no_mask1]*e[no_mask1], ez[no_mask1])
    
    nu[no_mask1] = ome_nu[no_mask1] - ome[no_mask1] #true anomaly
    #qui trani fa operazioni extra, forse gli servono per i termini post newtoniani
#     nu[m_comb3] = ome_nu[m_comb3] - ome[m_comb3] #true anomaly (inclined, prograde)
#     nu[m_comb4] = ome_nu[m_comb4] - ome[m_comb4] #true anomaly (inclined, retrograde)

    return a, T, e, inc, ome, Ome, nu


################################################################################################################################
################################################################################################################################
def comp_t_gw_in_Gyr(a, m1, m2):
    '''calcola tempo scala coalescenza BH, units: a[UA], m[M_sol]'''

    G = 4.30091e-3 * 206264.8062471 #AU (km/s)^2 M_sol^-1
    c = 299792.458 #km/s
    
    # (gli ultimi 2) il primo fattore e' per passare da sec a Gyr, il secondo e' per correggere un AU/km rimanente
    return ( ( 5.0 * c**5 * a**4 / (265.0 * G**3 * m1 * m2 * (m1 + m2)) ) * (3.171e-17 * 1.496e+8) )


################################################################################################################################
################################################################################################################################
