import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit

##routine to extract mean density fields and Onsager coefficient

##load parameters from coord
with h5py.File(f'coord.h5', 'r') as f:
    ##number of beads per polymer
    N=int(f['parameter/reference_Nbeads'][()]) 
    ##number of polymers
    n_polym=int(f['parameter/n_polymers'][()]) 
    ##box dimensions 
    lxyz=np.array([f['parameter/lxyz']]).flatten() 
    lx=lxyz[0]
    ly=lxyz[1]
    lz=lxyz[2]
    ##box cross section
    A=ly*lz 
    ##box volume
    V=lx*ly*lz
    ##average bead density
    rho0=n_polym*N/V
    #box discretization
    nxyz=np.array([f['parameter/nxyz']]).flatten() 
    nx=nxyz[0]
    ny=nxyz[1]
    nz=nxyz[2]
    

##load analysis data
with h5py.File(f'coord_ana.h5', 'r') as f:
    ##density field
    phi=np.array(f['density_field']).mean(axis=(0,3,4))
    ##numbers of polymers converted per reaction
    num_conversions=np.array([f['num_conversions']]) 
    delta_mc_num_conversions=np.array([f['num_conversions'].attrs["DeltaMC"]])
    ##average squared end to end distance
    Re=np.array(f['Re']) 
    delta_mc_Re=np.array([f['Re'].attrs["DeltaMC"]])
    ##MSD
    MSD=np.array(f['MSD'])
    

ReA=np.mean(Re[:,0])
MSDA=MSD[:,3]
MSDB=MSD[:,11]
##Nbar
Nbar=n_polym/V
##simulation duration in MC steps
t=delta_mc_num_conversions*len(num_conversions[0,:]) 
##A current
J=-num_conversions[0,-1,2]/(t*A)
##diffusion constant of A
DA=MSDA[-1]/(6*t)
##diffusion constant of B
DB=MSDB[-1]/(6*t)
##rouse time
tr=ReA/DA 

## Fitting function for density profile
##
## @param m Slope
## @return fit function
def f(x,m,b):
    return m*x+b


#### MAIN #####

##length of box/conversion zones
boxlen=3 
##normalize densities
phiA=phi[0,boxlen:nx-boxlen]
phiB=phi[1,boxlen:nx-boxlen]
phitot=phiA+phiB
phiA=phiA/phitot
phiB=phiB/phitot

##fit density to get gradient
xmin=boxlen/nx*lx
xmax=(nx-boxlen)/nx*lx
x=np.linspace(xmin,xmax,1000)
x_grid=np.linspace(xmin,xmax,nx-2*boxlen)
popt, pcov = curve_fit(f, x_grid, phiA)
phi_fit=f(x,popt[0],popt[1])
gradphi=popt[0] #density gradient


##gradient of chemical potential/kb*T
gradmu=gradphi/N*(1/phiA+1/(1-phiA))*Nbar
##lambda in units of 1/(tr*R_e)
lam=-J/gradmu


##plot
fig,ax=plt.subplots(2,2,figsize=(8,8))
##plot phi_A
ax[0,0].set_title("density $\phi_A$")
ax[0,0].set_ylabel("$\phi_A/(\phi_A+\phi_B)$")
ax[0,0].grid()
#ax[0].plot(x,f(x,popt))
ax[0,0].plot(x_grid,phiA,color='b')

##plot grad_mu
ax[0,1].set_title(r"$\nabla\mu$ ")
ax[0,1].set_ylabel(r"$\frac{\nabla\mu}{k_BT}$   $[\frac{1}{R_e}]$")
ax[0,1].grid()
ax[0,1].plot(x_grid,gradmu,color='b')

##plot lambda
ax[1,0].set_title("$\Lambda$")
ax[1,0].set_xlabel("$x$ $[R_e]$")
ax[1,0].set_ylabel(r"$\Lambda$   $[\frac{1}{MCstep\cdot R_e}]$")
ax[1,0].grid()
ax[1,0].plot(x_grid,lam,color='b')

##plot lambda/(N*phiA*phiB)
ax[1,1].set_title(r"$D_{A}=\frac{\Lambda}{N\phi_A \phi_B}$")
ax[1,1].set_xlabel("$x$ $[R_e]$")
ax[1,1].set_ylabel(r"$D_{A}$ $[\frac{1}{MCstep\cdot R_e}]$")
ax[1,1].grid()
ax[1,1].plot(x_grid,lam/(N*phiA*phiB),color='b')

plt.tight_layout()
plt.show() 


