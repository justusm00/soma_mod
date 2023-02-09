# -*- coding: utf-8 -*-
# @Author: Justus Multhaup
# @Date:   2022-12-15 11:21:50
# @Last Modified by:   Your name
# @Last Modified time: 2023-01-21 19:30:03
import h5py
import numpy as np
import matplotlib.pyplot as plt

#routine to plot mean density field at given x without area51



with h5py.File('coord.h5', 'r') as f:
    ##number of beads per polymer
    N=int(f['parameter/reference_Nbeads'][()]) 
    ##number of polymers
    n_polym=int(f['parameter/n_polymers'][()]) 
    n_poly_type=int(f['parameter/n_poly_type'][()]) 
    ##box dimensions 
    lxyz=np.array(f['parameter/lxyz'])
    ##box discretization
    nxyz=np.array(f['parameter/nxyz'])
    ##umbrella field (target density)
    phi_target=np.array(f['umbrella_field'])



##load analysis data
with h5py.File(f'coord_ana.h5', 'r') as f:
    ##density field averaged over time, y and z
    phi=np.array(f['density_field'])
    delta_mc_phi=np.array(f['density_field'].attrs["DeltaMC"])
    ##sum of density fields
    phitot=sum(phi)
    ##normalized densities
    phi_norm=phi/phitot
    #phi_norm=phi
    ##numbers of polymers converted per reaction
    num_conversions=np.array(f['num_conversions']) 
    delta_mc_num_conversions=np.array(f['num_conversions'].attrs["DeltaMC"])
    ##MSD
    MSD=np.array(f['MSD'])


scale=np.prod(nxyz)/(N*n_polym)
phi=(phi*scale)
phiA=phi[:,0]
phiB=phi[:,1]
phi_std=np.std(phi,axis=0)
cum_std=sum(phi_std[phi_target > 0])
print(cum_std)
