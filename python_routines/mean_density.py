import h5py
import numpy as np
import matplotlib.pyplot as plt

#routine to plot mean density field at given x


boxlen=3 #length of box/conversion zones

idx=1

#get density field
with h5py.File('coord_ana.h5','r') as anafile:
    phi=np.array(anafile['density_field']).mean(axis=(3,4))
    delta_mc_phi=np.array(anafile['density_field'].attrs["DeltaMC"])
num_steps=len(phi) #number of MC sweeps
tmax=num_steps*delta_mc_phi
t=np.arange(0,tmax,delta_mc_phi)
phiA=phi[:,0]
phiB=phi[:,1]
phitot=phiA+phiB
phiA=phiA/phitot
phiB=phiB/phitot
print(f"Mean: {np.mean(phiA[:,idx])}")
print(f"Variance: {np.std(phiA[:,idx])}")
plt.plot(t,phiA[:,idx])
plt.show()
