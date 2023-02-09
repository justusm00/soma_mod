import h5py
import numpy as np
import matplotlib.pyplot as plt

#plot mean square displacement over time

#load ana data
with h5py.File(f'coord_ana.h5', 'r') as f:
    MSD=np.array(f['MSD'])
    delta_mc_MSD=np.array([f['MSD'].attrs["DeltaMC"]])
   
MSDA=MSD[:,3]
t_end=len(MSDA)*delta_mc_MSD
t=np.arange(0,t_end,delta_mc_MSD)
plt.plot(t,MSDA)
plt.show()
