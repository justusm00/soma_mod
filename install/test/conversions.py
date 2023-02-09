import h5py
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


##load analysis data
with h5py.File(f'coord_ana.h5', 'r') as f:
    ##numbers of polymers converted per reaction
    num_conversions=np.array([f['num_conversions']]) 
    delta_mc_num_conversions=np.array([f['num_conversions'].attrs["DeltaMC"]])

num_conversions=num_conversions[0,:,2]
tmax=delta_mc_num_conversions*len(num_conversions)

dummy=deepcopy(num_conversions)
print(dummy)
for i in range(1,len(num_conversions)):
    num_conversions[i]=num_conversions[i]-dummy[i-1]

t=np.arange(0,tmax,delta_mc_num_conversions)
plt.plot(t,num_conversions)
plt.show()
