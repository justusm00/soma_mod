# -*- coding: utf-8 -*-
# @Author: Your name
# @Date:   2022-12-09 15:54:44
# @Last Modified by:   Your name
# @Last Modified time: 2023-01-21 21:34:07
import h5py
import numpy as np

phi_left_A=0.3
phi_right_A=0.7
phi_left_B=0.7
phi_right_B=0.3
phi_top_A=0.3
phi_top_B=0.7
#phi_bot_A=np.array([0.1,0.1,0.9,0.9,0.1,0.1,0.9,0.9,0.1,0.1])
phi_bot_A=np.array([0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9])
phi_top_A=phi_bot_A
#phi_bot_B=np.array([0.9,0.9,0.1,0.1,0.9,0.9,0.1,0.1,0.9,0.9])
phi_bot_B=np.array([0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9,0.1])
phi_top_B=phi_bot_A


##change coord
f= h5py.File(f'coord.h5', 'r+') 
nxyz=np.array(f['parameter/nxyz'])
nz=nxyz[2]
##lamellar floor
for z in range(nz):
    f['umbrella_field'][0,:,0,z] = phi_bot_A[z]
    f['umbrella_field'][1,:,0,z] = phi_bot_B[z]
    f['umbrella_field'][0,:,-1,z] = phi_bot_A[z]
    f['umbrella_field'][1,:,-1,z] = phi_bot_B[z]

""" f['umbrella_field'][0,0,:,:]=phi_left_A
f['umbrella_field'][0,-1,:,:]=phi_right_A
f['umbrella_field'][1,0,:,:]=phi_left_B
f['umbrella_field'][1,-1,:,:]=phi_right_B
f['umbrella_field'][0,:,-1,:]=phi_top_A
f['umbrella_field'][1,:,-1,:]=phi_top_B
f['umbrella_field'][0,:,0,:]=phi_bot_A
f['umbrella_field'][1,:,0,:]=phi_bot_B """
f.close()

"""
##change coord_ana
f= h5py.File(f'coord_ana.h5', 'r+') 
f['umbrella_field'][0,0,:,:]=phi_left
f['umbrella_field'][0,-1,:,:]=phi_right
f.close()
"""





