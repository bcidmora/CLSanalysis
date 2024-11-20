import h5py
import numpy as np
import os
import random
import set_of_functions as vf
        
location = os.path.expanduser('~') + '$YOUR_PATH_TO_THE_DATA$'

# ENSEMBLE FILES: D200
f = h5py.File(location + '/D200/cls21_D200_r000_isosinglet_Sm1_fwd.hdf5','r')
f1 = h5py.File(location + '/D200/cls21_D200_r000_single_fwd.hdf5','r')

weight_sarah = np.loadtxt(location + '/D200/sarah_reweighting.dat')

name = list(f.keys())
name1 = list(f1.keys())

### Number of gauge configurations
ncfgs = np.array(f[name[0]+'/data']).shape[0]

weight = np.array(vf.RW_NORMALIZATION(weight_sarah, ncfgs), dtype=np.float64)

## List of tmaxs used for the fitting procedure. 
#listTMaxSingleHads = [np.array(f1[name1[0]+'/data']).shape[-1]+1]*len(name1)
listTMaxSingleHads = [25,25,25,25,25,25,25,25,25,25, 25,25,25,25,25,25,25,25,25,25, 25,25,25,25,25,25,25,25,25]

# For the multihadron operators, I modified each of the entries. 
#list_tmax_multihads=[]
#for ix in range(0,len(name)):
    #list_tmax_multihads_ix=[]
    #for jx in range(0,np.array(f[name[ix]+'/data']).shape[1]):
        #list_tmax_multihads_ix.append(np.array(f[name[ix]+'/data']).shape[-1]+1)
    #list_tmax_multihads.append(list_tmax_multihads_ix)
    
listTMaxMultiHads = [[25,22,23,19,17],
                       [22,24,24,23,23,22,20,17],
                       [24,25,23],
                       [24,24,23,24,23,24,22,24,22],
                       [24,23,23,22,23,21,20],
                       [23,24,23,22,23,24,24,23,22,23,23,20,19,20,18],
                       [24,22,23,23,23,20],
                       [24,23,24,22,22,20],
                       [23,22,24,24,24,23,24,23,23,23,23,19,19]]


singleTMinsFitPlots = [8,10,10,10,10,10,8,10,10,10,10,10,10,10,10,8,8,8,6,8,8,8,8,8,6,8,8,8,10]
multiTMinsFitPlots = [[14, 13, 13, 11, 8],  #PSQ0 G1g
                    [14, 14, 14, 13, 11, 11, 11, 8], #PSQ0 G1u
                    [14, 14, 14], #PSQ0 Hu
                    [15, 15, 14, 14, 14, 13, 12, 12, 13], #PSQ1 G1
                    [14, 15, 14, 14, 12, 12, 10], #PSQ1 G2
                    [13, 14, 13, 13, 14, 14, 13, 13, 11, 12, 13, 10, 11, 12, 10], #PSQ2 G
                    [13, 13, 12, 12, 12, 12], #PSQ3 F1
                    [13, 13, 12, 13, 13, 10], #PSQ3 F2
                    [14, 13, 13, 12, 13, 12, 12, 13, 12, 13, 12, 11, 10]] #PSQ3 G 

# PARAMETERS OF THE LATTICE
aLat = np.float64(0.0633) # in fm
betaLat = np.float64(3.55)# 
fmToMev = np.float64(197.327) #1 fm = 197.3 Mev
LatSize = np.float64(64.)
