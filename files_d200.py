import h5py
import numpy as np
import os
import random
import set_of_functions as vf

### This is the source directory where the codes, data, etc.         
location = os.path.expanduser('~') + '$YOUR_PATH_TO_THE_DATA$'

### Names of the files with single correlators and multi hadrons
hdf5NameMulti = 'cls21_D200_r000_isosinglet_Sm1_fwd.hdf5'
hdf5NameSingle = 'cls21_D200_r000_single_fwd.hdf5'

### Reweighting factors
weight_sarah = np.loadtxt(location + '/data/D200/sarah_reweighting.dat')

### ENSEMBLE FILES: D200
f = h5py.File(location + 'data/D200/' + hdf5NameMulti, 'r')
f1 = h5py.File(location + 'data/D200/' + hdf5NameSingle, 'r')

### How Many Irreps has each?
name = list(f.keys())
name1 = list(f1.keys())

### Name of the quantum number to study
the_hadron_state='_'

if 'isosinglet' in hdf5NameMulti: 
    the_hadron_state = '_isosinglet_strange_fermionic_'
    the_non_interacting_levels = [] # Here goes a list of the possible nearby non-interacting levels (Look at the X451 for an example)
elif 'isodoublet' in hdf5NameMulti: 
    the_hadron_state = '_isodoublet_strange_fermionic_'
    the_non_interacting_levels = [] # Here goes a list of the possible nearby non-interacting levels (Look at the X451 for an example)
else: 
    the_hadron_state='_' #You can put here anything you want
    the_non_interacting_levels = []

### Number of gauge configurations
# ncfgs = np.array(f[name[0]+'/data']).shape[0]
ncfgs = 2000 

### Final reweighting factors
weight = np.array(vf.RW_NORMALIZATION(weight_sarah, ncfgs), dtype=np.float64)

### This is the list of chosen operators 
chosen_operators_list = []

### This part here is only in case you have a subset of the original data, for my case it was 100 cnfgs, every 20 starting on 20 until 2000.
if '-' in hdf5NameMulti:
    ncfgs = np.array(f[name[0]+'/data']).shape[0]
    pre_weight = weight
    cnfgs_list_string = vf.GETTING_MIN_MAX_CONFIGS(hdf5NameSingle) 
    cnfgs_list = np.arange(int(cnfgs_list_string[0])-1, int(cnfgs_list_string[1])+1, int(cnfgs_list_string[2]))
    worked_weight = []
    for ii in range(len(cnfgs_list)):
        worked_weight.append(pre_weight[cnfgs_list[ii]])   
    weight = np.array(vf.RW_NORMALIZATION(worked_weight,ncfgs), dtype=np.float64)

### List of tmaxs used for the fitting procedure. 
listTMaxSingleHads = [25,25,25,25,25,25,25,25,25,25, 25,25,25,25,25,25,25,25,25,25, 25,25,25,25,25,25,25,25,25]

### TMax used for the fits of correlation matrices    
listTMaxMultiHads = [[25,22,23,19,17],
                       [22,24,24,23,23,22,20,17],
                       [24,25,23],
                       [24,24,23,24,23,24,22,24,22],
                       [24,23,23,22,23,21,20],
                       [23,24,23,22,23,24,24,23,22,23,23,20,19,20,18],
                       [24,22,23,23,23,20],
                       [24,23,24,22,22,20],
                       [23,22,24,24,24,23,24,23,23,23,23,19,19]]

### Minimum time slices used for the fits of single hadron correlators.
singleTMinsFitPlots = [8,10,10,10,10,10, #P^{2} = 0
                       8,10,10,10,10,10, # P^{2} = 1
                       8,10,8,8,10,8, # P^{2} = 2
                       8,10,10,8,10,8, # P^{2} = 3
                       8,8,8,8,10] # P^{2} = 4


### Minimum time slices for the fits of multihadron correlators
multiTMinsFitPlots = [[14, 13, 13, 11, 8],  #PSQ0 G1g
                    [14, 14, 14, 13, 11, 11, 11, 8], #PSQ0 G1u
                    [14, 14, 14], #PSQ0 Hu
                    [15, 15, 14, 14, 14, 13, 12, 12, 13], #PSQ1 G1
                    [14, 15, 14, 14, 12, 12, 10], #PSQ1 G2
                    [13, 14, 13, 13, 14, 14, 13, 13, 11, 12, 13, 10, 11, 12, 10], #PSQ2 G
                    [13, 13, 12, 12, 12, 12], #PSQ3 F1
                    [13, 13, 12, 13, 13, 10], #PSQ3 F2
                    [14, 13, 13, 12, 13, 12, 12, 13, 12, 13, 12, 11, 10]] #PSQ3 G 

### PARAMETERS OF THE LATTICE
aLat = np.float64(0.0633) # in fm
betaLat = np.float64(3.55)# 
fmToMev = np.float64(197.327) #1 fm = 197.3 Mev
LatSize = np.float64(64.)
