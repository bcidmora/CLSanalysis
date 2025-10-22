import h5py
import numpy as np
import os
import set_of_functions as vf

### This is the source directory where the codes, data, etc. 
location = os.path.expanduser('~')+'$YOUR_PATH_TO_THE_DATA$'


hdf5NameMulti = 'cls21_X451_r001_isodoublet_Sm2_fwd_135_cnfgs.hdf5' 
hdf5NameSingle= 'cls21_X451_r001_single_hadrons_135_cnfgs.hdf5'

### ENSEMBLE FILES: X451
f = h5py.File(location+'data/X451/'+hdf5NameMulti,'r') 
f1 = h5py.File(location+'data/X451/'+hdf5NameSingle,'r') 

### Reweighting factors ###NEW
weight_raw = np.loadtxt(location + 'data/X451/X451r001.ms1.dat_ascii', unpack=True)[1]
weight_raw = np.array(weight_raw)


### How Many Irreps has each?
name = list(f.keys())
name1 = list(f1.keys())

### Name of the quantum number to study (Choose any name that's convenient to remember)
if 'isosinglet' in hdf5NameMulti: 
    the_hadron_state = '_isosinglet_strange_fermionic_'
    the_non_interacting_levels = [ ] # P3_F2
elif 'isodoublet' in hdf5NameMulti: 
    the_hadron_state = '_isodoublet_strange_fermionic_'
    the_non_interacting_levels = [ ['K(1)L(1)', 'P(1)X(1)'], # P0_G1g 
                                  [ 'K(0)S(0)', 'K(0)L(0)', 'P(0)X(0)' , 'K(1)L(1)', 'P(1)X(1)'], # P0_G1u
                                  ['K(1)L(1)', 'P(1)X(1)'], # P0_Hg
                                  ['P(0)X(0)', 'K(1)L(1)','P(1)X(1)'], # P0_Hu
                                  ['K(0)S(1)', 'K(0)L(1)', 'P(0)X(1)', 'K(1)L(0)', 'P(1)X(0)'], # P1_G1
                                  ['P(0)X(1)'], # P1_G2
                                  ['K(0)S(2)', 'K(0)L(2)', 'P(0)S(2)', 'K(1)L(1)', 'P(1)X(1)'], # P2_G
                                  ['K(0)S(3)', 'K(0)L(3)', 'P(0)X(3)'], # P3_G
                                  ['K(1)L(2)', 'P(1)S(2)'], # P3_F1
                                  ['K(1)L(2)', 'P(1)S(2)'], # P3_F2
                                  ['K(0)S(4)', 'K(0)L(4)', 'P(0)X(4)', 'K(1)S(1)', 'K(1)L(1)','P(1)X(1)']] # P4_G1 # The last ones are summed mom
else: 
    the_hadron_state='_' #You can put here anything you want
    the_non_interacting_levels = []

### Number of gauge configurations
### This is normalizing it to the total amount of gauge configs originally in this ensemble.
ncfgs = 2026

### Reweighting factors
pre_weight = np.array(vf.RW_NORMALIZATION(weight_raw,ncfgs), dtype=np.float128)

chosen_operators_list = [[], [], [], []]

### Number of gauge configurations
ncfgs = np.array(f1[name1[0]+'/data']).shape[0] 
cnfgs_list0 = list(np.arange(7,393,16)) 
cnfgs_list1 = list(np.arange(175,895,16)) 
cnfgs_list2 = list(np.arange(911,1903,16))
cnfgs_list3 = list(np.arange(1919,1952,16))
raw_cnfgs_list = (cnfgs_list0 + cnfgs_list1 + cnfgs_list2 + cnfgs_list3)
sorted_cnfgs_list = sorted(raw_cnfgs_list)
cnfgs_list = np.array(sorted_cnfgs_list)
    
# Now this is normalizing with respect to the configs that are effectively used. 
worked_weight = []
for ii in range(len(cnfgs_list)):
    worked_weight.append(pre_weight[cnfgs_list[ii]])    

### Final reweighting factors
weight = np.array(vf.RW_NORMALIZATION(worked_weight,ncfgs), dtype=np.float128)

### List of tmaxs used for the fitting procedure. 
listTMaxSingleHads = [24,24,23,21,23,25,24,25,24,21,25,23,25,24,24,22,22,24,22,23,22,21,22,20,20,21,19,20]

### TMax used for the fits of correlation matrices
listTMaxMultiHads=[]
for ix in range(0,len(name)):
    list_tmax_multihads_ix=[]
    for jx in range(0,np.array(f[name[ix]+'/data']).shape[1]):
        list_tmax_multihads_ix.append(np.array(f[name[ix]+'/data']).shape[-1])
    listTMaxMultiHads.append(list_tmax_multihads_ix)

### Minimum time slices used for the fits of single hadron correlators.
singleTMinsFitPlots = [11,10,15,14,15,15, # P^{2} = 0
                       11,12,12,12,12,11, # P^{2} = 1
                       12,12,11,12,12,12, # P^{2} = 2
                       12,12,11,12,12,11, # P^{2} = 3
                       12,12,11,11] # P^{2} = 4

### Minimum time slices for the fits of multihadron correlators
multiTMinsFitPlots = []
for ix in range(0,len(name)):
    multiTMinsFitPlots_ij=[]
    for jx in range(0,np.array(f[name[ix]+'/data']).shape[1]):
        multiTMinsFitPlots_ij.append(10)
    multiTMinsFitPlots.append(multiTMinsFitPlots_ij)

### PARAMETERS OF THE LATTICE
aLat = np.float64(0.0761) # in fm
betaLat = np.float64(3.46)
fmToMev = np.float64(197.327) #1 fm = 197.3 Mev
LatSize = np.float64(40.)
