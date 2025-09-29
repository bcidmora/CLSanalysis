import h5py
import numpy as np
import os
import set_of_functions as vf

### This is the source directory where the codes, data, etc. 
location = os.path.expanduser('~')+'$YOUR_PATH_TO_THE_DATA$'


hdf5NameMulti = 'cls21_X451_r001_isodoublet_Sm2_fwd_110_cnfgs.hdf5' 
hdf5NameSingle= 'cls21_X451_r001_single_hadrons_110_cnfgs.hdf5'

### ENSEMBLE FILES: X451
f = h5py.File(location+'data/X451/'+hdf5NameMulti,'r') 
f1 = h5py.File(location+'data/X451/'+hdf5NameSingle,'r') 

### Reweighting factors ###NEW
weight_raw = np.loadtxt(location + 'data/X451/X451r001.ms1.dat_ascii', unpack=True)[1]
weight_raw = np.array(weight_raw)


### How Many Irreps has each?
name = list(f.keys())
name1 = list(f1.keys())


### Number of gauge configurations
### This is normalizing it to the total amount of gauge configs originally in this ensemble.
ncfgs = 2026

### Reweighting factors
pre_weight = np.array(vf.RW_NORMALIZATION(weight_raw,ncfgs), dtype=np.float128)



chosen_operators_list = [[], [], [], []]


### Number of gauge configurations
ncfgs = np.array(f1[name1[0]+'/data']).shape[0] # SH
cnfgs_list1 = list(np.arange(175,894,16)) 
cnfgs_list2 = list(np.arange(911,1902,16))
cnfgs_list3 = list(np.arange(1919,1952,16))
cnfgs_list = np.array(cnfgs_list1 + cnfgs_list2 + cnfgs_list3)
    
# Now this is normalizing with respect to the configs that are effectively used. 
worked_weight = []
for ii in range(len(cnfgs_list)):
    worked_weight.append(pre_weight[cnfgs_list[ii]])    

### Final reweighting factors
weight = np.array(vf.RW_NORMALIZATION(worked_weight,ncfgs), dtype=np.float128)

### List of tmaxs used for the fitting procedure. 
listTMaxSingleHads = [np.array(f1[name1[0]+'/data']).shape[-1]]*len(name1)

### TMax used for the fits of correlation matrices
listTMaxMultiHads=[]
for ix in range(0,len(name)):
    list_tmax_multihads_ix=[]
    for jx in range(0,np.array(f[name[ix]+'/data']).shape[1]):
        list_tmax_multihads_ix.append(np.array(f[name[ix]+'/data']).shape[-1])
    listTMaxMultiHads.append(list_tmax_multihads_ix)

### Minimum time slices used for the fits of single hadron correlators.
singleTMinsFitPlots = [10]*len(name1)

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
