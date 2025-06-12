import h5py
import numpy as np
import os
import set_of_functions as vf

### This is the source directory where the codes, data, etc. 
location = os.path.expanduser('~')+'$YOUR_PATH_TO_THE_DATA$'

### lambdas_analysis: If True: It does the analysis over the lambda correlator matrices, meaning it will choose the normalization of the reweighting factors according to the gauge configurations that were actually computed. 
### lambdas_analysis: If False: It does the analysis over the single hadrons only, which were the first 20 gauge configs for the pion and 100 gauge configs 
### nucleon_analysis: If True, it only analyses the new nucleon (100 cnfgs) averaged over irrep row and over source time.

lambdas_analysis = False
nucleon_analysis = True

hdf5NameMulti = 'cls21_X451_r001_isosinglet_Sm1_fwd_t01.hdf5'


if nucleon_analysis: 
    hdf5NameSingle = 'cls21_X451_r001_singles_nucleon.hdf5'
else:
    hdf5NameSingle = 'cls21_X451_r001_singles.hdf5'

### ENSEMBLE FILES: N451
# f = h5py.File(location+'data/X451/cls21_X451_r001_isosinglet_Sm1_fwd.hdf5','r') #lambdas
# f = h5py.File(location+'data/X451/cls21_X451_r001_isosinglet_Sm1_fwd_t00.hdf5','r') #lambdas T00
# f = h5py.File(location+'data/X451/cls21_X451_r001_isosinglet_Sm1_fwd_t01.hdf5','r') #lambdas T01
# f1 = h5py.File(location+'data/X451/cls21_X451_r001_singles.hdf5','r') #SH
f = h5py.File(location+'data/X451/'+hdf5NameMulti,'r') #lambdas T01
f1 = h5py.File(location+'data/X451/'+hdf5NameSingle,'r') #SH

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

### Number of gauge configurations
if lambdas_analysis:
    ncfgs = np.array(f[name[0]+'/data']).shape[0] # lambdas
    # cnfgs_list = np.arange(39,2010,40) # lambdas T00 
    # cnfgs_list = np.arange(19,2010,40) # lambdas T01 
    cnfgs_list = np.arange(19,2010,20) # ALL t00 and t01 
else:
    ncfgs = np.array(f1[name1[0]+'/data']).shape[0] # SH
    if nucleon_analysis:
        cnfgs_list = np.arange(19,2010,20) # Single Hadrons
    else:
        cnfgs_list = np.arange(0,20,1) # Single Hadrons
    
# Now this is normalizing with respect to the configs that are effectively used. 
worked_weight = []
for ii in range(len(cnfgs_list)):
    worked_weight.append(pre_weight[cnfgs_list[ii]])    

### Final reweighting factors
weight = np.array(vf.RW_NORMALIZATION(worked_weight,ncfgs), dtype=np.float128)

### List of tmaxs used for the fitting procedure. 
if nucleon_analysis:
    listTMaxSingleHads = [25]
else:
    listTMaxSingleHads = [36, 21, 33]#[np.array(f1[name1[0]+'/data']).shape[-1]]*len(name1)

### TMax used for the fits of correlation matrices
listTMaxMultiHads=[]
for ix in range(0,len(name)):
    list_tmax_multihads_ix=[]
    for jx in range(0,np.array(f[name[ix]+'/data']).shape[1]):
        list_tmax_multihads_ix.append(np.array(f[name[ix]+'/data']).shape[-1])
    listTMaxMultiHads.append(list_tmax_multihads_ix)

### Minimum time slices used for the fits of single hadron correlators.
if nucleon_analysis:
    singleTMinsFitPlots = [14]
else:
    singleTMinsFitPlots = [13,9,15] #[10]*len(name1)

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
