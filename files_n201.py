import h5py
import numpy as np
import os
import set_of_functions as vf

### This is the source directory where the codes, data, etc. 
location = os.path.expanduser('~') + '$YOUR_PATH_TO_THE_DATA$'

### Names of the files with single correlators and multi hadrons
hdf5NameMulti = 'cls21_N201_r000_isotriplet_S0_multiple_fwd.hdf5'
hdf5NameSingle = 'cls21_N201_r000_isotriplet_S0_singles_fwd.hdf5'

### ENSEMBLE FILES: N201
f = h5py.File(location+'data/N201/'+hdf5NameMulti,'r')
f1 = h5py.File(location+'data/N201/'+hdf5NameSingle,'r')

### Reweighting factors
# weight_raw = np.array([1.0]*len(f[list(f.keys())[0]+'/data'])) ###OLD
weight_raw = np.array(np.loadtxt(location + 'data/N201/N201r001.ms1.dat_ascii', unpack=True))[1]

### How many irreps has each?
name = list(f.keys())
name1 = list(f1.keys())

### Nr. Gauge Configurations
ncfgs = np.array(f[name[0]+'/data']).shape[0]

### Final reweighting factors
weight = np.array(vf.RW_NORMALIZATION(weight_raw, ncfgs), dtype=np.float64)

### List of tmaxs used for the fitting procedure. 
listTMaxSingleHads = [np.array(f1[name1[0]+'/data']).shape[-1]]*len(name1)

listTMaxMultiHads=[]
for ix in range(0,len(name)):
    list_tmax_multihads_ix=[]
    for jx in range(0,np.array(f[name[ix]+'/data']).shape[1]):
        list_tmax_multihads_ix.append(np.array(f[name[ix]+'/data']).shape[-1])
    listTMaxMultiHads.append(list_tmax_multihads_ix)
    
singleTMinsFitPlots = [6]*len(name1)

multiTMinsFitPlots = []
for ix in range(0,len(name)):
    multiTMinsFitPlots_ij=[]
    for jx in range(0,np.array(f[name[ix]+'/data']).shape[1]):
        multiTMinsFitPlots_ij.append(6)
    multiTMinsFitPlots.append(multiTMinsFitPlots_ij)

###PARAMETERS OF THE LATTICE from arXiv:2112.06696 [hep-lat]
aLat = np.float64(0.0749) # in fm
betaLat = np.float64(3.55)
fmToMev = np.float64(197.327) #1 fm = 197.3 Mev
LatSize = np.float64(48.)
