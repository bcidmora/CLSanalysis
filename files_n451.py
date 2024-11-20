import h5py
import numpy as np
import os
import set_of_functions as vf

location = os.path.expanduser('~')+'$YOUR_PATH_TO_THE_DATA$'

# ENSEMBLE FILES: N451
f = h5py.File(location+'/N451/cls21_N451_r000_isotriplet_S0_multiple_fwd.hdf5','r')
f1 = h5py.File(location+'/N451/cls21_N451_r000_isotriplet_S0_singles_fwd.hdf5','r')

weight_raw = np.array(np.loadtxt(location + '/N451/N451r000_rw.dat', unpack=True))[3]

### How Many Irreps has each?
name = list(f.keys())
name1 = list(f1.keys())

### Number of gauge configurations
ncfgs = np.array(f[name[0]+'/data']).shape[0]

weight = np.array(vf.RW_NORMALIZATION(weight_raw,ncfgs), dtype=np.float64)

## List of tmaxs used for the fitting procedure. 
listTMaxSingleHads = [np.array(f1[name1[0]+'/data']).shape[-1]]*len(name1)

listTMaxMultiHads=[]
for ix in range(0,len(name)):
    list_tmax_multihads_ix=[]
    for jx in range(0,np.array(f[name[ix]+'/data']).shape[1]):
        list_tmax_multihads_ix.append(np.array(f[name[ix]+'/data']).shape[-1])
    listTMaxMultiHads.append(list_tmax_multihads_ix)


singleTMinsFitPlots = [12]*len(name1)

multiTMinsFitPlots = []
for ix in range(0,len(name)):
    multiTMinsFitPlots_ij=[]
    for jx in range(0,np.array(f[name[ix]+'/data']).shape[1]):
        multiTMinsFitPlots_ij.append(12)
    multiTMinsFitPlots.append(multiTMinsFitPlots_ij)

# ------------ SOME BASIC CONSTANTS TO DO ANALYSIS LATER -----------------------------

aLat = np.float64(0.0749) # in fm
betaLat = np.float64(3.46)
fmToMev = np.float64(197.327) #1 fm = 197.3 Mev
LatSize = np.float64(48.)
