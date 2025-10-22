import h5py
import numpy as np
import os
import set_of_functions as vf

### This is the source directory where the codes, data, etc. 
location = os.path.expanduser('~')+'$YOUR_PATH_TO_THE_DATA$'

hdf5NameMulti = 'cls21_N451_r000_isotriplet_S0_multiple_fwd.hdf5'
hdf5NameSingle = 'cls21_N451_r000_isotriplet_S0_singles_fwd.hdf5'

### ENSEMBLE FILES: N451
f = h5py.File(location+'data/N451/'+hdf5NameMulti,'r')
f1 = h5py.File(location+'data/N451/'+hdf5NameSingle,'r')

### Final reweighting factors
weight_raw = np.array(np.loadtxt(location + 'data/N451/N451r000.ms1.dat_ascii', unpack=True))[1]

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

### Number of gauge configurations available
ncfgs = np.array(f[name[0]+'/data']).shape[0]
# ncfgs = 1011

### Final reweighting factors (normalized to the total)
weight = np.array(vf.RW_NORMALIZATION(weight_raw,ncfgs), dtype=np.float64)

### This is the list of chosen operators 
chosen_operators_list = []

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
singleTMinsFitPlots = [12]*len(name1)

### Minimum time slices for the fits of multihadron correlators
multiTMinsFitPlots = []
for ix in range(0,len(name)):
    multiTMinsFitPlots_ij=[]
    for jx in range(0,np.array(f[name[ix]+'/data']).shape[1]):
        multiTMinsFitPlots_ij.append(12)
    multiTMinsFitPlots.append(multiTMinsFitPlots_ij)

### PARAMETERS OF THE LATTICE
aLat = np.float64(0.0749) # in fm
betaLat = np.float64(3.46)
fmToMev = np.float64(197.327) #1 fm = 197.3 Mev
LatSize = np.float64(48.)
