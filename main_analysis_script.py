import correlators_script as cs
import effective_masses_script as efs
import eigenvalues_script as evs
import fitting_script as fs
import rows_cols_script as rcs
import sys


### ------- WHAT IT IS DONE --------

runCorrs = True
runEigenvals = False
runRowsCols = False
runEffMass = False
runFits = False


### ------ MAIN VARIABLES ---------

myEns = str(sys.argv[1]).upper()
myWhichCorrelator = str(sys.argv[2]).lower()
myTypeRs = str(sys.argv[3]).lower()
myRebinOn = str(sys.argv[4]).lower()
myRb = 1
myVersion = 'test'
myKbt = 500
myNrIrreps = 1 #None # 2

### Fitting parameters
myTypeFit = '1' #'2'
myTypeCorrelation =  'Correlated' # 'Uncorrelated'

### GEVP parameters
myOneTMin = False # True
myOneT0 =  True # False
myT0 = 3
mySorting = 'eigenvals' # 'vecs_fix_norm' # 'vecs_var' # 'vecs_var_norm'# 'vecs_fix_norm'

### Oher parameters
myKbtSamples = None #np.array(np.loadtxt('bootstrap_samples.txt')) 
myEffMassDistance = 1 #None #2 #3
myOpAnalysis = False

if myRebinOn=='rb': 
    reBin = '_bin'+str(myRb)
else: 
    reBin = ''

### ----  GETTING THE INFO FROM THE FOLLOWING FILES -------

if myEns == 'N451': from files_n451 import *
elif myEns == 'N201': from files_n201 import * 
elif myEns == 'D200': from files_d200 import *
elif myEns == 'X451': from files_x451 import *

myLocation = vf.DIRECTORY_EXISTS(location + '$YOUR_OUTPUT_PATH$/%s/'%myEns)
myWeight = weight
myCnfgs = ncfgs # None # 20 # 100


### -------- PRINTING INFO OF ENSEMBLE ---------

vf.INFO_PRINTING(myWhichCorrelator, myEns)


### ------------ START --------------
##  Single Hadron correlators
if myWhichCorrelator =='s':
    myArchivo, myIrreps = f1, name1 
    
    ### Correlators analysis
    if runCorrs: 
        locationWorkedCorrelators = cs.SingleCorrelatorAnalysis(myArchivo, myLocation, myVersion, myTypeRs, myIrreps, myWeight, rebin_on = myRebinOn, rb = myRb, kbt = myKbt, number_cfgs = myCnfgs, nr_irreps=myNrIrreps, own_kbt_list = myKbtSamples)
    else:
        locationWorkedCorrelators = myLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion
    
    myCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    
    ### Effective Masses analysis
    if runEffMass: 
        efs.SingleCorrelatorEffectiveMass(myCorrelator, myTypeRs, dist_eff_mass=myEffMassDistance) 
        
    ### Fits analysis
    if runFits: 
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_SingleHadrons/')
        myFitCorrelator =  h5py.File(myFitsLocation + 'Single_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        fs.FitSingleCorrelators(myCorrelator, myFitCorrelator, myTypeRs, listTMaxSingleHads, one_tmin = myOneTMin, type_fit = myTypeFit, type_correlation = myTypeCorrelation)
        
        myFitCorrelator.close()
        
## Multi-Hadron correlators
elif myWhichCorrelator=='m':        
    myArchivo, myIrreps = f, name
    
    ### Correlators analysis
    if runCorrs: 
        locationWorkedCorrelators = cs.MultiCorrelatorAnalysis(myArchivo, myLocation, myVersion, myTypeRs, myIrreps, myWeight, rebin_on = myRebinOn, rb = myRb, kbt = myKbt, number_cfgs = myCnfgs, nr_irreps=myNrIrreps, own_kbt_list=myKbtSamples)
    else:
        locationWorkedCorrelators = myLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion
    
    myCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    
    ### GEVP calculation
    if runEigenvals or runRowsCols:
        myT0Min = int(input('T0 min: '))
        myT0Max = int(input('T0 max: ')) 
    
    if runEigenvals:
        evs.EigenvaluesExtraction(myCorrelator, myTypeRs, t0_min = myT0Min, t0_max = myT0Max, sorting=mySorting)
    
    ### Operators Analysis
    if runRowsCols:
        # rcs.REMOVING_COLS_ROWS(myCorrelator, myTypeRs, t0_min = myT0Min, t0_max = myT0Max)#, nr_irreps=myNrIrreps)
        rcs.ADDING_COLS_ROWS(myCorrelator, myTypeRs, t0_min = myT0Min, t0_max = myT0Max)
    
    ### Effective Masses analysis
    if runEffMass: 
        # if runRowsCols: myOpAnalysis = True
        efs.MultiCorrelatorEffectiveMass(myCorrelator, myTypeRs, dist_eff_mass = myEffMassDistance)#, op_analysis=myOpAnalysis)

    ### Fits analysis
    if runFits:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_Matrices/')
        myFitCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        fs.FitMultiCorrelators(myCorrelator, myFitCorrelator, myTypeRs, listTMaxMultiHads, type_fit = myTypeFit, type_correlation = myTypeCorrelation, one_tmin = myOneTMin, one_t0 = myOneT0, chosen_t0 = myT0)
        
        myFitCorrelator.close()
        

## Ratio Multi-Hadrons 
elif myWhichCorrelator=='mr':        
    locationWorkedCorrelators = myLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin + '_v%s.h5'%myVersion
    
    myRatioCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    
    ### Effective Masses analysis
    if runEffMass: 
        efs.MultiCorrelatorEffectiveMass(myRatioCorrelator, myTypeRs)
    
    ### Fits analysis            
    if runFits:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_Ratios/')
        myFitRatioCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        fs.FitMultiCorrelators(myRatioCorrelator, myFitRatioCorrelator, myTypeRs, listTMaxMultiHads, type_fit = myTypeFit, type_correlation = myTypeCorrelation, one_tmin = True, one_t0 = True, chosen_t0 = myT0, nr_irreps=myNrIrreps)
        
        myFitCorrelator.close()

else: 
    print('Not proper choice.')
    sys.exit()
    

### --------- PRINTS WHERE IT IS SAVED --------        

print('-'*(len(locationWorkedCorrelators)+1))
print('Correlator analysis saved as: \n' + locationWorkedCorrelators )
print('_'*(len(locationWorkedCorrelators)+1))


### -------- CLOSES ALL OTHER FILES --------
myArchivo.close()
myCorrelator.close()

