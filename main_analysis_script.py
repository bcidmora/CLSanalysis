import correlators_script as cs
import effective_masses_script as efs
import eigenvalues_script as evs
import fitting_script as fs
import add_rem_operators as rcs
import sys


### ------- WHAT IT IS DONE --------

# runReduceDataSet = False
runCorrs = False
runEigenvals = False
runSorting = False
runRowsCols = False
runEffMass = False
runFits = False


### ------ MAIN VARIABLES ---------

myEns = str(sys.argv[1]).upper()
myWhichCorrelator = str(sys.argv[2]).lower()
myTypeRs = str(sys.argv[3]).lower()
myRebinOn = str(sys.argv[4]).lower()
myRb = 1
myVersion = '_test' 
myKbt = 500

### This is the amount of irreps to compute or when to start and when to finish the analysis
myNrIrreps = None # 2 # 1
myFirstIrrep = None # 1 # 2
myLastIrrep = None

### Fitting parameters
myTypeFit = '1' #'2' #'g'
myTypeCorrelation =  'Correlated' # 'Uncorrelated'

### GEVP parameters
myOneTMin = False # True
myOneT0 =  True # False
myT0 = 2
mySorting = 'eigenvals' #'eigenvals' # 'vecs_fix' # 'vecs_fix_norm' # 'vecs_var' # 'vecs_var_norm'


### Oher parameters
myKbtSamples = None #np.array(np.loadtxt('bootstrap_samples.txt')) 
myEffMassDistance = 1 #None #2 #3
myOperatorAnalysisMethod = 'from_list' # 'adding' # 'removing' # 'from_list'

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
myChosenOpsList = chosen_operators_list


### -------- PRINTING INFO OF ENSEMBLE ---------

vf.INFO_PRINTING(myWhichCorrelator, myEns)


### ------------ START --------------

##  Single Hadron correlators
if myWhichCorrelator =='s':
    myArchivo, myIrreps = f1, name1 

     ### Correlators analysis
    if runCorrs: 
        locationWorkedCorrelators = cs.SingleCorrelatorAnalysis(myArchivo, myLocation, myVersion, myTypeRs, myIrreps, myWeight, rebin_on = myRebinOn, rb = myRb, kbt = myKbt, number_cfgs = myCnfgs, nr_irreps=myNrIrreps, own_kbt_list = myKbtSamples, first_irrep = myFirstIrrep , last_irrep = myLastIrrep)
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
        
        fs.FitSingleCorrelators(myCorrelator, myFitCorrelator, myTypeRs, listTMaxSingleHads, myIrreps, one_tmin = myOneTMin, type_fit = myTypeFit, type_correlation = myTypeCorrelation)
        
        myFitCorrelator.close()
        
## Multi-Hadron correlators
elif myWhichCorrelator=='m':        
    myArchivo, myIrreps = f, name
    
     ### Correlators analysis
    if runCorrs: 
        locationWorkedCorrelators = cs.MultiCorrelatorAnalysis(myArchivo, myLocation, myVersion, myTypeRs, myIrreps, myWeight, rebin_on = myRebinOn, rb = myRb, kbt = myKbt, number_cfgs = myCnfgs, nr_irreps=myNrIrreps, own_kbt_list=myKbtSamples, first_irrep = myFirstIrrep , last_irrep = myLastIrrep)
    else:
        locationWorkedCorrelators = myLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion
    
    myCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    
    ### GEVP calculation
    if runEigenvals or runRowsCols:
        myT0Min = int(input('T0 min: '))
        myT0Max = int(input('T0 max: ')) 
    
    if runEigenvals:
        evs.EigenvaluesExtraction(myCorrelator, myTypeRs, myIrreps, t0_min = myT0Min, t0_max = myT0Max, sorting=mySorting)
    
    ### Operators Analysis
    if runRowsCols:
        rcs.OperatorsAnalysis(myCorrelator, myTypeRs, myOperatorAnalysisMethod, myIrreps, t0_min = myT0Min, t0_max = myT0Max, ops_analysis_list = myChosenOpsList)
    
    ### Effective Masses analysis
    if runEffMass: 
        efs.MultiCorrelatorEffectiveMass(myCorrelator, myTypeRs, dist_eff_mass = myEffMassDistance)

    ### Fits analysis
    if runFits:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_Matrices/')
        myFitCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        fs.FitMultiCorrelators(myCorrelator, myFitCorrelator, myTypeRs, listTMaxMultiHads, myIrreps,  type_fit = myTypeFit, type_correlation = myTypeCorrelation, one_tmin = myOneTMin, one_t0 = myOneT0, chosen_t0 = myT0, gevp=True, operators_analysis = False, the_operator_analysis_method = myOperatorAnalysisMethod)
        
        myFitCorrelator.close()
        

## Ratio Multi-Hadrons (THIS PART IS STILL UNDER CONSTRUCTION)
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

