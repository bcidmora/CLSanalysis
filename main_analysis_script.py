import correlators_script as cs
import effective_masses_script as efs
import eigenvalues_script as evs
import fitting_script as fs
import sys



### ------- WHAT IT IS DONE --------

runCorrs = True
runEigenvals = True
runEffMass = True
runFits = True



### ------ MAIN VARIABLES ---------

myEns = str(sys.argv[1]).upper()
myWhichCorrelator = str(sys.argv[2]).lower()
myTypeRs = str(sys.argv[3]).lower()
myRebinOn = str(sys.argv[4]).lower()
myRb = 2
myVersion = 'test'
myKbt = 500
myTypeFit = '1' #'1'
myTypeCorrelation = 'Correlated' # 'Uncorrelated'
myOneTMin = True # False
myOneT0 =  True # False
myT0 = 4
myNrIrreps = None # 1
myKbtSamples = None #np.array(np.loadtxt('bootstrap_samples.txt')) 


if myRebinOn=='rb': 
    reBin = '_bin'+str(myRb)
else: 
    reBin = ''



### ----  GETTING THE INFO FROM THE FOLLOWING FILES -------

if myEns == 'N451': from files_n451 import *
elif myEns == 'N201': from files_n201 import * 
elif myEns == 'D200': from files_d200 import *

myLocation = vf.DIRECTORY_EXISTS(location + '$YOUR_OUTPUT_PATH$/%s/'%myEns)
myWeight = weight
myCnfgs = ncfgs



### -------- PRINTING INFO OF ENSEMBLE ---------

vf.INFO_PRINTING(myWhichCorrelator, myEns)


### ------------ START --------------
##  Single Hadron correlators
if myWhichCorrelator =='s':
    myArchivo, myIrreps = f1, name1 

    if runCorrs: 
        locationWorkedCorrelators = cs.SingleCorrelatorAnalysis(myArchivo, myLocation, myVersion, myTypeRs, myIrreps, myWeight, rebin_on = myRebinOn, rb = myRb, kbt = myKbt, number_cfgs = myCnfgs, nr_irreps=myNrIrreps, own_kbt_list = myKbtSamples)
    else:
        locationWorkedCorrelators = myLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion
    
    myCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    
    if runEffMass: 
        efs.SingleCorrelatorEffectiveMass(myCorrelator, myTypeRs) 
        
    if runFits: 
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_SingleHadrons/')
        myFitCorrelator =  h5py.File(myFitsLocation + 'Single_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        fs.FitSingleCorrelators(myCorrelator, myFitCorrelator, myTypeRs, listTMaxSingleHads, one_tmin = myOneTMin, type_fit = myTypeFit, type_correlation = myTypeCorrelation)
        
        
        myFitCorrelator.close()
        
## Multi-Hadron correlators
elif myWhichCorrelator=='m':        
    myArchivo, myIrreps = f, name
    
    if runCorrs: 
        locationWorkedCorrelators = cs.MultiCorrelatorAnalysis(myArchivo, myLocation, myVersion, myTypeRs, myIrreps, myWeight, rebin_on = myRebinOn, rb = myRb, kbt = myKbt, number_cfgs = myCnfgs, nr_irreps=myNrIrreps, own_kbt_list=myKbtSamples)
    else:
        locationWorkedCorrelators = myLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion
    
    myCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    
    if runEigenvals:
        myT0Min = int(input('T0 min: '))
        myT0Max = int(input('T0 max: ')) 
        evs.EigenvaluesExtraction(myCorrelator, myTypeRs, t0_min = myT0Min, t0_max = myT0Max)
    
    if runEffMass: 
        efs.MultiCorrelatorEffectiveMass(myCorrelator, myTypeRs)
        
    if runFits:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_Matrices/')
        myFitCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        fs.FitMultiCorrelators(myCorrelator, myFitCorrelator, myTypeRs, listTMaxMultiHads, type_fit = myTypeFit, type_correlation = myTypeCorrelation, one_tmin = myOneTMin, one_t0 = myOneT0, chosen_t0 = myT0)
        
        myFitCorrelator.close()
        

## Ratio Multi-Hadrons 
elif myWhichCorrelator=='mr':        
    locationWorkedCorrelators = myLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin + '_v%s.h5'%myVersion
    
    myRatioCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    
    if runEffMass: 
        efs.MultiCorrelatorEffectiveMass(myRatioCorrelator, myTypeRs)
        
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

