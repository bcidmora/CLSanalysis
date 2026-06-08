### MAIN FUNCTIONS CARRING OUT THE JOBS
import correlators_script as cs
import effective_masses_script as efs
import eigenvalues_script as evs
import fitting_script as fts
import add_rem_operators as rcs

### MAIN ANALYSIS FUNCTIONS
import bin_size_script as bs
import dispersion_relation_script as drs
import correlated_differences_script as cds

### LAYOUT AND EXTRA FUNCTIONS
import set_of_layout_functions as vfl
import set_of_analysis_functions as vfa

### EXTRA LIBRARIES
import ast
import numpy as np
import os
import h5py
import sys

import ensembles as ed

### ------ EXTRA VARIABLES ----- ### TO CHECK 
# runSorting = False
# runReduceDataSet = False
# runPlotBinOnly = False

### ------- WHAT IT IS DONE --------
myArgs = vfl.parse_args()
myRuns = vfl.WhichRuns(myArgs, ed.ensembles)

### Other parameters
myKbtSamples = None #np.array(np.loadtxt('bootstrap_samples.txt')) #None

### ------ CHECKING EVERYTHING ---------
vfl.VALIDATE_RUNS(myRuns)

### This is just naming schemes
reBin = f"_bin{myRuns.rb}" if myRuns.rebin else ""

### Information from the ensembles dictionary
myLocation = vfl.DIRECTORY_EXISTS(f'{ed.outputLocation}{myRuns.ensemble}/')

### ----  GETTING THE INFO FROM THE FOLLOWING FILES -------
myCnfgs = ed.ensembles[myRuns.ensemble]['ncfgs']
myWeight = vfa.REWEIGHTS(ed.ensembles[myRuns.ensemble]['weight_raw'], myCnfgs)

if not ed.ensembles[myRuns.ensemble]['allConfigs']:
    myTempCnfgs = ed.ensembles[myRuns.ensemble]['nfgsList']
    myCnfgs = len(myTempCnfgs)
    myWeight = np.asarray(vfa.RW_NORMALIZATION([myWeight[ii] for ii in myTempCnfgs], myCnfgs), dtype=np.float128)

myLatSize = ed.ensembles[myRuns.ensemble]['LatSize']

### -------- PRINTING INFO OF ENSEMBLE ---------
vfl.INFO_PRINTING(myRuns.correlator, myRuns.ensemble)

### ------------ START --------------
##  Single Hadron correlators
if myRuns.correlator =='s':
    myVersion =  f'{myRuns.ensemble}_singles_test' 
    myArchivo = h5py.File(ed.ensembles[myRuns.ensemble]['fs'], 'r')
    myIrreps = list(myArchivo.keys())        
    
    ### Correlators analysis
    if myRuns.corrs:         
        locationWorkedCorrelators = cs.SingleCorrelatorAnalysis(myArchivo, myLocation, myVersion, myRuns.rs_type, myIrreps, myWeight, rebin_on = myRuns.rebin, rb = myRuns.rb, kbt = myRuns.kbt, nr_irreps = myRuns.the_irreps.nr_irreps, own_kbt_list = myKbtSamples, first_irrep = myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep, number_cfgs = myRuns.the_configs.nr_configs)
    else:
        locationWorkedCorrelators = f'{myLocation}Single_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5'
    
     try:
        myCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    except FileNotFoundError:
        sys.exit(f'Cannot find the correlator file for further analysis. Check path: \n {locationWorkedCorrelators}')
        
    ### Effective Masses analysis
    if myRuns.effmass: 
        efs.SingleCorrelatorEffectiveMass(myCorrelator, myRuns.rs_type, dist_eff_mass = myRuns.dist_eff_mass) 
        
    ### Fits analysis
    if myRuns.fits: 
        myTMaxList = ed.ensembles[myRuns.ensemble]['singleTMaxFits']
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myLocation}Fits_SingleHadrons/')
        myFitCorrelator =  h5py.File(f'{myFitsLocation}Single_correlators_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'a')
        
        fts.FitSingleCorrelators(myCorrelator, myFitCorrelator, myRuns.rs_type, myTMaxList, myIrreps, myRuns.fit.type_fit, myRuns.fit.type_correlation, one_tmin = myRuns.fit.one_tmin, first_irrep = myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep)
        
        myFitCorrelator.close()
    
    ### Dispersion relation analysis
    if myRuns.disp:
        myDispMode = myRuns.disp_run.mode
        myTMinSPlot = myArchivoPre[f'singleTMinResults-{myRuns.fit_type}exp']
        myDispLocation = vfl.DIRECTORY_EXISTS(f'{ed.location}/Plots/{myRuns.ensemble}/Dispersion_Relation/')
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myLocation}Fits_SingleHadrons/')
        myDispFile  = h5py.File(f'{myFitsLocation}Single_correlators_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'a')
        
        drs.DISPERSION_RELATION(myDispFile, myRuns.rs_type, myTMinSPlot, myLatSize, myDispLocation, myDispMode, myRuns.ensemble)
        
        myDispFile.close()
        
        ### THIS NEEDS TO BE CHECKED
#     ### Binning analysis
#     # if any([myRuns.corrs, runBinSizeAnalysis]):
#         # vfl.SINGLE_INFO_PRINTING(myIrreps, runBinSizeAnalysis, myRuns.corrs)
#     if runBinSizeAnalysis:
#         myBinSizeIrrep = ast.literal_eval(input('Enter irrep for bin-size analysis: []'))
#         myBinSizeIrrepName = myIrreps[myBinSizeIrrep[0]]
#         myBinLocation = vfl.DIRECTORY_EXISTS(f'{myLocation}Bin_Size_Analysis_{myBinSizeIrrepName}_{myRuns.rs_type}/')
#         
#         bs.BinSizeAnalysis(myArchivo, myBinLocation, myRuns.rs_type, myBinSizeIrrep, myIrreps, myTMinFit, myMaxBinSize, myBinSizeFitRange, myChosenBinSize, myVersion, myWeight, one_tmin=myRuns.fit.one_tmin, type_fit=myRuns.fit.type_fit, type_correlation=myRuns.fit.type_correlation, kbt=myRuns.kbt, number_cfgs=myCnfgs, own_kbt_list=myKbtSamples, isospin_label=ed.ensembles[myRuns.ensemble]['iso_label'], ensemble=myRuns.ensemble, plots_only=runPlotBinOnly)
#     
    myArchivo.close()
    myCorrelator.close()
        
## Multi-Hadron correlators
elif myRuns.correlator=='m':  
    myIsospin = myRuns.isospin
    myChosenIsospin = ed.ensembles[myRuns.ensemble][myIsospin]['iso_tag']
    myArchivo = h5py.File(ed.ensembles[myRuns.ensemble][myIsospin]['fm'], 'r')
    myIrreps = list(myArchivo.keys())
    myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_test' 
    
    ### Correlators analysis
    if myRuns.corrs: 
        locationWorkedCorrelators = cs.MultiCorrelatorAnalysis(myArchivo, myLocation, myVersion, myRuns.rs_type, myIrreps, myWeight, rebin_on = myRuns.rebin, rb = myRuns.rb, kbt = myRuns.kbt, nr_irreps = myRuns.the_irreps.nr_irreps, own_kbt_list = myKbtSamples, first_irrep = myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep, number_cfgs = myRuns.the_configs.nr_configs)
    else:
        locationWorkedCorrelators = f'{myLocation}Matrix_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5'
    
    try:
        myCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    except FileNotFoundError:
        print('Cannot find the correlator file for the further analysis steps.')
        sys.exit()
    
    ### GEVP calculation
    if myRuns.eigenvals:        
        mySorting = myRuns.gevp.sorting
        myTD = myRuns.gevp.td
        myRsSorting = myRuns.gevp.rs_sorting
        
        myT0Min = myRuns.gevp.t0min
        myT0Max = myRuns.gevp.t0max

        evs.EigenvaluesExtraction(myCorrelator, myRuns.rs_type, myIrreps, myT0Min, myT0Max, sorting = mySorting, the_td = myTD, rs_sorting = myRsSorting)
    
    ### Reduced operator set
    if myRuns.ops:            
        mySorting = myRuns.gevp.sorting
        myTD = myRuns.gevp.td
        myRsSorting = myRuns.gevp.rs_sorting
        
        myT0Min = myRuns.gevp.t0min
        myT0Max = myRuns.gevp.t0max
        
        myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_reduced' 
        locationWorkedCorrelators = f'{myLocation}Matrix_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5'
        myChosenOpsList = ed.ensembles[myRuns.ensemble][myIsospin]['operatorsChoice']
        rcs.OperatorsAnalysis(myCorrelator, myRuns.rs_type, myIrreps, myT0Min, myT0Max, locationWorkedCorrelators, myChosenOpsList, sorting = mySorting, the_td = myTD, rs_sorting = myRsSorting)
    
    ### Effective Masses analysis
    if myRuns.effmass:
        if myRuns.ops_flag or myRuns.ops:
            myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_reduced' 
            locationWorkedCorrelators = f'{myLocation}Matrix_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5'
            myCorrelator = h5py.File(locationWorkedCorrelators, 'r+')     
            efs.MultiCorrelatorEffectiveMass(myCorrelator, myRuns.rs_type, dist_eff_mass = myRuns.dist_eff_mass)
        else:
            efs.MultiCorrelatorEffectiveMass(myCorrelator, myRuns.rs_type, dist_eff_mass = myRuns.dist_eff_mass)

    ### Fits analysis
    if myRuns.fits:
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myLocation}Fits_Matrices/')
        
        if myRuns.ops_flag or myRuns.ops:
            myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_reduced' 
            locationWorkedCorrelators = f'{myLocation}Matrix_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5'
            myCorrelator = h5py.File(locationWorkedCorrelators, 'r+')     
            myTMaxList = ed.ensembles[myRuns.ensemble][myIsospin]['multiTMaxFitsOpChoices']
        else:
            myTMaxList = ed.ensembles[myRuns.ensemble][myIsospin]['multiTMaxFits']
            
        myFitCorrelator = h5py.File(f'{myFitsLocation}Matrix_correlators_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'a')
        fts.FitMultiCorrelators(myCorrelator, myFitCorrelator, myRuns.rs_type, myTMaxList, myIrreps, myRuns.fit.type_fit, myRuns.fit.type_correlation, one_tmin = myRuns.fit.one_tmin, one_t0 = myRuns.fit.one_t0, chosen_t0 = myRuns.fit.t0, first_irrep = myRuns.the_irreps.first_irrep , last_irrep = myRuns.the_irreps.last_irrep)
        
        myFitCorrelator.close()
    
    ### Still under construction
    # myNonInteractingLevels = ed.ensembles[myRuns.ensemble][myIsospin]['nonInteractingLevels']
    # if runCorrDiff:
        # myNonInteractingFileDict = ed.ensemble_dictionary[myRuns.ensemble]['singleMesons']
        # myFitsLocation = vfa.DIRECTORY_EXISTS(myAnalysisLocation + 'Fits_Matrices/')
        # myFitCorrelator = h5py.File(
            # myFitsLocation + f'Matrix_correlators_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'r')
        # myCorrDiffLocation = vfa.DIRECTORY_EXISTS(myAnalysisLocation + 'Correlated_Differences/')
        # cds.CORRELATED_DIFFERENCES(myNonInteractingLevels, myTMinPlot, myTMinMPlot, myNonInteractingFileDict, myFitCorrelator, myCorrDiffLocation, myTMinFit, myLatSize,myRuns.the_irreps.nr_irreps, myNrOps, myRuns.fit.t0, myTD, myRuns.rs_type, myCorrDiffPlot, myThresholdMasses, myRuns.ensemble, myIsoLabel, myTwoBody, myThreeBody, fmToMev, myRuns.rb)
    
    myArchivo.close()
    myCorrelator.close()


### Ratio Multi-Hadrons
elif myRuns.correlator=='mr':     
    myIsospin = myRuns.isospin
    myNonInteractingLevels = ed.ensembles[myRuns.ensemble][myIsospin]['nonInteractingLevels']
    myChosenIsospin = ed.ensembles[myRuns.ensemble][myIsospin]['iso_tag']
    myTypeCorrelation = myRuns.fit.type_correlation
    
    ### Multi-hadron corrs
    if myRuns.ops_flag:
        myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_reduced' 
    else:
        myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_test' 
        
    ### Ratio of Correlators
    if myRuns.corrs:
        
        ### Multi-hadron corrs
        myArchivo = h5py.File(ed.ensembles[myRuns.ensemble][myIsospin]['fm'], 'r')
        myMultiIrreps = list(myArchivo.keys())
        myArchivo.close()
        myMHCorrelator = h5py.File(f'{myLocation}Matrix_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5', 'r+')
        
        ### Single-hadron corrs
        mySHVersion =  f'{myRuns.ensemble}_singles_test' 
        myArchivo = h5py.File(ed.ensembles[myRuns.ensemble]['fs'], 'r')
        mySingleIrreps = list(myArchivo.keys())
        myArchivo.close()
        mySHCorrelator = h5py.File(f'{myLocation}Single_correlators_{myRuns.rs_type}{reBin}_{mySHVersion}.h5', 'r+')               
        
        ### Taking the ratios
        locationWorkedCorrelators = cs.MultiCorrelatorAnalysisRatios(myMHCorrelator, mySHCorrelator, myLocation, myNonInteractingLevels, myVersion, myRuns.rs_type, myMultiIrreps, mySingleIrreps, rebin_on = myRuns.rebin, rb = myRuns.rb, nr_irreps = myRuns.the_irreps.nr_irreps, first_irrep = myRuns.the_irreps.first_irrep , last_irrep = myRuns.the_irreps.last_irrep)
        
        myMHCorrelator.close()
        mySHCorrelator.close()
    else:
        locationWorkedCorrelators = f'{myLocation}Ratio_matrix_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5'
        
    try:
        myRatioCorrelator = h5py.File(locationWorkedCorrelators, 'r+')
    except FileNotFoundError:
        sys.exit(f"File not found. Check directory:\n{locationWorkedCorrelators}")
    
    ### Effective Masses analysis
    if myRuns.effmass: 
        efs.RatioMultiCorrelatorEffectiveMass(myRatioCorrelator, myRuns.rs_type)
    
    ### Fits analysis
    if myRuns.fits:          
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myLocation}Fits_Ratios/')
        myTMaxList = ed.ensembles[myRuns.ensemble][myIsospin]['multiTMaxFitsRatios']
        
        ### Multi-hadron corrs
        myArchivo = h5py.File(ed.ensembles[myRuns.ensemble][myIsospin]['fm'], 'r')
        myMultiIrreps = list(myArchivo.keys())
        myArchivo.close()
        
        myFitRatioCorrelator = h5py.File(f'{myFitsLocation}Ratio_matrix_correlators_{myChosenIsospin}_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'a')
        
        fts.FitRatioMultiCorrelators(myRatioCorrelator, myFitRatioCorrelator, myRuns.rs_type, myTMaxList, myMultiIrreps, myRuns.fit.type_fit, myRuns.fit.type_correlation, one_tmin = myRuns.fit.one_tmin, one_t0 = myRuns.fit.one_t0, chosen_t0 = myRuns.fit.t0, first_irrep = myRuns.the_irreps.first_irrep , last_irrep = myRuns.the_irreps.last_irrep)
        
        myFitRatioCorrelator.close()
    
    myRatioCorrelator.close()    
else: 
    print('Not proper choice.')
    sys.exit()


### --------- PRINTS WHERE IT IS SAVED --------        

print('-'*(len(locationWorkedCorrelators)+1))
print('Correlator analysis saved as: \n' + locationWorkedCorrelators )
print('_'*(len(locationWorkedCorrelators)+1))
