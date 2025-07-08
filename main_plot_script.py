import plot_correlators_script as pcorr
import plot_effective_masses_script as peff
import plot_fits_script as pfit
import set_of_functions as vf
import sys
import os
import h5py
from PyPDF2 import PdfMerger

### ------- WHAT IT IS DONE --------

plotCorrs = False  ### Plot correlators
plotEffMass = True  ### Plot effective masses
plotFits = False ### Plot Fits
joinPlots = False ### Puts all the plots found in one pdf file


### ------ MAIN VARIABLES ---------

myEns = str(sys.argv[1]).upper()
myWhichCorrelator = str(sys.argv[2]).lower()
myTypeRs = str(sys.argv[3]).lower()
myRebinOn = str(sys.argv[4]).lower()

myRb = 1
myVersion = '_test'
myNrExponentials = '1'
myTypeCorrelation = 'Correlated' # 'Uncorrelated'
myOneTmin = True
myT0 = 4

myNrIrreps = None # None # 2 # 1
myFirstIrrep = None # 1 # 2
myLastIrrep = None

### This is for the multihadorn operatos:
myDiagonalCorrs = True ### Plots the diagonal of the correlators
myGevpFlag = False ### Plots the eigenvalues
myOperatorsFlag = False ### Plots the eigenvalues from the operators analysis
myOperatorsMethod = 'from_list' # 'adding' # 'removing'

### This is for all the fit plots
myZoomFit = True ### It does Zoom to the fitting range
myChiPlots = True ### Plotting the Chi^{2/dof}
myTotalChiPlots = True ### Plotting Total Chi^{2}
myDeltaChiPlots = True ### Plotting delta Chi^{2} form time slice to time slice

myDataLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH(SAME_THAN_CORRS_SCRIPT_OUTPUT)$/%s/'%myEns)

### Make sure you have all these variables defined in files_ens.py
if myEns == 'N451': from files_n451 import singleTMinsFitPlots, multiTMinsFitPlots, name, name1
elif myEns == 'N201': from files_n201 import singleTMinsFitPlots, multiTMinsFitPlots, name, name1 
elif myEns == 'D200': from files_d200 import singleTMinsFitPlots, multiTMinsFitPlots, name, name1
elif myEns == 'X451': from files_x451 import singleTMinsFitPlots, multiTMinsFitPlots, name, name1

### This is the info about the file, if it was rebinned
if myRebinOn=='rb': 
    reBin = '_bin'+str(myRb)
else: 
    reBin = ''

### This is the information about the resampling
if myTypeRs=='jk':
    myResamplingScheme='Jackknife'
elif myTypeRs=='bt':
    myResamplingScheme='Bootstrap' 


### -------- PRINTING INFO OF ENSEMBLE ---------

vf.INFO_PRINTING(myWhichCorrelator, myEns)

### ------------ START ----------------

if myWhichCorrelator =='s':
    ### Original list of irreps
    myIrreps = name1
    
    ### Correlators data
    mySingleCorrelatorData = h5py.File(myDataLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion,'r')
    
    ### Directory where the plots will be saved
    myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOU_PLOTS_DIRECTORY$/Plots/%s/SingleHadrons/'%myEns +  '%s/'%myResamplingScheme)
    
    ### Plots all the correlators. Look at the booleans here
    if plotCorrs:        
        pcorr.PlotSingleHadronCorrelators(mySingleCorrelatorData, myTypeRs, myVersion, myPlotLocation, reBin, nr_irreps=myNrIrreps, first_irrep=myFirstIrrep, last_irrep = myLastIrrep)
    
    ### Plots effective masses of single hadrons
    if plotEffMass: 
        peff.PlotSingleHadronsEffectiveMasses(mySingleCorrelatorData, myResamplingScheme, myVersion, myPlotLocation, reBin, nr_irreps=myNrIrreps, first_irrep=myFirstIrrep, last_irrep = myLastIrrep)
    
    ### Plots fits of single hadrons
    if plotFits: 
        myFitsLocation = vf.DIRECTORY_EXISTS(myDataLocation + 'Fits_SingleHadrons/')
        myFitCorrelator =  h5py.File(myFitsLocation + 'Single_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        pfit.PlotSingleHadronsFits(myFitCorrelator, myTypeCorrelation, myNrExponentials,  singleTMinsFitPlots, myVersion, myPlotLocation, reBin, myIrreps, first_irrep=myFirstIrrep, last_irrep = myLastIrrep, zoom_fit=myZoomFit, chi_plots=myChiPlots, total_chi=myTotalChiPlots, delta_chi=myDeltaChiPlots)
        
        myFitCorrelator.close()
        
    ### Puts all the plots in one PDF file. It checks if the file exists first
    if joinPlots:
        ### Loop over all the irreps in this ensemble
        irreps = list(mySingleCorrelatorData.keys())
        for aa in irreps:
            ##3 Loop over all the operators in this ensemble
            ops = list(mySingleCorrelatorData[aa+'/Operators'])
            x=[]
                
            ### Corrs
            if os.path.isfile(myPlotLocation + 'Correlator_' + aa[:4] +'_%s'%aa[-1] + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'Correlator_' + aa[:4] +'_%s'%aa[-1] + reBin + '_v%s.pdf'%myVersion)
            
            ### Corrs log-plots
            if os.path.isfile(myPlotLocation + 'Correlator_' + aa[:4] +'_%s'%aa[-1] + '_log' +  reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'Correlator_' + aa[:4] +'_%s'%aa[-1] + '_log' +  reBin + '_v%s.pdf'%myVersion)
                
            ### Histogram Corrs
            if os.path.isfile(myPlotLocation + 'Histogram_correlators_' + aa[:4] +'_%s'%aa[-1] + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'Histogram_correlators_' + aa[:4] +'_%s'%aa[-1] + reBin + '_v%s.pdf'%myVersion)
            
            ### Effective Masses Corrs
            if os.path.isfile(myPlotLocation + 'EffectiveMass_' + aa[:4] +'_%s'%aa[-1] + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'EffectiveMass_' + aa[:4] +'_%s'%aa[-1] + reBin + '_v%s.pdf'%myVersion)
            
            ### Fits Corrs
            if os.path.isfile(myPlotLocation + 'Tmin_Fits_' + aa[:4] +'_%s'%aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'Tmin_Fits_' + aa[:4] +'_%s'%aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)

            ### Zoom Fits Corrs
            if os.path.isfile(myPlotLocation + 'Tmin_Fits_Zoom_' + aa[:4] +'_%s'%aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'Tmin_Fits_Zoom_' + aa[:4] +'_%s'%aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
            
            ### Chi^{2} Fits Corrs
            if os.path.isfile(myPlotLocation + 'Tmin_Chisqr_' + aa[:4] +'_%s'%aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'Tmin_Chisqr_' + aa[:4] +'_%s'%aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
        
            ### Total Chi^{2} Fits Corrs
            if os.path.isfile(myPlotLocation + 'Tmin_TotalChisqr_' + aa[:4] +'_%s'%aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'Tmin_TotalChisqr_' + aa[:4] +'_%s'%aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
                
            ### Delta Chi^{2} Fits Corrs
            if os.path.isfile(myPlotLocation + 'Tmin_DeltaChisqr_' + aa[:4] +'_%s'%aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'Tmin_DeltaChisqr_' + aa[:4] +'_%s'%aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
                
            merger = PdfMerger()
            for pdf in x:
                merger.append(open(pdf, 'rb'))
            
            with open(myPlotLocation + myEns + "_%s"%aa +  reBin + "_%s"%myTypeRs + "_v%s.pdf"%myVersion, "wb") as fout:
                merger.write(fout)
                
        print('Now all the plots are in one file for each irrep')
    
    mySingleCorrelatorData.close()
    
elif myWhichCorrelator=='m':     
    ### The original list of irreps
    myIrreps = name
    
    ### The correlators data
    myMatrixCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
    
    ### Directory where the plots will be saved
    myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOU_PLOTS_DIRECTORY$/Plots/%s/Matrices/'%myEns +  '%s/'%myResamplingScheme)
    
    ### Plots the correlators. Look at the booleans
    if plotCorrs: 
        pcorr.PlotMultiHadronCorrelators(myMatrixCorrelatorData, myTypeRs, myVersion, myT0, myPlotLocation, reBin, nr_irreps=myNrIrreps, first_irrep=myFirstIrrep, last_irrep = myLastIrrep, diag_corrs= myDiagonalCorrs, gevp=myGevpFlag, ops_analysis=myOperatorsFlag)
    
    ### Plots the effective masses of the eigenvalues from the GEVP and/or the operators analysis
    if plotEffMass: 
        peff.PlotMultiHadronsEffectiveMasses(myMatrixCorrelatorData, myResamplingScheme, myVersion, myT0, myPlotLocation, reBin,nr_irreps=myNrIrreps, first_irrep=myFirstIrrep, last_irrep = myLastIrrep, diag_corrs= myDiagonalCorrs, gevp=myGevpFlag, ops_analysis=myOperatorsFlag)
        
    if plotFits:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myDataLocation + 'Fits_Matrices/')
        myFitCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        pfit.PlotMultiHadronsFits(myFitCorrelator, myTypeCorrelation, myNrExponentials, myTypeRs, multiTMinsFitPlots, myT0, myVersion, myPlotLocation, reBin, myIrreps, gevp=myGevpFlag, zoom_fit=myZoomFit, chi_plots=myChiPlots, total_chi=myTotalChiPlots, delta_chi=myDeltaChiPlots, ops_analysis=myOperatorsFlag, ops_analysis_method=myOperatorsMethod)
        
        
        myFitCorrelator.close()

    ### This part puts all the plots together in one pdf file. 
    if joinPlots:
         ### Loop over all the irreps in this ensemble
        irreps = list(myMatrixCorrelatorData.keys())
        for aa in irreps:
            ##3 Loop over all the operators in this ensemble
            ops = list(myMatrixCorrelatorData[aa+'/Operators'])
            x=[]
            for bb in range(len(ops)):
                
                ### Diagonal Corrs
                if os.path.isfile(myPlotLocation + 'DiagonalCorrelator_' + aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'DiagonalCorrelator_' + aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion)
                
                ### Diagonal corrs log-plot
                if os.path.isfile(myPlotLocation + 'DiagonalCorrelator_' + aa + '_%s_log'%str(bb) + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'DiagonalCorrelator_' + aa + '_%s_log'%str(bb) + reBin + '_v%s.pdf'%myVersion)
                
                ### Histogram diagonal corrs
                if os.path.isfile(myPlotLocation + 'Histogram_DiagCorrelator_' + aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'Histogram_DiagCorrelator_' + aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion)
               
                ### Effective Mass diagonal corrs
                if os.path.isfile(myPlotLocation + 'EffectiveMass_DiagonalCorrelators_'+ aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'EffectiveMass_DiagonalCorrelators_'+ aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion)
                
                ### Eigenvalues 
                if os.path.isfile(myPlotLocation + 'Eigenvalues_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'Eigenvalues_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion)
                
                ### Eigenvalues log-plots
                if os.path.isfile(myPlotLocation + 'Eigenvalues_' + aa +  '_%s_log'%str(bb) + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'Eigenvalues_' + aa +  '_%s_log'%str(bb) + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion)
                
                ### Histogram Eigenvlaues
                if os.path.isfile(myPlotLocation +'Histogram_Eigenvalues_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0) + reBin +  '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation +'Histogram_Eigenvalues_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0) + reBin +  '_v%s.pdf'%myVersion)
                
                ### Effective Mass Eigenvlaues
                if os.path.isfile(myPlotLocation + 'EffectiveMass_Eigenvalues_'+ aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'EffectiveMass_Eigenvalues_'+ aa + '_%s'%str(bb) + reBin + '_v%s.pdf'%myVersion)
                
                ### Fits Eigenvalues
                if os.path.isfile(myPlotLocation + 'Tmin_Fits_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0)  + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'Tmin_Fits_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0)  + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
                
                ### Chi^{2} Fits Eigenvalues 
                if os.path.isfile(myPlotLocation + 'Tmin_Chisqr_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0)  + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'Tmin_Chisqr_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0)  + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
                
                ### Zoom Chi^{2} Fits Eigenvalues 
                if os.path.isfile(myPlotLocation + 'Tmin_Chisqr_Zoom_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0)  + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'Tmin_Chisqr_Zoom_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0)  + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
                    
                 ### Total Chi^{2} Fits Eigenvalues 
                if os.path.isfile(myPlotLocation + 'Tmin_TotalChisqr_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0)  + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'Tmin_TotalChisqr_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0)  + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
                    
                ### Delta Chi^{2} Fits
                if os.path.isfile(myPlotLocation + 'Tmin_DeltaChisqr_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0)  + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'Tmin_DeltaChisqr_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0)  + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
                
            ### All Correlators together log-plot
            if os.path.isfile(myPlotLocation + 'ALLDiagonalCorrelators_'+ aa + '_log'+ reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'ALLDiagonalCorrelators_'+ aa + '_log'+ reBin + '_v%s.pdf'%myVersion)

            ### Effective Masses all diagonal correlators together
            if os.path.isfile(myPlotLocation + 'EffectiveMass_ALLDiagonalCorrelators_'+ aa + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'EffectiveMass_ALLDiagonalCorrelators_'+ aa + reBin + '_v%s.pdf'%myVersion)
            
            ### All eigenvalues together log-plots
            if os.path.isfile(myPlotLocation + 'ALLEigenvalues_'+ aa + '_log'+ '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'ALLEigenvalues_'+ aa + '_log'+ '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion)
        
            ### Effective Masses all eigenvalues together
            if os.path.isfile(myPlotLocation + 'EffectiveMass_ALLEigenvalues_'+ aa + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion):
                x.append(myPlotLocation + 'EffectiveMass_ALLEigenvalues_'+ aa + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion)
            
            ### If the operators analysis was performed, then it will look for the plots
            if 'Operators_Analysis' in list(myMatrixCorrelatorData[aa].keys()):
                
                ### What type of operators analysis is in
                if any('Ops_chosen' in the_keys for the_keys in list(myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys())):
                    
                    the_list_of_chosen_ops = list(filter(lambda x: 'Ops_chosen' in x, myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys()))
                    
                    ### Loop over all the operators analysis
                    for the_op_item in the_list_of_chosen_ops:
                        the_len_data = len(myMatrixCorrelatorData[aa+'/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Mean'%myT0])
                        ### Loop over the eigenvalues of this analysis
                        for bb in range(the_len_data):
                            if os.path.isfile(myPlotLocation + 'EffectiveMass_Eigenvalues_'+ the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion):
                                x.append(myPlotLocation + 'EffectiveMass_Eigenvalues_'+ the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion)
                    
                    ### All eigenvalues together in one plot
                    if os.path.isfile(myPlotLocation + 'EffectiveMass_ALLEigenvalues_'+ the_op_item +'_' + aa + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion):
                            x.append(myPlotLocation + 'EffectiveMass_ALLEigenvalues_'+ the_op_item +'_' + aa + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion)
                
            if any('Add_Op' in the_keys for the_keys in list(myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys())):
                
                
                the_list_of_chosen_ops = list(filter(lambda x: 'Add_Op' in x, myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys()))
                
                ### Loop over all the operators analysis
                for the_op_item in the_list_of_chosen_ops:
                    the_len_data = len(myMatrixCorrelatorData[aa+'/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Mean'%myT0])
                    
                    ### Loop over the eignvalues
                    for bb in range(the_len_data):
                        if os.path.isfile(myPlotLocation + 'EffectiveMass_Eigenvalues_'+ the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion):
                            x.append(myPlotLocation + 'EffectiveMass_Eigenvalues_'+ the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion)
                
                ### All the eigenvalues together in one plot
                if os.path.isfile(myPlotLocation + 'EffectiveMass_ALLEigenvalues_'+ the_op_item +'_' + aa + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion):
                        x.append(myPlotLocation + 'EffectiveMass_ALLEigenvalues_'+ the_op_item +'_' + aa + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion)
                    
            if any('Remove_Op' in the_keys for the_keys in list(myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys())):
                the_list_of_chosen_ops = list(filter(lambda x: 'Remove_Op' in x, myMatrixCorrelatorData[aa+'/Operators_Analysis'].keys()))
                
                ### Loop over all the operators of the analysis
                for the_op_item in the_list_of_chosen_ops:
                    the_len_data = len(myMatrixCorrelatorData[aa+'/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Mean'%myT0])
                    
                    ### Loop over the eigenvalues of this selection
                    for bb in range(the_len_data):
                        if os.path.isfile(myPlotLocation + 'EffectiveMass_Eigenvalues_'+ the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion):
                            x.append(myPlotLocation + 'EffectiveMass_Eigenvalues_'+ the_op_item +'_' + aa + '_%s'%str(bb) + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion)
                
                ### All eigenvalues together in one plot
                if os.path.isfile(myPlotLocation + 'EffectiveMass_ALLEigenvalues_'+ the_op_item +'_' + aa + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion):
                    x.append(myPlotLocation + 'EffectiveMass_ALLEigenvalues_'+ the_op_item +'_' + aa + '_t0_%s'%str(myT0) + reBin + '_v%s.pdf'%myVersion)

            merger = PdfMerger()
            for pdf in x:
                merger.append(open(pdf, 'rb'))
            
            with open(myPlotLocation + myEns + "_%s"%aa + "_t0%s"%str(myT0) + reBin + "_%s"%myTypeRs + "_v%s.pdf"%myVersion, "wb") as fout:
                merger.write(fout)

        print('Now all the plots are in one file for each irrep')
        
    myMatrixCorrelatorData.close()
    
    
elif myWhichCorrelator=='mr':
    myRatioCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
    myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOU_PLOTS_DIRECTORY$/Plots/%s/Matrices_Ratios/'%myEns +  '%s/'%myResamplingScheme)

    if plotCorrs: 
        pcorr.PlotMultiHadronCorrelators(myRatioCorrelatorData, myTypeRs, myVersion, myT0, myPlotLocation, reBin)

    if plotEffMass: 
        peff.PlotMultiHadronsEffectiveMasses(myRatioCorrelatorData, myResamplingScheme, myVersion, myT0, myPlotLocation,reBin)
        
    if plotFits:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myDataLocation + 'Fits_Ratios/')
        myFitCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        pfit.PlotMultiHadronsFits(myFitCorrelator, myTypeCorrelation, myNrExponentials,  multiTMinsFitPlots, myT0, myVersion, myPlotLocation, reBin)
        
        myFitCorrelator.close()

    if joinPlots:
        print('Now all the plots are in one file')
        
    myRatioCorrelatorData.close()

### --------- PRINTS WHERE IT IS SAVED ---------
        
print('-'*(len(myPlotLocation)+1))
print('Correlator Plots saved : \n' + myPlotLocation)
print('_'*(len(myPlotLocation)+1))

