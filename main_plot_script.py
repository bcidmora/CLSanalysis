import plot_correlators_script as pcorr
import plot_effective_masses_script as peff
import plot_fits_script as pfit
import set_of_functions as vf
import sys
import os
import h5py
from PyPDF2 import PdfMerger

### ------- WHAT IT IS DONE --------

plotCorrs = False
plotEffMass = False
plotFits = False
joinPlots = False


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


myDataLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH(SAME_THAN_CORRS_SCRIPT_OUTPUT)$/%s/'%myEns)

if plotFits or joinPlots: # SEE if this can be moved down
    if myEns == 'N451': from files_n451 import singleTMinsFitPlots, multiTMinsFitPlots, name, name1
    elif myEns == 'N201': from files_n201 import singleTMinsFitPlots, multiTMinsFitPlots, name, name1 
    elif myEns == 'D200': from files_d200 import singleTMinsFitPlots, multiTMinsFitPlots, name, name1
    elif myEns == 'X451': from files_x451 import singleTMinsFitPlots, multiTMinsFitPlots, name, name1
    
if myRebinOn=='rb': 
    reBin = '_bin'+str(myRb)
else: 
    reBin = ''


if myTypeRs=='jk':
    myResamplingScheme='Jackknife'
elif myTypeRs=='bt':
    myResamplingScheme='Bootstrap' 


### -------- PRINTING INFO OF ENSEMBLE ---------

vf.INFO_PRINTING(myWhichCorrelator, myEns)

### ------------ START ----------------

if myWhichCorrelator =='s':
    mySingleCorrelatorData = h5py.File(myDataLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion,'r')
    myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOU_PLOTS_DIRECTORY$/Plots/%s/SingleHadrons/'%myEns +  '%s/'%myResamplingScheme)
    
    if plotCorrs:        
        pcorr.PlotSingleHadronCorrelators(mySingleCorrelatorData, myTypeRs, myVersion, myPlotLocation, reBin)
    
    if plotEffMass: 
        peff.PlotSingleHadronsEffectiveMasses(mySingleCorrelatorData, myResamplingScheme, myVersion, myPlotLocation, reBin)
    
    if plotFits: 
        myFitsLocation = vf.DIRECTORY_EXISTS(myDataLocation + 'Fits_SingleHadrons/')
        myFitCorrelator =  h5py.File(myFitsLocation + 'Single_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        pfit.PlotSingleHadronsFits(myFitCorrelator, myTypeCorrelation, myNrExponentials,  singleTMinsFitPlots, myVersion, myPlotLocation, reBin)
        
        myFitCorrelator.close()
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
                
    mySingleCorrelatorData.close()
    
elif myWhichCorrelator=='m':        
    myMatrixCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
    myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOU_PLOTS_DIRECTORY$/Plots/%s/Matrices/'%myEns +  '%s/'%myResamplingScheme)
    
    if plotCorrs: 
        pcorr.PlotMultiHadronCorrelators(myMatrixCorrelatorData, myTypeRs, myVersion, myT0, myPlotLocation, reBin, nr_irreps=myNrIrreps, first_irrep=myFirstIrrep, last_irrep = myLastIrrep)
    
    if plotEffMass: 
        peff.PlotMultiHadronsEffectiveMasses(myMatrixCorrelatorData, myResamplingScheme, myVersion, myT0, myPlotLocation, reBin,nr_irreps=myNrIrreps, first_irrep=myFirstIrrep, last_irrep = myLastIrrep)
        
    if plotFits:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myDataLocation + 'Fits_Matrices/')
        myFitCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        pfit.PlotMultiHadronsFits(myFitCorrelator, myTypeCorrelation, myNrExponentials, multiTMinsFitPlots, myT0, myVersion, myPlotLocation, reBin)
        
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

        print('Now all the plots are in one file')
        
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
print('Correlator analysis saved : \n' + myPlotLocation)
print('_'*(len(myPlotLocation)+1))

