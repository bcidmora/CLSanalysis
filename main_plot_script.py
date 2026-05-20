###     PLOTTING SCRIPTS
import plot_correlators_script as pcorr
import plot_effective_masses_script as peff
import plot_fits_script as pfit
import plot_fitted_eff_masses as pfem

### LAYOUT AND EXTRA FUNCTIONS
import set_of_layout_functions as vfl
import set_of_analysis_functions as vfa
import set_of_plot_functions as vfp

### EXTRA LIBRARIES
import sys
import os
import h5py
from PyPDF2 import PdfMerger

import ensembles as ed


### ------- WHAT IT IS DONE --------

myArgs = vfp.parse_args()
myRuns = vfp.WhichRuns(myArgs, ed.ensembles)


### ------ MAIN VARIABLES ---------
reBin = f"_bin{myRuns.rb}" if myRuns.rebin else ""

### This is the information about the resampling
if myRuns.rs_type=='jk':
    myResamplingScheme='Jackknife'
elif myRuns.rs_type=='bt':
    myResamplingScheme='Bootstrap' 
    

### Information from the ensembles dictionary
myDataLocation = vfl.DIRECTORY_EXISTS(f'{ed.outputLocation}{myRuns.ensemble}/')

### -------- PRINTING INFO OF ENSEMBLE ---------

vfl.INFO_PRINTING(myRuns.correlator, myRuns.ensemble)

### ------------ START ----------------

if myRuns.correlator =='s':
    if not myRuns.ib_corr:
        myVersion =  f'{myRuns.ensemble}_singles_fwd' 
        myArchivoPre = ed.ensembles[myRuns.ensemble]
    else:
        myVersion =  f'{myRuns.ensemble}_omega' 
        myArchivoPre = ed.ensembles[myRuns.ensemble]['ib']        
    
    ### Original list of irreps
    myArchivo = h5py.File(myArchivoPre['fs'], 'r')
    myIrreps = list(myArchivo.keys())
    myArchivo.close()

    ### Correlators data
    mySingleCorrelatorData = h5py.File(f'{myDataLocation}/Single_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5','r')
    
    ### Directory where the plots will be saved
    myPlotLocation = vfl.DIRECTORY_EXISTS(f'{ed.location}/Plots/{myRuns.ensemble}/SingleHadrons/{myResamplingScheme}/'
    
    if myRuns.corrs:    
        pcorr.PlotSingleHadronCorrelators(mySingleCorrelatorData, myRuns.rs_type, myVersion, myPlotLocation, reBin, nr_irreps = myRuns.the_irreps.nr_irreps, first_irrep = myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep)
    
    if myRuns.effmass: 
        peff.PlotSingleHadronsEffectiveMasses(mySingleCorrelatorData, myResamplingScheme, myVersion, myPlotLocation, reBin, nr_irreps=myRuns.the_irreps.nr_irreps, first_irrep=myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep)
    
    if myRuns.fits:
        mySingleTmins = myArchivoPre[f'singleTMinResults-{myRuns.fit_param.type_fit}exp']
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myDataLocation}/Fits_SingleHadrons/')
        myFitCorrelator =  h5py.File(f'{myFitsLocation}Single_correlators_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'a')
        
        pfit.PlotSingleHadronsFits(myFitCorrelator, myRuns.fit_param.type_correlation, myRuns.fit_param.type_fit, mySingleTmins, myVersion, myPlotLocation, reBin, myIrreps, first_irrep = myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep, zoom_fit = myRuns.zoom_fit, chi_plots = myRuns.plot_chi, total_chi = myRuns.total_chi_plot, delta_chi = myRuns.delta_chi_plot, all_fits_comparison = myRuns.all_fits, nr_irreps=myRuns.the_irreps.nr_irreps)
        
        myFitCorrelator.close()
    
    if myRuns.fitmass:        
        mySingleTmins = myArchivoPre[f'singleTMinResults-{myRuns.fit_param.type_fit}exp']
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myDataLocation}Fits_SingleHadrons/')
        myFitCorrelator =  h5py.File(f'{myFitsLocation}Single_correlators_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'a')
        
        pfem.PlotSingleHadronsEffectiveMassesFits(myFitCorrelator, mySingleCorrelatorData, myResamplingScheme, myRuns.fit_param.type_correlation, myRuns.fit_type, mySingleTmins, myVersion, myPlotLocation, reBin, myIrreps, first_irrep=myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep)
    
    ### Puts all the plots in one PDF file. It checks if the file exists first
    if myRuns.join:
        ### Loop over all the irreps in this ensemble
        irreps = list(mySingleCorrelatorData.keys())
        all_ef_mass_x = []
        for aa in irreps:
            ### Corrs
            the_corr_plot = f'{myPlotLocation}Correlator_{aa[:4]}_{aa[-1]}{reBin}_{myVersion}.pdf'
            
            ### Corrs log-plots
            the_corr_log_plot = f'{myPlotLocation}Correlator_{aa[:4]}_{aa[-1]}_log{reBin}_{myVersion}.pdf'
            
            ### Histogram Corrs
            the_hist_plot = f'{myPlotLocation}Histogram_correlators_{aa[:4]}_{aa[-1]}{reBin}_{myVersion}.pdf'
            
            ### Effective Masses Corrs
            the_eff_mass_plot = f'{myPlotLocation}EffectiveMass_{aa[:4]}_{aa[-1]}{reBin}_{myVersion}.pdf'
            
            ### Fits Corrs
            the_tmin_plot = f'{myPlotLocation}Tmin_Fits_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
            
            ### Zoom Fits Corrs
            the_tmin_zoom_plot = f'{myPlotLocation}Tmin_Fits_Zoom_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
            
            ### Chi^{2} Fits Corrs
            the_chi_plot = f'{myPlotLocation}Tmin_Chisqr_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
            
            ### Total Chi^{2} Fits Corrs
            the_total_chi_plot = f'{myPlotLocation}Tmin_TotalChisqr_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
            
            ### Delta Chi^{2} Fits Corrs
            the_delta_chi_plot = f'{myPlotLocation}Tmin_DeltaChisqr_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
            
            ### All fits in one plot
            the_diff_fits_plot = f'{myPlotLocation}Different_Fits_{aa[:4]}_{aa[-1]}_{reBin}_{myVersion}.pdf'
            
            ### Fitted effective masses
            the_fitted_mass_plot = f'{myPlotLocation}Fitted_Effective_Masses_{aa[:4]}_{aa[-1]}_{myRuns.fit_param.type_fit}exp_{reBin}_{myVersion}.pdf'
            
            # the_plots = [the_corr_plot, the_corr_log_plot, the_hist_plot, the_eff_mass_plot, the_tmin_plot, the_tmin_zoom_plot, the_chi_plot, the_total_chi_plot, the_delta_chi_plot, the_diff_fits_plot, the_fitted_mass_plot]
            
            the_plots = [the_tmin_plot, the_chi_plot, the_fitted_mass_plot]
        
            x=[]
            for item in the_plots:
                if os.path.isfile(item): x.append(item)
                if 'Fitted_Effective_Masses' in item: all_ef_mass_x.append(item)
                
            merger = PdfMerger()
            for pdf in x:
                merger.append(open(pdf, 'rb'))
            
            with open(f'{myPlotLocation}{myRuns.ensemble}_{aa}{reBin}_{myRuns.rs_type}_{myVersion}.pdf', "wb") as fout:
                merger.write(fout)
        
        merger_eff = PdfMerger()
        for pdf_eff in all_ef_mass_x:
            merger_eff.append(open(pdf_eff, 'rb'))
        with open(f'{myPlotLocation}{myRuns.ensemble}_All_Fitted_Masses{reBin}_{myRuns.rs_type}_{myVersion}.pdf', "wb") as fout:
            merger_eff.write(fout)
        print('Now all the plots are in one file for each irrep')
    
    mySingleCorrelatorData.close()
    
elif myRuns.correlator=='m':
    myIsospin = myRuns.isospin
    myChosenIsospin = ed.ensembles[myRuns.ensemble][myIsospin]['iso_tag']    
    myArchivo = h5py.File(ed.ensembles[myRuns.ensemble][myIsospin]['fm'], 'r')
    myIrreps = list(myArchivo.keys())
    myArchivo.close()
    print(myRuns.all_corr)
    
    if myRuns.ops_flag:
        myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_reduced' 
    else:
        myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_test' 
    
    myMatrixCorrelatorData = h5py.File(f'{myDataLocation}Matrix_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5','r')
    
    ### Directory where the plots will be saved
    myPlotLocation = vfl.DIRECTORY_EXISTS(f'{ed.location}/Plots/{myRuns.ensemble}/Matrices/{myResamplingScheme}/')
    
     ### Plots the correlators. Look at the booleans
    if myRuns.corrs: 
        pcorr.PlotMultiHadronCorrelators(myMatrixCorrelatorData, myChosenIsospin, myRuns.rs_type, myVersion, myRuns.fit_param.t0, myPlotLocation, reBin, nr_irreps = myRuns.the_irreps.nr_irreps, first_irrep = myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep, diag_corrs = myRuns.diag_flag, all_corr = myRuns.all_corr)
    
    ### Plots the effective masses of the eigenvalues from the GEVP and/or the operators analysis
    if myRuns.effmass: 
        peff.PlotMultiHadronsEffectiveMasses(myMatrixCorrelatorData, myChosenIsospin, myResamplingScheme, myVersion, myRuns.fit_param.t0, myPlotLocation, reBin, nr_irreps=myRuns.the_irreps.nr_irreps, first_irrep=myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep, diag_corrs= myRuns.diag_flag,all_corr = myRuns.all_corr)

        
    if myRuns.fits:     
        multiTMinsFitPlots = ed.ensembles[myRuns.ensemble][myIsospin]['multiTMinResults']
        myFitsLocation = vfl.DIRECTORY_EXISTS(myDataLocation + 'Fits_Matrices/')
        myFitCorrelator = h5py.File(f'{myFitsLocation}Matrix_correlators_{myChosenIsospin}_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'a')
                             
        pfit.PlotMultiHadronsFits(myFitCorrelator, myChosenIsospin, myRuns.fit_param.type_correlation, myRuns.fit_param.type_fit, myRuns.rs_type, multiTMinsFitPlots, myRuns.fit_param.t0, myVersion, myPlotLocation, reBin, myIrreps, gevp=myRuns.gevp_flag, zoom_fit=myRuns.zoom_fit, chi_plots=myRuns.plot_chi, first_irrep=myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep, total_chi=myRuns.total_chi_plot, delta_chi=myRuns.delta_chi_plot, ops_analysis=myRuns.ops_flag, ops_analysis_method=myOperatorsMethod)
        
        myFitCorrelator.close()
    
    ### Plot the fitted effective masses of the eigenvalues
    if myRuns.fitmass:
        multiTMinsFitPlots = ed.ensembles[myRuns.ensemble][myIsospin]['multiTMinResults']
        myFitsLocation = vf.DIRECTORY_EXISTS(f'{myDataLocation}Fits_Matrices/')
        myFitCorrelator = h5py.File(f'{myFitsLocation}Matrix_correlators_{myChosenIsospin}_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'a')
        
        pfem.PlotMultiHadronsEffectiveMassesFits(myFitCorrelator, myMatrixCorrelatorData, myChosenIsospin, myResamplingScheme, myRuns.fit_param.type_correlation, myRuns.fit_param.type_fit, multiTMinsFitPlots, myRuns.fit_param.t0, myVersion, myPlotLocation, reBin, myIrreps) 
        myFitCorrelator.close()

    ### This part puts all the plots together in one pdf file. 
    if myRuns.join:
        ### Loop over all the irreps in this ensemble
        irreps = list(myMatrixCorrelatorData.keys())
        for aa in irreps:
            ##3 Loop over all the operators in this ensemble
            ops = list(myMatrixCorrelatorData[aa+'/Operators'])
            x=[]
            for bb in range(len(ops)):
                
                ## Diagonal Corrs
                the_corr_plot = f'{myPlotLocation}DiagonalCorrelator_{myChosenIsospin}_{aa}_{bb}{reBin}_{myVersion}.pdf'
                
                ### Diagonal corrs log-plot
                the_corr_log_plot = f'{myPlotLocation}DiagonalCorrelator_{myChosenIsospin}_{aa}_{bb}_log{reBin}_{myVersion}.pdf'
                
                ### Histogram Corrs
                the_hist_plot = f'{myPlotLocation}Histogram_DiagCorrelator_{myChosenIsospin}_{aa}_{bb}{reBin}_{myVersion}.pdf'
                
                ### Effective Mass diagonal corrs
                the_eff_mass_plot = f'{myPlotLocation}EffectiveMass_DiagonalCorrelators_{myChosenIsospin}_{aa}_{bb}{reBin}_{myVersion}.pdf'
                
                ### Eigenvalues 
                the_eigens_plot = f'{myPlotLocation}Eigenvalues_{myChosenIsospin}_{aa}_{bb}_t0_{myRuns.fit_param.t0}{reBin}_{myVersion}.pdf'
                
                ### Eigenvalues log-plots
                the_eigens_log_plot = f'{myPlotLocation}Eigenvalues_{myChosenIsospin}_{aa}_{bb}_log_t0_{myRuns.fit_param.t0}{reBin}_{myVersion}.pdf'
                
                ### Histogram Eigenvlaues
                the_eigens_hist_plot = f'{myPlotLocation}Histogram_Eigenvalues_{myChosenIsospin}_{aa}_{bb}_t0_{myRuns.fit_param.t0}{reBin}_{myVersion}.pdf'
                
                ### Effective Mass Eigenvlaues
                the_eigens_eff_mass_plot = f'{myPlotLocation}EffectiveMass_Eigenvalues_{myChosenIsospin}_{aa}_{bb}_t0_{myRuns.fit_param.t0}{reBin}_{myVersion}.pdf'
                
                ### Fits Eigenvalues
                the_tmin_plot = f'{myPlotLocation}Tmin_Fits_{myChosenIsospin}_{aa}_{bb}_t0_{myRuns.fit_param.t0}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
                
                ### Chi^{2} Fits Eigenvalues 
                the_chi_plot = f'{myPlotLocation}Tmin_Chisqr_{myChosenIsospin}_{aa}_{bb}_t0_{myRuns.fit_param.t0}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
                
                ### Total Chi^{2} Fits Eigenvalues 
                the_total_chi_plot = f'{myPlotLocation}Tmin_TotalChisqr_{myChosenIsospin}_{aa}_{bb}_t0_{myRuns.fit_param.t0}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
                
                ### Delta Chi^{2} Fits Corrs
                the_delta_chi_plot = f'{myPlotLocation}Tmin_DeltaChisqr_{myChosenIsospin}_{aa}_{bb}_t0_{myRuns.fit_param.t0}_{myRuns.fit_param.type_fit}exp{reBin}_{myVersion}.pdf'
                
                ### All fits in one plot
                # the_diff_fits_plot = f'{myPlotLocation}Different_Fits_{aa[:4]}_{aa[-1]}_{reBin}_{myVersion}.pdf'
                
                ### Fitted Effective Masses of Eigenvalues
                the_fitted_mass_plot = f'{myPlotLocation}Fitted_Effective_Masses_{myChosenIsospin}_{aa}_{bb}_t0_{myRuns.fit_param.t0}{reBin}_{myVersion}.pdf'
                
                the_plots = [the_corr_plot, the_corr_log_plot, the_hist_plot, the_eff_mass_plot, the_eigens_plot, the_eigens_log_plot, the_eigens_hist_plot, the_eigens_eff_mass_plot, the_tmin_plot, the_chi_plot, the_total_chi_plot, the_delta_chi_plot, the_fitted_mass_plot]
                
                for item in the_plots:
                    if os.path.isfile(item): x.append(item)
                  
            ### All Correlators together log-plot
            the_all_diag_corr_log_plot = f'{myPlotLocation}ALLDiagonalCorrelators_{myChosenIsospin}_{aa}_log{reBin}_{myVersion}.pdf'
            
            ### Effective Masses all diagonal correlators together
            the_all_diag_corr_effmass_plot = f'{myPlotLocation}EffectiveMass_ALLDiagonalCorrelators_{myChosenIsospin}_{aa}{reBin}_{myVersion}.pdf'

            ### All eigenvalues together log-plots
            the_all_eigens_log_plot = f'{myPlotLocation}ALLEigenvalues_{myChosenIsospin}_{aa}_log_t0_{myRuns.fit_param.t0}{reBin}_{myVersion}.pdf'
        
            ### Effective Masses all eigenvalues together
            the_all_eigens_eff_mass_plot = f'{myPlotLocation}EffectiveMass_ALLEigenvalues_{myChosenIsospin}_{aa}_t0_{myRuns.fit_param.t0}{reBin}_{myVersion}.pdf'
            
            the_all_plots = [the_all_diag_corr_log_plot, the_all_diag_corr_effmass_plot, the_all_eigens_log_plot, the_all_eigens_eff_mass_plot]
            
            for item in the_all_plots: 
                if os.path.isfile(item): x.append(item)
                        
            merger = PdfMerger()
            for pdf in x:
                merger.append(open(pdf, 'rb'))
            
            with open(f'{myPlotLocation}{myRuns.ensemble}_{myChosenIsospin}_{aa}_t0{myRuns.fit_param.t0}{reBin}_{myRuns.rs_type}_{myVersion}.pdf', "wb") as fout:
                merger.write(fout)
        
        print('Now all the plots are in one file for each irrep')
        
    myMatrixCorrelatorData.close()
    
    
elif myRuns.correlator=='mr':
    myIsospin = myRuns.isospin
    myNonInteractingLevels = ed.ensembles[myRuns.ensemble][myRuns.isospin]['nonInteractingLevels']
    myChosenIsospin = ed.ensembles[myRuns.ensemble][myIsospin]['iso_tag']
    
    if myRuns.ops_flag:
        myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_reduced' 
    else:
        myVersion =  f'{myRuns.ensemble}_{myChosenIsospin}_test'
        
    myRatioCorrelatorData = h5py.File(f'{myDataLocation}Ratio_matrix_correlators_{myRuns.rs_type}{reBin}_{myVersion}.h5','r')
    
    myPlotLocation = vfl.DIRECTORY_EXISTS(f'{ed.location}/Plots/{myRuns.ensemble}/Matrices_Ratios/{myResamplingScheme}/')

    if myRuns.corrs: 
        pcorr.PlotRatioHadronCorrelators(myRatioCorrelatorData, myChosenIsospin, myRuns.rs_type, myVersion, myRuns.fit_param.t0, myPlotLocation, reBin, first_irrep = myRuns.the_irreps.first_irrep, last_irrep = myRuns.the_irreps.last_irrep)

    if myRuns.effmass: 
        peff.PlotRatioHadronsEffectiveMasses(myRatioCorrelatorData, myChosenIsospin, myRuns.rs_type, myVersion, myRuns.fit_param.t0, myPlotLocation, reBin)
    
    ### NOT WORKED OUT YET
    if myRuns.fits:  
        multiTMinsFitPlots = ed.ensembles[myRuns.ensemble][myIsospin]['multiTMinResultsRatios']
        myFitsLocation = vfl.DIRECTORY_EXISTS(f'{myDataLocation}Fits_Ratios/')
        myFitCorrelator = h5py.File(f'{myFitsLocation}Ratio_matrix_correlators_{myChosenIsospin}_{myRuns.rs_type}{reBin}_fits_{myVersion}.h5', 'a')
        
        pfit.PlotRatioMultiHadronsFits(myFitCorrelator, myChosenIsospin, myRuns.fit_param.type_correlation, myRuns.fit_param.type_fit, myRuns.rs_type, multiTMinsFitPlots, myRuns.fit_param.t0, myVersion, myPlotLocation, reBin, chi_plots = myRuns.plot_chi )
        
        myFitCorrelator.close()

    if myRuns.join:
        print('Now all the plots are in one file')
        
    myRatioCorrelatorData.close()

### --------- PRINTS WHERE IT IS SAVED ---------
        
print('-'*(len(myPlotLocation)+1))
print('Correlator Plots saved : \n' + myPlotLocation)
print('_'*(len(myPlotLocation)+1))

