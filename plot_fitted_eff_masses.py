import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_functions as vf


def PlotSingleHadronsEffectiveMassesFits(the_single_fit_data, the_single_correlator_data, the_rs_scheme, the_type_fit, the_nr_exps, the_tmins, the_version, the_location, the_rebin, the_irreps, **kwargs):    
    
    the_eff_mass_color = '#5d83d5'
    the_fit_color = '#b90f22'
    
    ### These are the irreps in this file
    s_irreps = list(the_single_fit_data.keys())
    
    ### If not all the irreps are wanted t be plotted
    if kwargs.get('nr_irreps')!=None:
        the_first_irrep = 0
        the_last_irrep = int(kwargs.get('nr_irreps'))
    else:
        the_first_irrep = 0
        the_last_irrep = len(s_irreps)
    ### This one checks for an irrep in particular
    if kwargs.get('first_irrep')!=None and kwargs.get('last_irrep')!=None:
        the_first_irrep = int(kwargs.get('first_irrep'))
        the_last_irrep = int(kwargs.get('last_irrep'))
    elif kwargs.get('last_irrep')!=None and kwargs.get('first_irrep')==None:
        the_first_irrep = 0
        the_last_irrep = int(kwargs.get('last_irrep'))
    elif kwargs.get('first_irrep')!=None and kwargs.get('last_irrep')==None:
        the_first_irrep = int(kwargs.get('first_irrep'))
        the_last_irrep = len(s_irreps)
    
    s_irreps = s_irreps[the_first_irrep:the_last_irrep]
    
    
    ### Loop over the irreps found in the fits file.
    for the_irrep in s_irreps:
        
        print("---------------------------------------------------------------------------")
        print("Irrep: ", the_irrep)
        
        ### Central values of the correlators
        the_mean_corr = np.array(the_single_correlator_data[the_irrep + '/Effective_masses/Mean'])
        
        ### Statistical errors of the correlators
        the_sigmas_corr = np.array(the_single_correlator_data[the_irrep + '/Effective_masses/Sigmas'])
        
        ### Time extent
        the_nt_corr = np.array(the_single_correlator_data[the_irrep + '/Time_slices'])
        the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
        the_nt_ticks = np.arange(the_nt_corr[0]+1, the_nt_corr[-1], int(len(the_nt_corr)/5))
        
        ### The SH operator that appears in the plot
        the_op = list(the_single_correlator_data[the_irrep+'/Operators'])[0]
        OperatorNamePlot = vf.OPERATORS_SH(the_op.decode('utf-8'))
        
        OperatorNamePlot = vf.PLOT_SINGLE_HADRON_NAMES(OperatorNamePlot)
        
        da_irrep = vf.IrrepInfo(the_irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot
        
        ### This is the info of the fits for this irrep
        dis_set =  np.array(the_single_fit_data[the_irrep + '/%sexp'%the_nr_exps + '/Tmin/%s/Mean'%the_type_fit])
        
        ### This is the range of min time slices which the fits were performed
        the_nt_fit = [int(x) for x in dis_set[0]]
        
        ### Max time slice used for the fit
        the_nt_max = int(dis_set[1][0])
        
        ### This is the tmin chosen for this fit. It can be changed in "file_ens.py"
        the_chosen_tmin = the_nt_fit.index(the_tmins[the_irreps.index(the_irrep)])
        
        ### Mean values of the fits
        the_fit_data = dis_set[2]
        
        ### Statistical errors of those central values
        the_fit_sigmas = dis_set[3]
        
        ### The Chi^2
        the_chi_corr = dis_set[4]    
        
        ### The title of the Plot
        the_title = OperatorNamePlot + ' (%s): '%MomentumIrrep + r' $\to$ %s '%NameIrrepPlot + '(%s-exp)'%the_nr_exps
       
        ### This is just to write the errors properly in the plot
        the_mean_fit_string = str('{:.5f}'.format(np.round(the_fit_data[the_chosen_tmin], 5)))
        the_error_string = vf.WRITTING_ERRORS_PLOTS(the_fit_sigmas[the_chosen_tmin],5)
        the_sigmas_fit_string = the_error_string[0]

        if the_error_string[1]==False:
            the_mean_fit_string = str(f'{np.round(the_fit_data[the_chosen_tmin], the_error_string[2]):.{the_error_string[2]}f}')   
        
        print("Plotting Fitted Effective Masses...")
        
        fit_fig = plt.figure()     
        ### THIS IS THE LABEL WITH TMIN AND TMAX
        # the_label = r'$t_{\mathrm{min}} = %s$'%str(int(the_nt_fit[the_chosen_tmin])) + '\n' +  r'$t_{\mathrm{max}} = %s$'%str(int(the_nt_max))  + '\n' + r'$\chi^{2}/\mathrm{d.o.f} = %s$'%np.round(the_chi_corr[the_chosen_tmin],3) + '\n' + r'$E_{\mathrm{fit}} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
        
        the_label = r'$\chi^{2}/\mathrm{d.o.f} = %s$'%np.round(the_chi_corr[the_chosen_tmin],3) + '\n' + r'$E_{\mathrm{fit}} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
        
        
        vf.PLOT_FITTED_EFF_MASSES(the_nt, the_mean_corr, the_sigmas_corr, the_fit_data, the_fit_sigmas, the_chosen_tmin, the_rs_scheme + ' data', the_label, the_title, the_nt_ticks, the_eff_mass_color, the_fit_color)
        
        fit_fig.savefig(the_location + 'Fitted_Effective_Masses_' + the_irrep[:4] +'_%s'%the_irrep[-1] + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version, bbox_inches='tight')
        
        
        
def PlotMultiHadronsEffectiveMassesFits(the_multi_hadrons_fit_data, the_matrix_correlator_data, the_quantum_number, the_rs_scheme, the_type_fit, the_nr_exps, the_tmins, the_t0, the_version, the_location, the_rebin, the_irreps, **kwargs):  
    
    m_irreps = list(the_multi_hadrons_fit_data.keys())
    
    the_eff_mass_color = '#5d83d5'
    the_fit_color = '#b90f22'

    the_ops_analysis_flag = kwargs.get('ops_analysis')
    
    ### If not all the irreps are wanted t be plotted
    if kwargs.get('nr_irreps')!=None:
        the_first_irrep = 0
        the_last_irrep = int(kwargs.get('nr_irreps'))
    else:
        the_first_irrep = 0
        the_last_irrep = len(m_irreps)
    ### This one checks for an irrep in particular
    if kwargs.get('first_irrep')!=None and kwargs.get('last_irrep')!=None:
        the_first_irrep = int(kwargs.get('first_irrep'))
        the_last_irrep = int(kwargs.get('last_irrep'))
    elif kwargs.get('last_irrep')!=None and kwargs.get('first_irrep')==None:
        the_first_irrep = 0
        the_last_irrep = int(kwargs.get('last_irrep'))
    elif kwargs.get('first_irrep')!=None and kwargs.get('last_irrep')==None:
        the_first_irrep = int(kwargs.get('first_irrep'))
        the_last_irrep = len(m_irreps)
    
    m_irreps = m_irreps[the_first_irrep:the_last_irrep]
    

    ### Loop over the irreps of this file
    for the_irrep in m_irreps:
        
        print("---------------------------------------------------------------------------")
        print("Irrep: ", the_irrep)
        
        ### Searching if the fit was done and the gevp plots must be included
        if '%sexp'%the_nr_exps in list(the_multi_hadrons_fit_data[the_irrep].keys()): 
            
            ### Retrieving the data
            the_data_fit = the_multi_hadrons_fit_data[the_irrep + '/%sexp/'%the_nr_exps + 't0_%s/Tmin/'%the_t0 + '%s/Mean'%the_type_fit]
            
            the_data = np.array(the_matrix_correlator_data[the_irrep + '/GEVP/t0_%s/Effective_masses/Mean'%the_t0])
            
            the_data_sigmas = np.array(the_matrix_correlator_data[the_irrep + '/GEVP/t0_%s/Effective_masses/Sigmas'%the_t0])
            
            ### Loop over the eigenvalues of this irrep
            for bb in range(len(list(the_data_fit.keys()))):
                
                ### The effective masses of the eigenvalues
                the_mean_corr = the_data[bb]
                
                ### The sigmas of these eff masses
                the_sigmas_corr = the_data_sigmas[bb]
                
                ### The time slices range
                the_nt_corr = np.array(the_matrix_correlator_data[the_irrep + '/Time_slices'])
                
                ### This is the time slices shifted for the plot
                the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
                
                ### bb-th Eigenvalue
                dis_set = np.array(the_data_fit.get('lambda_%s'%bb))
                
                ### The central values of the diagonalized correlator
                the_fit_data = dis_set[2]
                
                ### The statistical errors
                the_fit_sigmas = dis_set[3]
                
                ### The chi^{2} of the fit
                the_chi_corr = dis_set[4]
                
                ### The statistical error of the chi^{2}
                the_chi_sigmas = dis_set[5]
                
                ## The minimum time slices that the fit was performed
                the_nt_fit = [int(x) for x in dis_set[0]]
                
                ### These are the ticks that appear in the plot
                the_nt_ticks = np.arange(int(the_nt_corr[-1]/5), the_nt_corr[-1], int(the_nt_corr[-1]/5))

                ### Information about the irre
                da_irrep = vf.IrrepInfo(the_irrep)
                MomentumIrrep = da_irrep.TotalMomPlot
                NameIrrepPlot = da_irrep.NamePlot
                
                ### Getting the position of the chosen tmin for the plots
                the_chosen_tmin = the_nt_fit.index(the_tmins[the_irreps.index(the_irrep)][bb])
                
                the_mean_fit_string = str('{:.5f}'.format(np.round(the_fit_data[the_chosen_tmin], 5)))
                the_error_string = vf.WRITTING_ERRORS_PLOTS(the_fit_sigmas[the_chosen_tmin],5)
                the_sigmas_fit_string = the_error_string[0]
                
                if the_error_string[1]==False:
                    the_mean_fit_string = str(f'{np.round(the_fit_data[the_chosen_tmin], the_error_string[2]):.{the_error_string[2]}f}')
                
                the_title = NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb)  + r' ($t_{0} = %s$)'%str(the_t0)

                ### THIS IS THE LABEL WITH TMIN AND TMAX
                # the_label = r'$t_{\mathrm{min}} = %s$'%str(int(the_nt_fit[the_chosen_tmin])) + '\n' + r'$t_{\mathrm{max}}$ = %s'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/\mathrm{d.o.f} = %s$'%np.round(the_chi_corr[the_chosen_tmin],3) + '\n' + r'$E_{\mathrm{fit}} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
                
                the_label = r'$\chi^{2}/\mathrm{d.o.f} = %s$'%np.round(the_chi_corr[the_chosen_tmin],3) + '\n' + r'$E_{\mathrm{fit}} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
                
                fit_fig = plt.figure()                
                vf.PLOT_FITTED_EFF_MASSES(the_nt, the_mean_corr, the_sigmas_corr, the_fit_data, the_fit_sigmas, the_chosen_tmin, the_rs_scheme + ' data', the_label, the_title, the_nt_ticks, the_eff_mass_color, the_fit_color)
                
                fit_fig.savefig(the_location + 'Fitted_Effective_Masses' + the_quantum_number + the_irrep + '_%s'%str(bb) + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
                    
        if 'Operators_Analysis' in list(the_multi_hadrons_fit_data[the_irrep].keys()) and the_ops_analysis_flag:
            if kwargs.get('ops_analysis_method')=='from_list':
                the_method = 'Ops_chosen_'
            elif kwargs.get('ops_analysis_method')=='adding':
                the_method = 'Add_Op_'
            elif kwargs.get('ops_analysis_method')=='removing':
                the_method = 'Remove_Op_'
            else: sys.exit("No method for the operators analysis chosen")
            
            the_list_of_chosen_ops = list(filter(lambda x: the_method in x, the_data[the_irrep+'/Operators_Analysis'].keys()))
            
                ### Loop over those elements
            for the_op_item in the_list_of_chosen_ops:
                
                the_data_fit = the_multi_hadrons_fit_data[the_irrep + '/Operators_Analysis/'+the_op_item+ '/t0_%s/Tmin/'%the_t0 + '%s/Mean'%the_type_fit]
                
                the_data = np.array(the_matrix_correlator_data[the_irrep + '/Operators_Analysis/' + the_op_item + '/t0_%s/Effective_masses/Mean'%the_t0])
            
                the_data_sigmas = np.array(the_matrix_correlator_data[the_irrep + '/Operators_Analysis/' + the_op_item +'/t0_%s/Effective_masses/Sigmas'%the_t0])
            
                ### This is the information of the selected operators (caption of the plot)
                the_chosen_ops_string = 'Ops:'
                for ss in the_op_item[11:]:
                    ### The opertaor written in a better form
                    the_op = the_op_list[int(ss)]
                    OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                    the_chosen_ops_string+=' '+OperatorNamePlot
                        
                        
                ### Loop over the eigenvalues of this irrep
                for bb in range(len(list(the_data.keys()))):
                    
                    the_mean_corr = the_data[bb]
                
                    the_sigmas_corr = the_data_sigmas[bb]
                
                    ### bb-th Eigenvalue
                    dis_set = np.array(the_data_fit.get('lambda_%s'%bb))
                    
                    ### The central values of the diagonalized correlator
                    the_fit_data = dis_set[2]
                    
                    ### The statistical errors
                    the_fit_sigmas = dis_set[3]
                    
                    ### The chi^{2} of the fit
                    the_chi_corr = dis_set[4]
                    
                    ### The tmins for the fits
                    the_nt_fit = [int(x) for x in dis_set[0]]
                    
                    ### The statistical error of the chi^{2}
                    the_chi_sigmas = dis_set[5]
                    
                    ### These are the ticks that appear in the plot
                    the_nt_ticks = np.arange(the_nt[0], the_nt[-1], 2)

                    ### Information about the irre
                    da_irrep = vf.IrrepInfo(the_irrep)
                    MomentumIrrep = da_irrep.TotalMomPlot
                    NameIrrepPlot = da_irrep.NamePlot
                    
                    the_chosen_tmin = the_nt_fit.index(the_tmins[the_irreps.index(the_irrep)][bb])
                    
                    the_mean_fit_string = str('{:.5f}'.format(np.round(the_fit_data[the_chosen_tmin], 5)))
                    the_error_string = vf.WRITTING_ERRORS_PLOTS(the_fit_sigmas[the_chosen_tmin],5)
                    the_sigmas_fit_string = the_error_string[0]
                        
                    if the_error_string[1]==False:
                        the_mean_fit_string = str(f'{np.round(the_fit_data[the_chosen_tmin], the_error_string[2]):.{the_error_string[2]}f}')
                
                
                    the_label = r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_chosen_tmin],3) + '\n' + r'$E_{fit} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
                    
                    the_title = NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb)  + r' ($t_{0} = %s$)'%str(the_t0)
                    
                    
                    print("Plotting Fitted Effective Masses...")
                    fit_fig = plt.figure() 
                    vf.PLOT_FITTED_EFF_MASSES(the_nt, the_mean_corr, the_sigmas_corr, the_fit_data, the_fit_sigmas, the_chosen_tmin, the_rs_scheme + ' data', the_label, the_title, the_nt_ticks, the_eff_mass_color, the_fit_color)
                    # plt.show()
                    fit_fig.savefig(the_location + 'Fitted_Effective_Masses' + the_quantum_number + the_op_item + '_' + the_irrep + '_%s'%str(bb) + '_t0_%s'%str(the_t0)  + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version)
        
                



### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------




if __name__=="__main__":
    
      myEns = str(sys.argv[1]).upper()
    myWhichCorrelator = str(sys.argv[2]).lower()
    myTypeRs = str(sys.argv[3]).lower()
    myRebinOn = str(sys.argv[4]).lower()
    
    myTypeFit = 'Correlated'
    myNrExponentials =  '1'
    myRb = 2
    myVersion = 'test'
    myT0 = 4 
    
    myDataLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH(SAME_THAN_CORRS_SCRIPT_OUTPUT)$/%s/'%myEns)
    
    if myEns == 'N451': from files_n451 import singleTMinsFitPlots, multiTMinsFitPlots
    elif myEns == 'N201': from files_n201 import singleTMinsFitPlots, multiTMinsFitPlots 
    elif myEns == 'D200': from files_d200 import singleTMinsFitPlots, multiTMinsFitPlots
    elif myEns == 'X451': from files_x451 import singleTMinsFitPlots, multiTMinsFitPlots
    
    vf.INFO_PRINTING(myWhichCorrelator, myEns)
    
    if myRebinOn=='rb': 
        reBin = '_bin'+str(myRb)
    else:
        reBin = '' 
        
    if myTypeRs=='jk':
        myResamplingScheme='Jackknife'
    elif myTypeRs=='bt':
        myResamplingScheme='Bootstrap'    

    if myWhichCorrelator=='s':
        mySingleCorrelatorData = h5py.File(myDataLocation + 'Fits_SingleHadrons/Single_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion,'r')
        
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/SingleHadrons/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotSingleHadronsFits(mySingleCorrelatorData, myTypeFit, myNrExponentials,  singleTMinsFitPlots, myVersion, myPlotLocation, reBin)
        
        mySingleCorrelatorData.close()
    
    elif myWhichCorrelator=='m':
        myMatrixCorrelatorData = h5py.File(myDataLocation + 'Fits_Matrices/Matrix_correlators_' + myTypeRs + reBin +'_fits_v%s.h5'%myVersion,'r')
        
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/Matrices/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotMultiHadronsFits(myMatrixCorrelatorData, myTypeFit, myNrExponentials, multiTMinsFitPlots, myT0, myVersion, myPlotLocation, reBin)
        
        myMatrixCorrelatorData.close()
        
    elif myWhichCorrelator=='mr':
        myRatioCorrelatorData = h5py.File(myDataLocation + 'Fits_Ratios/Matrix_correlators_ratios_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
        
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/Matrices_Ratios/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotMultiHadronsFits(myRatioCorrelatorData, myTypeFit, myNrExponentials, multiTMinsFitPlots, myT0, myVersion, myPlotLocation, reBin)
        
        myRatioCorrelatorData.close()
         

         
