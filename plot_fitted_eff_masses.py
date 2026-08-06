import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_plot_functions as vfp

def PlotSingleHadronsEffectiveMassesFits(the_single_fit_data, the_single_correlator_data, the_rs_scheme, the_type_fit, the_nr_exps, the_tmins, the_version, the_location, the_rebin, the_irreps, **kwargs):    
    
    the_eff_mass_color = '#3535B2'
    the_fit_color = '#b90f22'
    
    ### These are the irreps in this file
    s_irreps = list(the_single_fit_data.keys())
    
    ### If not all the irreps are wanted t be plotted
    ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(the_irreps)
    
    s_irreps = s_irreps[the_first_irrep:the_last_irrep]
    
    ### Loop over the irreps found in the fits file.
    for the_irrep in s_irreps:
        
        print("---------------------------------------------------------------------------")
        print(f"Irrep: {the_irrep}")
        
        ### Central values of the correlators
        the_mean_corr = np.asarray(the_single_correlator_data[f'{the_irrep}/Effective_masses/Mean'])
        
        ### Statistical errors of the correlators
        the_sigmas_corr = np.asarray(the_single_correlator_data[f'{the_irrep}/Effective_masses/Sigmas'])
        
        ### Time extent
        the_nt_corr = np.asarray(the_single_correlator_data[f'{the_irrep}/Time_slices'])
        the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
        
        the_nt_ticks = np.arange(the_nt_corr[0], the_nt_corr[-1], 3)

        ### The SH operator that appears in the plot
        the_op = list(the_single_correlator_data[f'{the_irrep}/Operators'])[0]
        theOperatorNamePlot = vfp.OPERATORS_SH(the_op.decode('utf-8'))
        
        da_irrep = vfp.IrrepInfo(the_irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        nrMomentumIrrep = da_irrep.Momentum
        NameIrrepPlot = da_irrep.NamePlot 
    
        OperatorNamePlot = vfp.SH_OPERATORS_RELABEL(theOperatorNamePlot, NameIrrepPlot, nrMomentumIrrep )
        
        ### This is the info of the fits for this irrep
        dis_set =  np.asarray(the_single_fit_data[f'{the_irrep}/{the_nr_exps}exp/Tmin/{the_type_fit}/Mean'])
        
        ### This is the range of min time slices which the fits were performed
        the_nt_fit = [int(x) for x in dis_set[0]]
        
        ### Max time slice used for the fit
        the_nt_max = int(dis_set[1][0])
        
        ### This is the tmin chosen for this fit. 
        the_chosen_tmin = the_nt_fit.index(the_tmins[the_irreps.index(the_irrep)])
        
        the_plot_fit_data_index = the_nt_fit.index(the_chosen_tmin)
        
        ### Mean values of the fits
        the_fit_data = dis_set[2]
        
        ### Statistical errors of those central values
        the_fit_sigmas = dis_set[3]
        
        ### The Chi^2
        the_chi_corr = dis_set[4]    
        
        ### The title of the Plot
        the_title = f'{OperatorNamePlot} ({the_nr_exps}-exp)'
       
        ### This is just to write the errors properly in the plot        
        the_mean_fit_string = f"{the_fit_data[the_chosen_tmin]:.5f}"
        the_error_string = vfp.WRITTING_ERRORS_PLOTS(the_fit_sigmas[the_chosen_tmin],5)
        the_sigmas_fit_string = the_error_string[0]

        if the_error_string[1]==False:
            the_mean_fit_string = f"{the_fit_data[the_chosen_tmin]:.{the_error_string[2]}f}"
        
        print("Plotting Fitted Effective Masses...")
        
        fit_fig = plt.figure()#plt.figure(figsize=(5.5,4.)) #plt.figure(figsize=(5.5,4.5))             
        the_label = r'$\chi^{2}/\mathrm{d.o.f} = %s$'%f"{the_chi_corr[the_chosen_tmin]:.3f}" + '\n' + r'$E_{\mathrm{fit}} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
        
        vfp.PLOT_FITTED_EFF_MASSES(the_nt, the_mean_corr, the_sigmas_corr, the_fit_data, the_fit_sigmas, the_chosen_tmin, f'{the_rs_scheme} data', the_label, the_title, the_nt_ticks, the_eff_mass_color, the_fit_color)
        plt.show()
        fit_fig.savefig(f'{the_location}Fitted_Effective_Masses_{the_irrep[:4]}_{the_irrep[-1]}_{the_nr_exps}exp_{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        
        
        
def PlotMultiHadronsEffectiveMassesFits(the_multi_hadrons_fit_data, the_matrix_correlator_data, the_quantum_number, the_rs_scheme, the_type_fit, the_nr_exps, the_tmins, the_t0, the_version, the_location, the_rebin, the_irreps, **kwargs):  
    
    m_irreps = list(the_multi_hadrons_fit_data.keys())
    
    the_eff_mass_color = '#3535B2'
    the_fit_color = '#b90f22'
    
    ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(m_irreps)
    
    m_irreps = m_irreps[the_first_irrep:the_last_irrep]     

    ### Loop over the irreps of this file
    for the_irrep in m_irreps:
        
        print("---------------------------------------------------------------------------")
        print(f"Irrep: {the_irrep}")
        
        ### Searching if the fit was done and the gevp plots must be included
        if f'{the_nr_exps}exp' in list(the_multi_hadrons_fit_data[the_irrep].keys()): 
            
            ### Retrieving the data
            the_data_fit = the_multi_hadrons_fit_data[f'{the_irrep}/{the_nr_exps}exp/t0_{the_t0}/Tmin/{the_type_fit}/Mean']
            
            the_data = np.asarray(the_matrix_correlator_data[f'{the_irrep}/GEVP/t0_{the_t0}/Effective_masses/Mean'])
            
            the_data_sigmas = np.asarray(the_matrix_correlator_data[f'{the_irrep}/GEVP/t0_{the_t0}/Effective_masses/Sigmas'])

            ### Loop over the eigenvalues of this irrep
            for bb in range(len(list(the_data_fit.keys()))):
                
                the_mean_corr = the_data[bb]
                
                the_sigmas_corr = the_data_sigmas[bb]
                
                the_nt_corr = np.array(the_matrix_correlator_data[f'{the_irrep}/Time_slices'])
                
                the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
                
                ### bb-th Eigenvalue
                dis_set = np.array(the_data_fit.get(f'lambda_{bb}'))
                
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
                
                the_title = fr'{NameIrrepPlot} ({MomentumIrrep}): $\to \;\lambda_{bb}$ ($t_{{0}} = {the_t0}a$)'

                ### THIS IS THE LABEL WITH TMIN AND TMAX
                # the_label = r'$t_{\mathrm{min}} = %s$'%str(int(the_nt_fit[the_chosen_tmin])) + '\n' + r'$t_{\mathrm{max}}$ = %s'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/\mathrm{d.o.f} = %s$'%np.round(the_chi_corr[the_chosen_tmin],3) + '\n' + r'$E_{\mathrm{fit}} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
                
                the_label = r'$\chi^{2}/\mathrm{d.o.f} = %s$'%np.round(the_chi_corr[the_chosen_tmin],3) + '\n' + r'$E_{\mathrm{fit}} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
                
                fit_fig = plt.figure()                
                vf.PLOT_FITTED_EFF_MASSES(the_nt, the_mean_corr, the_sigmas_corr, the_fit_data, the_fit_sigmas, the_chosen_tmin, the_rs_scheme + ' data', the_label, the_title, the_nt_ticks, the_eff_mass_color, the_fit_color)
                
                fit_fig.savefig(f'{the_location}Fitted_Effective_Masses_{the_quantum_number}_{the_irrep}_{bb}_t0_{the_t0}{the_rebin}_{the_version}.pdf')
                    
        
                    
                    
        
                



### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------



if __name__=="__main__":
    print("Nothing to do here.")
         
