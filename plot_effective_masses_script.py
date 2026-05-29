import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_layout_functions as vfl
import set_of_analysis_functions as vfa
import set_of_plot_functions as vfp


def PlotSingleHadronsEffectiveMasses(the_single_correlator_data, the_rs_scheme, the_version, the_location, the_rebin, **kwargs):
    
    s_irreps = list(the_single_correlator_data.keys())
    
    # ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(s_irreps)
    
    s_irreps = s_irreps[the_first_irrep:the_last_irrep]
    
    ### Loop over the SH irreps
    for the_irrep in s_irreps:
        ### Central values of the correlators
        the_mean_corr = np.array(the_single_correlator_data[f'{the_irrep}/Effective_masses/Mean'])
        
        ### Statistical errors of the correlators
        the_sigmas_corr = np.array(the_single_correlator_data[f'{the_irrep}/Effective_masses/Sigmas'])
        
        ##3 Time extent
        the_nt_corr = np.array(the_single_correlator_data[f'{the_irrep}/Time_slices'])
        the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
        
        the_nt_ticks = np.arange(the_nt_corr[2], the_nt_corr[-1], 2)
        
        ### The SH operator that appears in the plot
        the_op = list(the_single_correlator_data[f'{the_irrep}/Operators'])[0]
        theOperatorNamePlot = vfp.OPERATORS_SH(the_op.decode('utf-8'))
        
        da_irrep = vfp.IrrepInfo(the_irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        nrMomentumIrrep = da_irrep.Momentum
        NameIrrepPlot = da_irrep.NamePlot 
    
        OperatorNamePlot = vfp.SH_OPERATORS_RELABEL(theOperatorNamePlot, NameIrrepPlot, nrMomentumIrrep )
        
        print('Effective Mass plot in progress...')
        # the_efm_fig = plt.figure(figsize=(4,3))
        the_efm_fig = plt.figure()
        vfp.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o',  OperatorNamePlot)
        plt.show()
        the_efm_fig.savefig(f'{the_location}EffectiveMass_{the_irrep}{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        
        
def PlotMultiHadronsEffectiveMasses(the_matrix_correlator_data, the_quantum_number, the_rs_scheme, the_version, the_t0, the_location, the_rebin, **kwargs):
    ### Getting all the irreps in this ensemble
    m_irreps = list(the_matrix_correlator_data.keys())
    
    ### These variables are to plot the GEVP or the operators analysis eigenvalues
    the_diagonal_corrs_flag = kwargs.get('diag_corrs')
    all_corr_flag = kwargs.get('all_corr')
    
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
    
    ### Loop over the irreps
    for the_irrep in m_irreps:
        
        ### Operator list of that specific irrep
        the_op_list = list(the_matrix_correlator_data[f'{the_irrep}/Operators'])
        
        ### Time slices 
        the_nt = np.asarray(the_matrix_correlator_data[f'{the_irrep}/Time_slices'])
            
        ### The new time slices range for the effective masses
        the_nt_corr_efm = np.arange(the_nt[0]+0.5, the_nt[-1]+0.5, 1)
        
        ### These are the ticks in the x-axis
        the_nt_ticks = np.arange(the_nt[0]+1, the_nt[-1], int(len(the_nt)/5))
        
        ### Information about the irrep written for the plots
        da_irrep = vfp.IrrepInfo(the_irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot
        
        if the_diagonal_corrs_flag:

            ### Effective masses of correlators
            the_data_corr = np.asarray(the_matrix_correlator_data[f'{the_irrep}/Correlators/Real/Effective_masses'])
            
            ### Statistical error of the effective masses of the central values of the correlators
            the_data_sigmas_corr = np.asarray(the_matrix_correlator_data[f'{the_irrep}/Correlators/Real/Effective_masses_sigmas'])

            ### Loop over the size of the correlation matrix
            for bb in range(len(the_op_list)):
                ### The effective masses of the mean values of the diagonal of the correlators
                the_mean_efm = the_data_corr[bb]
                
                ### Their sigmas
                the_sigmas_efm = the_data_sigmas_corr[bb]
                
                ### The operator of this dataset
                the_op = the_op_list[bb]
                
                ### Convenient name for the plots
                theOperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
                OperatorNamePlot = vfp.MH_OPERATORS_RELABEL(theOperatorNamePlot)
                
                print('Effective mass diagonal correlators plots in process...')
                
                ### Checking for the data
                the_ymin = vfp.CHOOSING_YMIN_PLOT(the_mean_efm)
                
                efm_corr_fig = plt.figure(figsize=(6,4))
                vfp.PLOT_CORRELATORS(the_nt_corr_efm[1:], the_mean_efm[1:], the_sigmas_efm[1:], the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o', f'{NameIrrepPlot} ({MomentumIrrep}) ' + rf'$\to$ {OperatorNamePlot}', ymin=the_ymin)
                plt.show()
                efm_corr_fig.savefig(f'{the_location}EffectiveMass_DiagonalCorrelators_{the_quantum_number}_{the_irrep}_{bb}{the_rebin}_{the_version}.pdf', bbox_inches='tight')
            
            if all_corr_flag:
                ### Here all the diagonal of the correlators are put together
                efm_corr_all_fig = plt.figure()
                print('Effective mass ALL diagonal correlators together plot in process...')
                
                ### Loop over the operators of the correlation matrix
                for bb in range(len(the_op_list)):  
                    ### The diagonal of the correlator
                    the_mean_efm = the_data_corr[bb]
                    the_sigmas_efm = the_data_sigmas_corr[bb]

                    the_op = the_op_list[bb]
                    theOperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
                    OperatorNamePlot = vfp.MH_OPERATORS_RELABEL(theOperatorNamePlot)
                    
                    plt.errorbar(the_nt_corr_efm, the_mean_efm, yerr = the_sigmas_efm, marker=vfp.the_markers_list[bb], ls='None', ms=8, markeredgewidth=2, lw=1.5, elinewidth=1.5, capsize=6, label = OperatorNamePlot, color=vfp.the_colors[bb])
                plt.xlabel(r'$t\,/\,a$', fontsize=36)
                plt.ylabel(r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=36)

                plt.title( f'{NameIrrepPlot} ({MomentumIrrep}) ', fontsize=28)
                plt.xticks(the_nt_ticks,fontsize=18)
                plt.yticks(fontsize=18)
                plt.legend(fontsize=16, ncol=2, handletextpad=0.1, columnspacing=0.1)
                plt.ylim([0.3,1.5])
                plt.grid(True, alpha=0.2)
                plt.tight_layout()
                plt.show()
                efm_corr_all_fig.savefig(f'{the_location}EffectiveMass_ALLDiagonalCorrelators_{the_quantum_number}_{the_irrep}{the_rebin}_{the_version}.pdf', bbox_inches='tight')

        ### If GEVP was performed, the eigenvalues are also going to be plotted.
        if 'GEVP' in list(the_matrix_correlator_data[the_irrep].keys()):
            
            the_data = np.asarray(the_matrix_correlator_data[f'{the_irrep}/GEVP/t0_{the_t0}/Effective_masses/Mean'])
            the_data_sigmas = np.asarray(the_matrix_correlator_data[f'{the_irrep}/GEVP/t0_{the_t0}/Effective_masses/Sigmas'])

            ### Loop over the eigenvalues
            for bb in range(len(the_data)):
                the_mean_corr = the_data[bb]
                the_sigmas_corr = the_data_sigmas[bb]
                
                print('Effective mass eigenvalues plots in process...')
                
                 ### Checking for the data
                the_ymin = vfp.CHOOSING_YMIN_PLOT(the_mean_corr)
                
                efm_fig = plt.figure()
                vfp.PLOT_CORRELATORS(the_nt_corr_efm, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o',  f'{NameIrrepPlot} ({MomentumIrrep}) ' + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %sa$)'%str(the_t0), ymin=the_ymin)
                plt.show()
                efm_fig.savefig(f'{the_location}EffectiveMass_Eigenvalues_{the_quantum_number}_{the_irrep}_{bb}_t0_{the_t0}{the_rebin}_{the_version}.pdf', bbox_inches='tight')
            
            if all_corr_flag:
                ### All Eigenvalues in only one plot
                efm_corr_all_fig = plt.figure()
                
                print('Effective mass ALL eigenvalues together plot in process...')

                the_min_position = np.where(the_data[0] == min(the_data[0][:-3]))
                the_max_position = np.where(the_data[0] == max(the_data[-1][:-3]))
        
                the_min_y = (the_data[0][the_min_position]-the_data_sigmas[0][the_min_position])*.95
                the_max_y= (the_data[-1][the_max_position]+the_data_sigmas[-1][the_max_position])*1.05
                
                ### Loop over the eigenvalues
                for bb in range(len(the_data)):   
                    the_mean_efm = the_data[bb]
                    the_sigmas_efm = the_data_sigmas[bb]

                    the_op = the_op_list[bb]
                    theOperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
                    OperatorNamePlot = vfp.MH_OPERATORS_RELABEL(theOperatorNamePlot)
                    
                    plt.errorbar(the_nt_corr_efm[1:], the_mean_efm[1:], yerr = the_sigmas_efm[1:], marker=vfp.the_markers_list[bb], ls='None', ms=7, markeredgewidth=2, lw=1.5, elinewidth=1.5, capsize=6, label = rf'$\lambda_{bb}$', color=vfp.the_colors[bb])
                plt.xlabel(r'$t\,/\,a$', fontsize=20)
                plt.ylabel(r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=20)
                plt.title( fr'{NameIrrepPlot} ({MomentumIrrep}) $\to\, \lambda_{{i}}(t,t_{{0}} = {the_t0}a)$', fontsize=20)
                plt.xticks(the_nt_ticks, fontsize=12)
                plt.yticks(fontsize=12)
                if len(the_data)>6: the_n_cols = int(len(the_data)/3)
                else: the_n_cols = int(len(the_data)/2)

                plt.legend(fontsize=14, ncol=the_n_cols, handletextpad=0.05, columnspacing=0.3)
                plt.grid(True, alpha=0.2)
                plt.tight_layout()
                plt.show()
                efm_corr_all_fig.savefig(f'{the_location}EffectiveMass_ALLEigenvalues_{the_quantum_number}_{the_irrep}_t0_{the_t0}{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        
            
            
def PlotRatioHadronsEffectiveMasses(the_ratio_correlator_data, the_quantum_number, the_rs_scheme, the_version, the_t0, the_location, the_rebin, **kwargs):
    mr_irreps = list(the_ratio_correlator_data.keys())
    
    ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(mr_irreps)
    
    mr_irreps = mr_irreps[the_first_irrep:the_last_irrep]    
    
    for irrep in mr_irreps:
        
        ### The list of operatos
        the_op_list = list(the_ratio_correlator_data[f'{irrep}/Operators'])
        
        ### The correlators time slices
        the_nt = np.asarray(the_ratio_correlator_data[f'{irrep}/Time_slices'])
        
        ### The non-interacting levels
        the_non_int_hads = np.asarray(the_ratio_correlator_data[f'{irrep}/Single_hadron_corrs'])
        
        ### The effective masses array
        the_data = np.asarray(the_ratio_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Effective_masses/Mean']) # Shape (N, nr. SH, nt)
        
        ### The effective masses errors
        the_data_sigmas = np.asarray(the_ratio_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Effective_masses/Sigmas']) # Shape (N, nr. SH, nt)
        
        ### Just some helper variable
        the_data_shape = the_data.shape
        
        ### The new time slices range for the effective masses
        the_nt_corr_efm = np.arange(the_nt[0]+0.5, the_nt[-1]+0.5, 1)
        the_nt_ticks = np.arange(the_nt[0]+1, the_nt[-1], int(len(the_nt)/5))

        ### Loop over each eigenvalue
        for bb in range(the_data_shape[0]):
            
            the_op = the_op_list[bb]
            OperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
            da_irrep = vfp.IrrepInfo(irrep)
            MomentumIrrep = da_irrep.TotalMomPlot
            NameIrrepPlot = da_irrep.NamePlot
            NameIrrep = da_irrep.Name
            
            print('Effective Mass correlators ratio plot in progress...')
            efm_fig = plt.figure()
            
            ### Loop over the number of non-interacting levels
            for nn in range(the_data_shape[1]):
                the_mean_corr = the_data[bb, nn]
                the_sigmas_corr = the_data_sigmas[bb,nn]
                
                the_non_int = vfp.NON_INTERACTING_LABELS(the_non_int_hads[nn].decode('utf-8'))
                
                plt.errorbar(the_nt_corr_efm, the_mean_corr, yerr = the_sigmas_corr, marker=vfp.the_markers_list[nn], ls='None', ms=5.5, markeredgewidth=1.5, lw=1.5, elinewidth=1.5, zorder=3, capsize=3.5, label = the_non_int, color=vfp.the_colors[nn])
                
            plt.xlabel(r'$t\,/\,a$', fontsize=24)
            plt.ylabel(r'$a \;m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=24)
            plt.title(fr'{NameIrrepPlot} ({MomentumIrrep}) $\to \;\lambda_{bb}$ ($t_{{0}} = {the_t0}$)', fontsize=20)
            plt.xticks(the_nt_ticks, fontsize=14)
            plt.yticks(fontsize=14)
            if the_data_shape[1]>4: the_n_cols = int(the_data_shape[1]/3)
            else: the_n_cols = int(the_data_shape[1]/2)
            plt.legend(fontsize=14, ncol=the_n_cols, handletextpad=0.01)
            plt.tight_layout()
            plt.grid(True, alpha=0.2)
            plt.show()
            efm_fig.savefig(f'{the_location}EffectiveMass_Eigenvalues_ratios_{the_quantum_number}_{irrep}_{bb}_t0_{the_t0}{the_rebin}_{the_version}.pdf', bbox_inches='tight')




### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------


### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------


if __name__=="__main__":
    print("Nothing to do here.")
