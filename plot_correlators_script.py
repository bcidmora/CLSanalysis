import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys

### LAYOUT AND EXTRA FUNCTIONS
import set_of_layout_functions as vfl
import set_of_analysis_functions as vfa
import set_of_plot_functions as vfp

import warnings
warnings.filterwarnings('ignore')
  

def PlotSingleHadronCorrelators(the_single_correlator_data, the_type_rs, the_version, the_location, the_rebin, **kwargs):
    
    s_irreps = list(the_single_correlator_data.keys())
    
    ### How many irreps do you want to study        
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
    
    if the_type_rs=='jk':
        the_rs_scheme=r'Jackknife'
    elif the_type_rs=='bt':
        the_rs_scheme=r'Bootstrap'   
    
    for irrep in s_irreps:
        ### This is the central values of the correlator
        the_mean_corr = np.asarray(the_single_correlator_data[f'{irrep}/Correlators/Real/Mean'])
        
        ### These are the statistical errors of the correlators
        the_sigmas_corr = np.asarray(the_single_correlator_data[f'{irrep}/Correlators/Real/Sigmas'])
        
        ### This is the time extent
        the_nt = np.asarray(the_single_correlator_data[f'{irrep}/Time_slices'])
        the_nt_ticks = np.arange(the_nt[0]+1, the_nt[-1], int(len(the_nt)/5))

        ### Information about the operators
        the_op_list = list(the_single_correlator_data[f'{irrep}/Operators'])[0]
        theOperatorNamePlot = vfp.OPERATORS_SH(the_op_list.decode('utf-8'))        

        ### Information about the irrep
        da_irrep = vfp.IrrepInfo(irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        nrMomentumIrrep = da_irrep.Momentum
        NameIrrepPlot = da_irrep.NamePlot
        
        OperatorNamePlot = vfp.SH_OPERATORS_RELABEL(theOperatorNamePlot, NameIrrepPlot, nrMomentumIrrep )
        
        print('Correlator plot in progress...')
        the_corr_fig = plt.figure()
        vfp.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$\mathbb{Re}\;C(t)$', 'o', OperatorNamePlot)
        plt.show()
        the_corr_fig.savefig(f'{the_location}Correlator_{irrep[:4]}_{irrep[-1]}{the_rebin}_{the_version}.pdf', bbox_inches='tight')
    
        print('Correlator Log-plot in process...')
        the_log_corr_fig = plt.figure()
        vfp.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$\log \mathbb{Re}\;C(t)$', 'o', OperatorNamePlot, yscale='log')
        plt.show()
        the_log_corr_fig.savefig(f'{the_location}Correlator_{irrep[:4]}_{irrep[-1]}_log{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        
        print('Correlator histogram in process...')
        tt = int(len(the_nt)/3)+1
        the_gauss_fig = plt.figure()
        the_nt_mean = the_mean_corr[tt]
        the_rs = np.asarray(the_single_correlator_data[f'{irrep}/Correlators/Real/Resampled'])[tt]
        the_nr_bins = int(len(the_rs)*.05)

        the_mean_rs = np.mean(the_rs)
        the_means_dif = np.abs(the_nt_mean - the_mean_rs)
        the_stat_error = the_sigmas_corr[tt]

        vfp.PLOT_HISTOGRAMS(the_rs, r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error), the_mean_rs, r'$ \bar{C}_{%s}(t) = $'%the_type_rs + r' $%s$'%'{:.10e}'.format(the_mean_rs), the_nt_mean, r'$ \bar{C}(t) = %s$'%'{:.10e}'.format(the_nt_mean), OperatorNamePlot + rf' $\to$ $t = {tt+the_nt[0]} a$', the_nr_bins,  r'Correlator')
        plt.show()
        the_gauss_fig.savefig(f'{the_location}Histogram_correlators_{irrep[:4]}_{irrep[-1]}{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        

def PlotMultiHadronCorrelators(the_matrix_correlator_data, the_quantum_number, the_type_rs, the_version, the_t0, the_location, the_rebin, **kwargs):
    
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
        
    ### Choosing the resampling scheme to put in the legend
    if the_type_rs=='jk':
        the_rs_scheme='Jackknife'
    elif the_type_rs=='bt':
        the_rs_scheme='Bootstrap'
    
    ### Loop over the irreducible representations
    for irrep in m_irreps:
        ### The list of operators of this irrep
        the_op_list = list(the_matrix_correlator_data[f'{irrep}/Operators'])
        
        ### The correaltor dataset
        the_data_corr = vfa.RESHAPING_CORRELATORS(np.asarray(the_matrix_correlator_data[f'{irrep}/Correlators/Real/Mean']))
        
        ### The sigmas of this correlator dataset
        the_data_sigmas_corr = np.asarray(the_matrix_correlator_data[f'{irrep}/Correlators/Real/Sigmas'])
        
        ### This is the time interval
        the_nt = np.asarray(the_matrix_correlator_data[f'{irrep}/Time_slices'])

        ### These are going to be the ticks in the x-label of the plots
        the_nt_ticks = np.arange(the_nt[0]+1, the_nt[-1], int(len(the_nt)/5))
        
        ### Plotting eigenvalues in the full time range or only starting from t0
        if kwargs.get('full_range_nt')==None: the_start_nt = the_t0-the_nt[0]
        else: the_start_nt = the_nt[0]
        
        ### Information of the irrep to write it properly in the plots.
        da_irrep = vfp.IrrepInfo(irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot 
        NameIrrep = da_irrep.Name
        
        if the_diagonal_corrs_flag:
            
            ### Loop over the number of operators of this matrix
            for bb in range(len(the_op_list)):      
                
                ### The specific operator
                the_op = the_op_list[bb]
                theOperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
                OperatorNamePlot = vfp.MH_OPERATORS_RELABEL(theOperatorNamePlot)
                print('Correlator plots in process...')
                    
                ### Plotting the diagonal correlators
                corr_fig = plt.figure()
                vfp.PLOT_CORRELATORS(the_nt, the_data_corr[bb,bb], the_data_sigmas_corr[bb], the_rs_scheme, the_nt_ticks, r'$t$', r'$\mathbb{Re}\;C(t)$', 'o', f'{NameIrrepPlot} ({MomentumIrrep}) ' + r' $\to \;C_{%s}$'%(str(bb)+str(bb)) + f'= {OperatorNamePlot}')
                plt.show()
                corr_fig.savefig(f'{the_location}DiagonalCorrelator_{the_quantum_number}_{irrep}_{bb}{the_rebin}_{the_version}.pdf')
                
                ## Plotting the log of the diagonal correlators.
                print('Correlator Log-plots in progress...')
                corr_fig = plt.figure()
                vfp.PLOT_CORRELATORS(the_nt, the_data_corr[bb,bb], the_data_sigmas_corr[bb], the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$\log\mathbb{Re}\;C(t)$', 'o', f'{NameIrrepPlot} ({MomentumIrrep}) ' + r' $\to \;C_{%s}$'%(str(bb)+str(bb)) + f'= {OperatorNamePlot}', yscale='log')
                plt.show()
                corr_fig.savefig(f'{the_location}DiagonalCorrelator_{the_quantum_number}_{irrep}_{bb}_log{the_rebin}_{the_version}.pdf')
                
                ### Plotting the histogram at a certain time slice t
                print('Correlator histogram in progress...')
                tt = int(len(the_nt)/3)+1
                the_gauss_fig = plt.figure()
                the_nt_mean = the_data_corr[bb,bb,tt]
                the_rs = vfa.RESHAPING_CORRELATORS_RS_NT(np.asarray(the_matrix_correlator_data[f'{irrep}/Correlators/Real/Resampled']))[bb,bb,tt]
                the_nr_bins = int((np.asarray(the_matrix_correlator_data[f'{irrep}/Correlators/Real/Resampled']).shape[1])*.05)
                
                ### Here the mean value of the sampling data and the mean value are compared to check the quality of the resampled data
                the_mean_rs = np.mean(the_rs)
                the_means_dif = np.abs(the_nt_mean - the_mean_rs)
                the_stat_error = the_data_sigmas_corr[bb,tt]
                
                ### Plotting the histogram now
                vfp.PLOT_HISTOGRAMS(the_rs, r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error), the_mean_rs, r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r' $%s$'%'{:.10e}'.format(the_mean_rs), the_nt_mean, r'$ \bar{C}(t) = %s$'%'{:.10e}'.format(the_nt_mean), f'{NameIrrepPlot} ({MomentumIrrep}) ' +  r'$\to \;C_{%s}$'%(str(bb) + str(bb)) + r' $(t = %sa) = $'%(tt+the_nt[0]) + f' {OperatorNamePlot}', the_nr_bins, r'$C(t)$')
                plt.show()
                the_gauss_fig.savefig(f'{the_location}Histogram_DiagCorrelator_{the_quantum_number}_{irrep}_{bb}{the_rebin}_{the_version}.pdf')
            
            ### The Diagonal of the correlators are plotted all together with their errors to compare them directly. 
            if all_corr_flag:
                corr_fig = plt.figure()
                print('ALL Correlators Log-plot in progress...')
                ### Loop over each of the entries of the diagonal of the correlation matrix
                for bb in range(len(the_op_list)):
                    
                    ### Name of this operator
                    the_op = the_op_list[bb]
                    theOperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                    OperatorNamePlot = vf.MH_OPERATORS_RELABEL(theOperatorNamePlot)           
                    
                    plt.errorbar(the_nt, the_data_corr[bb][bb], the_data_sigmas_corr[bb],  marker=vfp.the_markers_list[bb], ls='None', ms=6, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, capsize=5, label = OperatorNamePlot, color=vfp.the_colors[bb])
                plt.xlabel(r'$t\,/\,a$',fontsize=36)
                plt.ylabel(r'$\log\mathbb{Re}\;C(t)$', fontsize=36)
                plt.title( f'{NameIrrepPlot} ({MomentumIrrep}) ' ,fontsize=26)
                # the_n_cols = int(len(the_op_list)/2)
                plt.legend(fontsize=16, ncol= 2, columnspacing=0.1, handletextpad=0.01)
                plt.yscale('log')
                plt.tight_layout()
                # if len(the_op_list)>7: the_n_cols = int(len(the_op_list)/3)
                # else: the_n_cols = int(len(the_op_list)/2)
                plt.xticks(the_nt_ticks,fontsize=18)
                plt.yticks(fontsize=18)
                plt.show()
                corr_fig.savefig(f'{the_location}ALLDiagonalCorrelators_{the_quantum_number}_{irrep}_log{the_rebin}_{the_version}.pdf')
            
        ### Here the Eigenvalues are plotted all together too in a log-plot
        if 'GEVP' in list(the_matrix_correlator_data[irrep].keys()):
            
            the_data = np.asarray(the_matrix_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Eigenvalues/Mean'])
            the_data_sigmas = np.asarray(the_matrix_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Eigenvalues/Covariance_matrix'])
            
            for bb in range(len(the_data)):
                ### The mean value of the eigenvalue_{i}
                the_mean_corr = the_data[bb]
                
                ### The corresponding sigmas of this eigenvalue
                the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                
                print("..................................................")
                print(f'Eigenvalue = {bb} plot in process...')
                
                corr_fig = plt.figure()
                ### Plotting the eigenvalues one by one
                vfp.PLOT_CORRELATORS(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], the_sigmas_corr[the_start_nt:], the_rs_scheme, the_nt_ticks, r'$t$', r'$\lambda_{i}(t)$', 'o',  f'{NameIrrepPlot} ({MomentumIrrep}): ' + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %s$)'%str(the_t0))
                plt.show()
                corr_fig.savefig(f'{the_location}Eigenvalues_{the_quantum_number}_{irrep}_{bb}_t0_{the_t0}{the_rebin}_{the_version}.pdf')
                
                ### Plotting the eigenvalues log-plots one by one
                print('Eigenvalue = %s Log-plot in progress...'%str(bb))
                corr_fig = plt.figure()
                vfp.PLOT_CORRELATORS(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], the_sigmas_corr[the_start_nt:], the_rs_scheme, the_nt_ticks, r'$t\,/\,a$', r'$\log\,(\lambda_{i}(t))$', 'o',  f'{NameIrrepPlot} ({MomentumIrrep}) ' + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %sa$)'%str(the_t0), yscale='log')
                plt.show()
                corr_fig.savefig(f'{the_location}Eigenvalues_{the_quantum_number}_{irrep}_{bb}_log_t0_{the_t0}{the_rebin}_{the_version}.pdf')
                
                ### Plotting the histograms of these eigenvalues
                print(f'Eigenvalue = {bb} histogram in progress...')
                tt = int(len(the_nt)/3)+1
                the_gauss_fig = plt.figure()
                the_nt_mean = the_mean_corr[tt]
                the_rs = np.asarray(the_matrix_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Eigenvalues/Resampled'])[bb].transpose()[tt]
                the_nr_bins = int(np.asarray(the_matrix_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Eigenvalues/Resampled']).shape[1] * .05)

                the_mean_rs = np.mean(the_rs)
                the_means_dif = np.abs(the_nt_mean - the_mean_rs)
                the_stat_error = the_sigmas_corr[tt]
                
                vfp.PLOT_HISTOGRAMS(the_rs, r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error), the_mean_rs, r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r' $%s$'%'{:.10e}'.format(the_mean_rs), the_nt_mean, r'$ \bar{C}(t) = %s$ '%'{:.10e}'.format(the_nt_mean), f'{NameIrrepPlot} ({MomentumIrrep}) ' +  r'$\to\; \lambda_{%s}$'%str(bb) + r'$(t_{0} = %sa)$'%str(the_t0) +  r' $[t = %sa]$'%(tt+the_nt[0]), the_nr_bins, r'Eigenvalue ($\lambda_{%s}$)'%str(bb))
                plt.show()
                the_gauss_fig.savefig(f'{the_location}Histogram_Eigenvalues_{the_quantum_number}_{irrep}_{bb}_t0_{the_t0}{the_rebin}_{the_version}.pdf')
            
            if all_corr_flag:
                print("..................................................\n")
                
                print('ALL Eigenvalues Log-plot in progress...')
                corr_fig = plt.figure()           
                for bb in range(len(the_data)):
                    the_mean_corr = the_data[bb]
                    the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                    
                    plt.errorbar(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], yerr = the_sigmas_corr[the_start_nt:], marker=vfp.the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=3.5, label = r'$\lambda_{%s}$'%str(bb), color=vfp.the_colors[bb])
                plt.xlabel(r'$t\,/\,a$',fontsize=28)
                plt.ylabel(r'$\log\,(\lambda_{i}(t))$',fontsize=28)
                plt.title( f'{NameIrrepPlot} ({MomentumIrrep}) ' + r'$\to\;\lambda_{i} (t_{0} = %sa$)'%str(the_t0), fontsize=18)
                plt.yscale('log')
                plt.legend(fontsize=16, ncol= 2, handletextpad=0.01)
                plt.tight_layout()
                # if len(the_data)>6: the_n_cols = int(len(the_data)/3)
                # else: the_n_cols = int(len(the_data)/2)
                # plt.legend(fontsize=12, ncol=the_n_cols, handletextpad=0.3)
                plt.xticks(the_nt_ticks)
                plt.show()
                corr_fig.savefig(f'{the_location}ALLEigenvalues_{the_quantum_number}_{irrep}_log_t0_{the_t0}{the_rebin}_{the_version}.pdf')



def PlotRatioHadronCorrelators(the_ratio_correlator_data, the_quantum_number, the_type_rs, the_version, the_t0, the_location, the_rebin, **kwargs):
    
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
    
    # the_nr_bins = 25
    for irrep in mr_irreps:
        
        ### The operators list
        the_op_list = list(the_ratio_correlator_data[f'{irrep}/Operators'])
        
        ### the non-int levels
        the_non_int = list(the_ratio_correlator_data[f'{irrep}/Single_hadron_corrs'])
        
        ### the time slices
        the_nt = np.asarray(the_ratio_correlator_data[f'{irrep}/Time_slices'])
        
        ### The data with shape (N, nr. non-int, nt)
        the_data = np.asarray(the_ratio_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Eigenvalues/Mean'])
        
        ### The sigmas
        the_data_sigmas = np.asarray(the_ratio_correlator_data[f'{irrep}/GEVP/t0_{the_t0}/Eigenvalues/Covariance_matrix'])
        
        ### The shape of the data
        the_data_shape = the_data.shape
        
        ### Loop over the eigenvalues
        for bb in range(the_data_shape[0]):
            
            print('Correlator plot in process...')
            the_corr_fig = plt.figure()
            
            ### Loop over the non-interacting levels
            for nn in range(the_data_shape[1]):
                
                ### The data
                the_mean_corr = the_data[bb, nn]
                
                ### the non-interacting levels
                the_non_int_n = vfp.NON_INTERACTING_LABELS(the_non_int[nn].decode('utf-8'))
                
                ### The errors
                the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb,nn]))
                
                the_nt_ticks = np.arange(the_nt[0], the_nt[-1], 5)

                the_op = the_op_list[bb]
                OperatorNamePlot = vfp.OPERATORS_MH(the_op.decode('utf-8'))
                da_irrep = vfp.IrrepInfo(irrep)
                
                MomentumIrrep = da_irrep.TotalMomPlot
                NameIrrepPlot = da_irrep.NamePlot
                NameIrrep = da_irrep.Name
                
                plt.errorbar(the_nt, the_mean_corr, yerr = the_sigmas_corr, marker=vfp.the_markers_list[nn], ls='None', ms=4.5, markeredgewidth=1.5, lw=1.5, elinewidth=1.5, zorder=3, capsize=3.5, label = the_non_int_n, color=vfp.the_colors[nn])
            plt.xlabel(r'$t/a$', fontsize=24)
            plt.ylabel(r'$\mathbb{Re}\;C(t)$', fontsize=24)
            plt.title( fr'{NameIrrepPlot} ({MomentumIrrep}) $\to \;\lambda_{bb}$ ($t_{{0}} = {the_t0}$)', fontsize=20)
            plt.xticks(the_nt_ticks, fontsize=14)
            plt.yticks(fontsize=14)
            if the_data_shape[1]>4: the_n_cols = int(the_data_shape[1]/3)
            else: the_n_cols = int(the_data_shape[1]/2)
            plt.legend(fontsize=14, ncol=the_n_cols, handletextpad=0.01)
            plt.tight_layout()
            plt.show()
            the_corr_fig.savefig(f'{the_location}/Eigenvalues_ratios_{the_quantum_number}_{irrep}_{bb}_{the_rebin}_{the_version}.pdf')
                



### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------




if __name__=="__main__":
    print("Nothing to do here.")
