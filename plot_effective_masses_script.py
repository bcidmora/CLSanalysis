import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_functions as vf


def PlotSingleHadronsEffectiveMasses(the_single_correlator_data, the_rs_scheme, the_version, the_location, the_rebin, **kwargs):
    
    s_irreps = list(the_single_correlator_data.keys())
    
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
    
    ### Loop over the SH irreps
    for the_irrep in s_irreps:
        ### Central values of the correlators
        the_mean_corr = np.array(the_single_correlator_data[the_irrep + '/Effective_masses/Mean'])
        
        ### Statistical errors of the correlators
        the_sigmas_corr = np.array(the_single_correlator_data[the_irrep + '/Effective_masses/Sigmas'])
        
        ##3 Time extent
        the_nt_corr = np.array(the_single_correlator_data[the_irrep + '/Time_slices'])
        the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
        the_nt_ticks = np.arange(int(the_nt_corr[-1]/5), the_nt_corr[-1], int(the_nt_corr[-1]/5))
        # the_nt_ticks = np.arange(the_nt_corr[2], the_nt_corr[-1], 2)
        
        ### The SH operator that appears in the plot
        the_op = list(the_single_correlator_data[the_irrep+'/Operators'])[0]
        OperatorNamePlot = vf.OPERATORS_SH(the_op.decode('utf-8'))
        
        da_irrep = vf.IrrepInfo(the_irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot
        
        print('Effective Mass plot in progress...')
        the_efm_fig = plt.figure()
        vf.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t$', r'$a_{t} \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o',  OperatorNamePlot + ' (%s): '%MomentumIrrep + r' $\to$ %s'%NameIrrepPlot)
        the_efm_fig.savefig(the_location + 'EffectiveMass_' + the_irrep[:4] +'_%s'%the_irrep[-1] + the_rebin + '_v%s.pdf'%the_version)
        
        
def PlotMultiHadronsEffectiveMasses(the_matrix_correlator_data, the_quantum_number, the_rs_scheme, the_version, the_t0, the_location, the_rebin, **kwargs):
    ### Getting all the irreps in this ensemble
    m_irreps = list(the_matrix_correlator_data.keys())
    
    ### These variables are to plot the GEVP or the operators analysis eigenvalues
    the_diagonal_corrs_flag = kwargs.get('diag_corrs')
    the_gevp_flag = kwargs.get('gevp')
    the_operators_analysis_flag = kwargs.get('ops_analysis')    
    
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
    
    ### Putting a lot of markers in case of a big matrix
    the_markers_list = ['o','v','s','p','^','*','x','d','>','D', '<','8','P','h','1','o','v','s','p','^','*','x']
    
    ### Loop over the irreps
    for the_irrep in m_irreps:
        
        ### Operator list of that specific irrep
        the_op_list = list(the_matrix_correlator_data[the_irrep+'/Operators'])
        
        ### Time slices 
        the_nt = np.array(the_matrix_correlator_data[the_irrep + '/Time_slices'])
            
        ### Effective masses of correlators
        the_data_corr = np.array(the_matrix_correlator_data[the_irrep + '/Correlators/Real/Effective_masses'])
        
        ### Statistical error of the effective masses of the central values of the correlators
        the_data_sigmas_corr = np.array(the_matrix_correlator_data[the_irrep + '/Correlators/Real/Effective_masses_sigmas'])

        ### The new time slices range for the effective masses
        the_nt_corr_efm = np.arange(the_nt[0]+0.5, the_nt[-1]+0.5, 1)
        
        ### These are the ticks in the x-axis
        the_nt_ticks = np.arange(the_nt[0]+1, the_nt[-1], int(len(the_nt)/5))
        
        ### Information about the irrep written for the plots
        da_irrep = vf.IrrepInfo(the_irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot
        
        if the_diagonal_corrs_flag:
            ### Loop over the size of the correlation matrix
            for bb in range(len(the_op_list)):
                ### The effective masses of the mean values of the diagonal of the correlators
                the_mean_efm = the_data_corr[bb]
                
                ### Their sigmas
                the_sigmas_efm = the_data_sigmas_corr[bb]
                
                ### The operator of this dataset
                the_op = the_op_list[bb]
                
                ### Convenient name for the plots
                OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                
                print('Effective mass diagonal correlators plots in process...')
                
                ### Checking for the data
                the_ymin = vf.CHOOSING_YMIN_PLOT(the_mean_efm)
                
                ### Plotting
                efm_corr_fig = plt.figure()
                vf.PLOT_CORRELATORS(the_nt_corr_efm, the_mean_efm, the_sigmas_efm, the_rs_scheme, the_nt_ticks, r'$t$', r'$a_{t} \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o', NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to \;C_{%s}$'%(str(bb)+str(bb)) + ' = '+OperatorNamePlot, ymin=the_ymin)
                # plt.show()
                efm_corr_fig.savefig(the_location + 'EffectiveMass_DiagonalCorrelators' + the_quantum_number + the_irrep + '_%s'%str(bb) + the_rebin + '_v%s.pdf'%the_version)
            
            ### Here all the diagonal of the correlators are put together
            efm_corr_all_fig = plt.figure()
            print('Effective mass ALL diagonal correlators together plot in process...')
            
            ### Loop over the operators of the correlation matrix
            for bb in range(len(the_op_list)):  
                ### The diagonal of the correlator
                the_mean_efm = the_data_corr[bb]
                the_sigmas_efm = the_data_sigmas_corr[bb]

                the_op = the_op_list[bb]
                OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                
                plt.errorbar(the_nt_corr_efm, the_mean_efm, yerr = the_sigmas_efm, marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=2, lw=1.5, elinewidth=1.5, zorder=3, capsize=2.75, label = OperatorNamePlot)
            plt.xlabel(r'$t$', fontsize=14)
            plt.ylabel(r'$a_{t} \;m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=14)
            plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to\;C_{ii}(t)$'+ ' [' + the_rs_scheme + ']')
            plt.xticks(the_nt_ticks)
            if len(the_op_list)>10: the_n_cols = int(len(the_op_list)/3)
            else: the_n_cols = int(len(the_op_list)/2)
            plt.ylim([0.3,1.5])
            plt.legend(fontsize=12, ncol=the_n_cols, handletextpad=0.3)
            plt.tight_layout()
            # plt.show()
            efm_corr_all_fig.savefig(the_location + 'EffectiveMass_ALLDiagonalCorrelators' + the_quantum_number + the_irrep + the_rebin + '_v%s.pdf'%the_version)

        ### If GEVP was performed, the eigenvalues are also going to be plotted.
        if 'GEVP' in list(the_matrix_correlator_data[the_irrep].keys()) and the_gevp_flag:
            
            the_data = np.array(the_matrix_correlator_data[the_irrep + '/GEVP/t0_%s/Effective_masses/Mean'%the_t0])
            the_data_sigmas = np.array(the_matrix_correlator_data[the_irrep + '/GEVP/t0_%s/Effective_masses/Sigmas'%the_t0])

            ### Loop over the eigenvalues
            for bb in range(len(the_data)):
                the_mean_corr = the_data[bb]
                the_sigmas_corr = the_data_sigmas[bb]
                
                print('Effective mass eigenvalues plots in process...')
                
                 ### Checking for the data
                the_ymin = vf.CHOOSING_YMIN_PLOT(the_mean_corr)
                
                efm_fig = plt.figure()
                vf.PLOT_CORRELATORS(the_nt_corr_efm, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t$', r'$a_{t} \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o',  NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %s$)'%str(the_t0), ymin=the_ymin)
                # plt.show()
                efm_fig.savefig(the_location + 'EffectiveMass_Eigenvalues' + the_quantum_number + the_irrep + '_%s'%str(bb) + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
        
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
                OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                
                plt.errorbar(the_nt_corr_efm, the_mean_efm, yerr = the_sigmas_efm, marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=2, lw=1.5, elinewidth=1.5, zorder=3, capsize=2.75, label = r'$\lambda_{%s}$'%str(bb))
            plt.xlabel(r'$t$', fontsize=14)
            plt.ylabel(r'$a_{t} \;m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=14)
            plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\lambda_{i}(t)$'+ ' [' + the_rs_scheme + ']' + r'$\to t_{0} = %s$'%str(the_t0))
            plt.xticks(the_nt_ticks)
            # plt.ylim([the_max_y,the_min_y])
            if len(the_data)>10: the_n_cols = int(len(the_data)/3)
            else: the_n_cols = int(len(the_data)/2)
            plt.ylim([0.3,1.5])
            plt.legend(fontsize=12, ncol=the_n_cols, handletextpad=0.3)
            plt.tight_layout()
            # plt.show()
            efm_corr_all_fig.savefig(the_location + 'EffectiveMass_ALLEigenvalues' + the_quantum_number + the_irrep + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
        
        ### If the operator analysis was performed, the eigenvalues are also going to be plotted.
        if 'Operators_Analysis' in list(the_matrix_correlator_data[the_irrep].keys()) and the_operators_analysis_flag:
            
            ### If the Operators were chosen by hand, then the plots are also generated.
            if any('Ops_chosen' in the_keys for the_keys in the_matrix_correlator_data[the_irrep +'/Operators_Analysis'].keys()):
                
                ### These are all the operators that the analysis was made on
                the_list_of_chosen_ops = list(filter(lambda x: 'Ops_chosen' in x, the_matrix_correlator_data[the_irrep +'/Operators_Analysis'].keys()))
                
                ### For each of these analysis, it will plot the eigenvalues.
                for the_op_item in the_list_of_chosen_ops:
                    
                    the_data = np.array(the_matrix_correlator_data[the_irrep + '/Operators_Analysis/' + the_op_item + '/t0_%s/Effective_masses/Mean'%the_t0])
                    the_data_sigmas = np.array(the_matrix_correlator_data[the_irrep + '/Operators_Analysis/' + the_op_item +'/t0_%s/Effective_masses/Sigmas'%the_t0])
                    
                    ### This is the information of the selected operators (caption of the plot)
                    the_chosen_ops_string = 'Ops:'
                    for ss in the_op_item[11:]:
                        ### The opertaor written in a better form
                        the_op = the_op_list[int(ss)]
                        OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                        the_chosen_ops_string+=' '+OperatorNamePlot
                    
                    for bb in range(len(the_data)):
                        the_mean_corr = the_data[bb]
                        the_sigmas_corr = np.array(the_data_sigmas[bb])
                        
                         ### Checking for the data
                        the_ymin = vf.CHOOSING_YMIN_PLOT(the_mean_efm)
                    
                        efm_fig = plt.figure()
                        vf.PLOT_CORRELATORS(the_nt_corr_efm, the_mean_corr,  the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t$', r'$a_{t} \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o', NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;\lambda_{%s}$'%str(bb) +r' $(t_{0} = %s$)'%str(the_t0), ymin=the_ymin)
                        
                        efm_fig.text(0.5, 0.01, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                        efm_fig.savefig(the_location + 'EffectiveMass_Eigenvalues' + the_quantum_number + the_op_item +'_' + the_irrep  + '_%s'%str(bb) + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
                        
                    efm_fig = plt.figure()
                    the_ymin = vf.CHOOSING_YMIN_PLOT(the_data[0])
                    the_ymax = vf.CHOOSING_YMAX_PLOT(the_data[0])
                    for bb in range(len(the_data)):
                        the_mean_corr = the_data[bb]
                        the_sigmas_corr = np.array(the_data_sigmas[bb])
                        
                        plt.errorbar(the_nt[the_t0-the_nt[0]:-1], the_mean_corr[the_t0-the_nt[0]:], yerr = the_sigmas_corr[the_t0-the_nt[0]:], marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.1, lw=1.5, elinewidth=1.5, zorder=3, capsize=3.5, label = r'$\lambda_{%s}^{mod}$'%bb)
                    plt.xlabel(r'$t$')
                    plt.ylabel(r'$a_{t} \;m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=14)
                    plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;\lambda_{i} (t_{0} = %s$)'%str(the_t0))
                    plt.tight_layout()
                    plt.legend()
                    plt.xticks(the_nt_ticks)
                    plt.ylim(ymin=the_ymin, ymax=the_ymax)
                    efm_fig.text(0.5, 0.01, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                    # plt.show()
                    efm_fig.savefig(the_location + 'EffectiveMass_ALLEigenvalues' + the_quantum_number + the_op_item +'_' + the_irrep  +  '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
                    
                    
            ### If the Operators were chosen by hand, then the plots are also generated.
            if any('Add_Op' in the_keys for the_keys in the_matrix_correlator_data[the_irrep+'/Operators_Analysis'].keys()):
                ### These are all the operators that the analysis was made on
                the_list_of_chosen_ops = list(filter(lambda x: 'Add_Op' in x, the_matrix_correlator_data[the_irrep +'/Operators_Analysis'].keys()))
                
                ### For each of these analysis, it will plot the eigenvalues.
                for the_op_item in the_list_of_chosen_ops:
                    
                    the_data = np.array(the_matrix_correlator_data[the_irrep + '/Operators_Analysis/' + the_op_item + '/t0_%s/Effective_masses/Mean'%the_t0])
                    the_data_sigmas = np.array(the_matrix_correlator_data[the_irrep + '/Operators_Analysis/' + the_op_item +'/t0_%s/Effective_masses/Sigmas'%the_t0])
                    
                    the_chosen_ops_string = 'Ops: ' + vf.OPERATORS_MH(the_op_list[0].decode('utf-8'))
                    for ss in range(1,int(the_op_item[7:])+1):
                        ### The operator written in a better form
                        the_op = the_op_list[ss]
                        OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                        the_chosen_ops_string+=' '+OperatorNamePlot
                    
                    for bb in range(len(the_data)):
                        the_mean_corr = the_data[bb]
                        the_sigmas_corr = np.array(the_data_sigmas[bb])
                        
                            ### Checking for the data
                        the_ymin = vf.CHOOSING_YMIN_PLOT(the_mean_efm)

                        efm_fig = plt.figure()
                        vf.PLOT_CORRELATORS(the_nt_corr_efm, the_mean_corr,  the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t$', r'$a_{t} \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o', NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;\lambda_{%s}$'%str(bb)+r' $(t_{0} = %s$)'%str(the_t0), ymin=the_ymin)
                        
                        efm_fig.text(0.5, 0.01, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                        efm_fig.savefig(the_location + 'EffectiveMass_Eigenvalues' + the_quantum_number + the_op_item +'_' + the_irrep  + '_%s'%str(bb) + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
                        
                    efm_fig = plt.figure()
                    the_ymin = vf.CHOOSING_YMIN_PLOT(the_data[0])
                    the_ymax = vf.CHOOSING_YMAX_PLOT(the_data[0])
                    for bb in range(len(the_data)):
                        the_mean_corr = the_data[bb]
                        the_sigmas_corr = np.array(the_data_sigmas[bb])
                        
                        plt.errorbar(the_nt[the_t0-the_nt[0]:-1], the_mean_corr[the_t0-the_nt[0]:], yerr = the_sigmas_corr[the_t0-the_nt[0]:], marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.1, lw=1.5, elinewidth=1.5, zorder=3, capsize=3.5, label = r'$\lambda_{%s}^{mod}$'%bb)
                    plt.xlabel(r'$t$')
                    plt.ylabel(r'$a_{t} \;m_{\mathrn{eff}}(t+\frac{1}{2})$', fontsize=14)
                    plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;\lambda_{i} (t_{0} = %s$)'%str(the_t0))
                    plt.tight_layout()
                    plt.legend()
                    plt.xticks(the_nt_ticks)
                    plt.ylim(ymin=the_ymin, ymax=the_ymax)
                    efm_fig.text(0.5, 0.01, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                    # plt.show()
                    efm_fig.savefig(the_location + 'EffectiveMass_ALLEigenvalues' + the_quantum_number + the_op_item +'_' + the_irrep  + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
                
            ### If the Operators were chosen by hand, then the plots are also generated.
            if any('Remove_Op' in the_keys for the_keys in the_matrix_correlator_data[the_irrep+'/Operators_Analysis'].keys()):
                
                ### These are all the operators that the analysis was made on
                the_list_of_chosen_ops = list(filter(lambda x: 'Remove_Op' in x, the_matrix_correlator_data[the_irrep +'/Operators_Analysis'].keys()))
                
                ### For each of these analysis, it will plot the eigenvalues.
                for the_op_item in the_list_of_chosen_ops:
                    
                    the_data = np.array(the_matrix_correlator_data[the_irrep + '/Operators_Analysis/' + the_op_item + '/t0_%s/Effective_masses/Mean'%the_t0])
                    the_data_sigmas = np.array(the_matrix_correlator_data[the_irrep + '/Operators_Analysis/' + the_op_item +'/t0_%s/Effective_masses/Sigmas'%the_t0])
                    
                    ### This is the information of the selected operators (caption of the plot)
                    the_chosen_ops_string = 'Ops: '
                    for ss in range(len(the_op_list)):
                        if ss==int(the_op_item[10:]): continue
                        else:
                            ### The operator written in a better form
                            the_op = the_op_list[ss]
                            OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                            the_chosen_ops_string+=' '+OperatorNamePlot
                    
                    for bb in range(len(the_data)):
                        the_mean_corr = the_data[bb]
                        the_sigmas_corr = np.array(the_data_sigmas[bb])
                        
                        ### Checking for the data
                        the_ymin = vf.CHOOSING_YMIN_PLOT(the_mean_efm)
                            
                        efm_fig = plt.figure()
                        vf.PLOT_CORRELATORS(the_nt_corr_efm, the_mean_corr,  the_sigmas_corr, the_rs_scheme, the_nt_ticks, r'$t$', r'$a_{t} \;m_{\mathrm{eff}}(t+\frac{1}{2})$', 'o', NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;\lambda_{%s}$'%str(bb)+r' $(t_{0} = %s$)'%str(the_t0), ymin=the_ymin)
                        
                        efm_fig.text(0.5, 0.01, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                        efm_fig.savefig(the_location + 'EffectiveMass_Eigenvalues' + the_quantum_number + the_op_item +'_' + the_irrep  + '_%s'%str(bb) + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
                        
                    efm_fig = plt.figure()
                    the_ymin = vf.CHOOSING_YMIN_PLOT(the_data[0])
                    the_ymax = vf.CHOOSING_YMAX_PLOT(the_data[0])
                    for bb in range(len(the_data)):
                        the_mean_corr = the_data[bb]
                        the_sigmas_corr = np.array(the_data_sigmas[bb])
                        the_op = the_op_list[bb]
                        OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                        
                        plt.errorbar(the_nt[the_t0-the_nt[0]:-1], the_mean_corr[the_t0-the_nt[0]:], yerr = the_sigmas_corr[the_t0-the_nt[0]:], marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.1, lw=1.5, elinewidth=1.5, zorder=3, capsize=3.5, label = r'$\lambda_{%s}^{mod}$'%bb)
                    plt.xlabel(r'$t$')
                    plt.ylabel(r'$a_{t} \;m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=14)
                    plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;\lambda_{i} (t_{0} = %s$)'%the_t0)
                    plt.ylim(ymin=the_ymin, ymax=the_ymax)
                    plt.tight_layout()
                    plt.legend()
                    plt.xticks(the_nt_ticks)
                    efm_fig.text(0.5, 0.01, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                    # plt.show()
                    efm_fig.savefig(the_location + 'EffectiveMass_ALLEigenvalues' + the_quantum_number + the_op_item +'_' + the_irrep + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)

        
            
            
def PlotRatioHadronsEffectiveMasses(the_ratio_correlator_data, the_rs_scheme, the_version, the_t0, the_location, the_rebin):
    mr_irreps = list(the_ratio_correlator_data.keys())
    for irrep in mr_irreps:
        the_op_list = list(the_ratio_correlator_data[irrep+'/Operators'])
        the_data = np.array(the_ratio_correlator_data[irrep + '/GEVP/t0_%s/Effective_masses/Mean'%the_t0])
        the_data_sigmas = np.array(the_ratio_correlator_data[irrep + '/GEVP/t0_%s/Effective_masses/Sigmas'%the_t0])
        for bb in range(len(the_op_list)):
            the_mean_corr = the_data[bb][the_t0-2:]
            the_sigmas_corr = the_data_sigmas[bb][the_t0-2:]
            the_nt_corr = np.array(the_ratio_correlator_data[irrep + '/Time_slices'])[the_t0-2:]
            the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
            the_nt_ticks = np.arange(5, the_nt_corr[-1], 5)

            the_op = the_op_list[bb]
            # OperatorNamePlot = vf.OPERATORS_SH(the_op.decode('utf-8'))
            OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
            da_irrep = vf.IrrepInfo(irrep)
            MomentumIrrep = da_irrep.TotalMomPlot
            NameIrrepPlot = da_irrep.NamePlot
            NameIrrep = da_irrep.Name
            
            print('Effective Mass plot in progress...')
            efm_fig = plt.figure()
            plt.errorbar(the_nt, the_mean_corr, yerr = the_sigmas_corr, marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5, label = '%s'%the_rs_scheme)
            plt.xlabel('t')
            plt.ylabel(r'$a_{t} \;m_{eff}(t+\frac{1}{2})$')
            plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0)
            plt.xticks(the_nt_ticks)
            plt.tight_layout()
            #plt.show()
            efm_fig.savefig(the_location + 'EffectiveMass_Eigenvalues_ratios_' + irrep + '_%s'%bb + the_rebin + '_v%s.pdf'%the_version)



### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------


if __name__=="__main__":
    
    myEns = str(sys.argv[1]).upper()
    myWhichCorrelator = str(sys.argv[2]).lower()
    myTypeRs = str(sys.argv[3]).lower()
    myRebinOn = str(sys.argv[4]).lower()
    
    myRb = 2
    myVersion = 'test'
    myT0 = 4
    
    myDataLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH(SAME_THAN_CORRS_SCRIPT_OUTPUT)$/%s/'%myEns)
    
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
        mySingleCorrelatorData = h5py.File(myDataLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/SingleHadrons/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotSingleHadronsEffectiveMasses(mySingleCorrelatorData, myResamplingScheme, myVersion, myPlotLocation, reBin)
        
        mySingleCorrelatorData.close()
    
    elif myWhichCorrelator=='m':
        myMatrixCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/Matrices/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotMultiHadronsEffectiveMasses(myMatrixCorrelatorData, myResamplingScheme, myVersion, myT0, myPlotLocation, reBin)
        
        myMatrixCorrelatorData.close()
        
    elif myWhichCorrelator=='mr':
        myRatioCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/Matrices/Ratios/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotRatioHadronsEffectiveMasses(myRatioCorrelatorData, myResamplingScheme, myVersion, myT0, myPlotLocation, reBin)
        
        myRatioCorrelatorData.close()
    
