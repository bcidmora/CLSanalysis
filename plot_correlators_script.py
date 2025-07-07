import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_functions as vf
import warnings
warnings.filterwarnings('ignore')
  

def PlotSingleHadronCorrelators(the_single_correlator_data, the_type_rs, the_version, the_location, the_rebin):
    s_irreps = list(the_single_correlator_data.keys())
    
    if the_type_rs=='jk':
        the_rs_scheme='Jackknife'
    elif the_type_rs=='bt':
        the_rs_scheme='Bootstrap'   
    
    the_nr_bins = 25
    for irrep in s_irreps:
        the_mean_corr = np.array(the_single_correlator_data[irrep + '/Correlators/Real/Mean'])
        the_sigmas_corr = np.array(the_single_correlator_data[irrep + '/Correlators/Real/Sigmas'])
        the_nt = np.array(the_single_correlator_data[irrep + '/Time_slices'])
        the_nt_ticks = np.arange(the_nt[0]+1, the_nt[-1], int(len(the_nt)/5))

        the_op_list = list(the_single_correlator_data[irrep+'/Operators'])[0]
        OperatorNamePlot = vf.OPERATORS_SH(the_op_list.decode('utf-8'))
        
        da_irrep = vf.IrrepInfo(irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot
        NameIrrep = da_irrep.Name
        
        print('Correlator plot in progress...')
        the_corr_fig = plt.figure()
        vf.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, 't', r'$\mathbb{Re}\;C(t)$', 'o', OperatorNamePlot + ' (%s): '%MomentumIrrep + r' $\to$ %s'%NameIrrepPlot)
        # plt.show()
        the_corr_fig.savefig(the_location + 'Correlator_' + irrep[:4] +'_%s'%irrep[-1] + the_rebin + '_v%s.pdf'%the_version)
        
        print('Correlator Log-plot in process...')
        the_log_corr_fig = plt.figure()
        vf.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, 't', r'$\mathbb{Re}\;C(t)$', 'o', OperatorNamePlot + ' (%s): '%MomentumIrrep + r' $\to$ %s'%NameIrrepPlot, yscale='log')
        # plt.show()
        the_log_corr_fig.savefig(the_location + 'Correlator_' + irrep[:4] + '_%s'%irrep[-1] + '_log' + the_rebin + '_v%s.pdf'%the_version)
        
        print('Correlator histogram in process...')
        tt = int(len(the_nt)/2)+1
        the_gauss_fig = plt.figure()
        the_nt_mean = the_mean_corr[tt]
        the_rs = np.array(the_single_correlator_data[irrep + '/Correlators/Real/Resampled'])[tt]
        the_nr_samples = np.array(the_single_correlator_data[irrep + '/Correlators/Real/Resampled']).shape[1]

        the_mean_rs = np.mean(the_rs)
        the_means_dif = np.abs(the_nt_mean - the_mean_rs)
        the_stat_error = the_sigmas_corr[tt]
        
        vf.PLOT_HISTOGRAMS(the_rs, r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error), the_mean_rs, r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r' $%s$'%the_mean_rs, the_nt_mean, r'$ \bar{C}(t) = $ %s'%the_nt_mean, OperatorNamePlot+ ' (%s): '%MomentumIrrep + ' t = %s'%(tt+the_nt[0]), the_nr_bins,  'Correlator')
        # plt.show()
        the_gauss_fig.savefig(the_location + 'Histogram_correlators_' + irrep[:4] + '_' + irrep[-1] + the_rebin + '_v%s.pdf'%the_version)
        
        

def PlotMultiHadronCorrelators(the_matrix_correlator_data, the_type_rs, the_version, the_t0, the_location, the_rebin, **kwargs):
    
    ### Getting all the irreps in this ensemble
    m_irreps = list(the_matrix_correlator_data.keys())
    
    ### If not all the irreps are wanted t be plotted
    if kwargs.get('nr_irreps')!=None:
        the_first_irrep = 0
        the_last_irrep = int(kwargs.get('nr_irreps'))
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
        
    ### Choosing the resampling scheme to put in the legend
    if the_type_rs=='jk':
        the_rs_scheme='Jackknife'
    elif the_type_rs=='bt':
        the_rs_scheme='Bootstrap'

    ### All different types of amrkers to plot more than one dataset at the same time
    the_markers_list = ['o','v','s','p','^','*','x','d','>','D', '<','8','P','h','1']
    
    ### This is the nr. of bins to plot the histograms
    the_nr_bins = 25
    
    ### Loop over the irreducible representations
    for irrep in m_irreps:
        ### The list of operators of this irrep
        the_op_list = list(the_matrix_correlator_data[irrep+'/Operators'])
        
        ### The correaltor dataset
        the_data_corr = vf.RESHAPING_CORRELATORS(np.array(the_matrix_correlator_data[irrep + '/Correlators/Real/Mean']))
        
        ### The sigmas of this correlator dataset
        the_data_sigmas_corr = np.array(the_matrix_correlator_data[irrep + '/Correlators/Real/Sigmas'])
        
        ### This is the time interval
        the_nt = np.array(the_matrix_correlator_data[irrep + '/Time_slices'])
        
        ### These are going to be the ticks in the x-label of the plots
        the_nt_ticks = np.arange(the_nt[0]+1, the_nt[-1], int(len(the_nt)/5))
        
        ### Plotting eigenvalues in the full time range or only starting from t0
        if kwargs.get('full_range_nt')==None: the_start_nt = the_t0-the_nt[0]
        else: the_start_nt = the_nt[0]
        
        ### Information of the irrep to write it properly in the plots.
        da_irrep = vf.IrrepInfo(irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot 
        NameIrrep = da_irrep.Name
        
        ### Loop over the number of operators of this matrix
        for bb in range(len(the_op_list)):      
            
            ### The specific operator
            the_op = the_op_list[bb]
            OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
            
            print('Correlator plots in process...')
                
            ### Plotting the diagonal correlators
            corr_fig = plt.figure()
            vf.PLOT_CORRELATORS(the_nt, the_data_corr[bb][bb], the_data_sigmas_corr[bb], the_rs_scheme, the_nt_ticks, 't', r'$\mathbb{Re}\;C(t)$', 'o', NameIrrepPlot + ' (%s) '%MomentumIrrep + r' $\to \;C_{%s}$'%(str(bb)+str(bb)) + '= ' + OperatorNamePlot)
            corr_fig.savefig(the_location + 'DiagonalCorrelator_' + irrep + '_%s'%str(bb) + the_rebin + '_v%s.pdf'%the_version)
            
            ### Plotting the log of the diagonal correlators.
            print('Correlator Log-plots in progress...')
            corr_fig = plt.figure()
            vf.PLOT_CORRELATORS(the_nt, the_data_corr[bb][bb], the_data_sigmas_corr[bb], the_rs_scheme, the_nt_ticks, 't', r'$\log\mathbb{Re}\;C(t)$', 'o', NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to \;C_{%s}$'%(str(bb)+str(bb)) + '= ' + OperatorNamePlot, yscale='log')
            corr_fig.savefig(the_location + 'DiagonalCorrelator_' + irrep  + '_%s_log'%str(bb) + the_rebin + '_v%s.pdf'%the_version)
            
            ### Plotting the histogram at a certain time slice t
            print('Correlator histogram in progress...')
            tt = int(len(the_nt)/2)+1
            the_gauss_fig = plt.figure()
            the_nt_mean = the_data_corr[bb][bb][tt]
            the_rs = vf.RESHAPING_CORRELATORS_RS_NT(np.array(the_matrix_correlator_data[irrep + '/Correlators/Real/Resampled']))[bb][bb][tt]
            the_nr_samples = np.array(the_matrix_correlator_data[irrep + '/Correlators/Real/Resampled']).shape[1]
            
            ### Here the mean value of the sampling data and the mean value are compared to check the quality of the resampled data
            the_mean_rs = np.mean(the_rs)
            the_means_dif = np.abs(the_nt_mean - the_mean_rs)
            the_stat_error = the_data_sigmas_corr[bb][tt]
            
            ### Plotting the histogram now
            vf.PLOT_HISTOGRAMS(the_rs, r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error), the_mean_rs, r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r' $%s$'%the_mean_rs, the_nt_mean, r'$ \bar{C}(t) = $ %s'%the_nt_mean, NameIrrepPlot + ' (%s) '%MomentumIrrep +  r'$\to \;C_{%s}$'%(str(bb) + str(bb)) + ' (t = %s)'%(tt+the_nt[0]) + ' ' + OperatorNamePlot, the_nr_bins, r'$Diag(C(t))_{%s}$'%(str(bb)+str(bb)))
            # plt.show()
            the_gauss_fig.savefig(the_location + 'Histogram_DiagCorrelator_' + irrep + '_%s'%str(bb) + the_rebin + '_v%s.pdf'%the_version)
        
        ### The Diagonal of the correlators are plotted all together with their errors to compare them directly. 
        corr_fig = plt.figure()
        print('ALL Correlators Log-plot in progress...')
        ### Loop over each of the entries of the diagonal of the correlation matrix
        for bb in range(len(the_op_list)):
            
            ### Name of this operator
            the_op = the_op_list[bb]
            OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))           
            
            plt.errorbar(the_nt, the_data_corr[bb][bb], the_data_sigmas_corr[bb],  marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=3.5, label = r'$C_{%s}$ = '%(str(bb)+str(bb)) + OperatorNamePlot)
        plt.xlabel('t')
        plt.ylabel(r'$\log\mathbb{Re}\;C(t)$')
        plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;C_{ii}$(t)')
        plt.yscale('log')
        plt.tight_layout()
        plt.legend()
        plt.xticks(the_nt_ticks)
        # plt.show()
        corr_fig.savefig(the_location + 'ALLDiagonalCorrelators_' + irrep  + '_log' + the_rebin + '_v%s.pdf'%the_version)
            
        ### Here the Eigenvalues are plotted all together too in a log-plot
        if 'GEVP' in list(the_matrix_correlator_data[irrep].keys()):
            
            the_data = np.array(the_matrix_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Mean'%the_t0])
            the_data_sigmas = np.array(the_matrix_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Covariance_matrix'%the_t0])
            
            for bb in range(len(the_data)):
                ### The mean value of the eigenvalue_{i}
                the_mean_corr = the_data[bb]
                
                ### The corresponding sigmas of this eigenvalue
                the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                
                print("..................................................\n")
                print('Eigenvalue = %s plot in process...'%str(bb))
                
                corr_fig = plt.figure()
                ### Plotting the eigenvalues one by one
                vf.PLOT_CORRELATORS(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], the_sigmas_corr[the_start_nt:], the_rs_scheme, the_nt_ticks, 't', r'$\lambda_{i}(t)$', 'o',  NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %s$)'%str(the_t0))
                corr_fig.savefig(the_location + 'Eigenvalues_' + irrep + '_%s'%str(bb)  +'_t0_%s'%str(the_t0)+ the_rebin + '_v%s.pdf'%the_version)
                
                ### Plotting the eigenvalues log-plots one by one
                print('Eigenvalue = %s Log-plot in progress...'%str(bb))
                corr_fig = plt.figure()
                vf.PLOT_CORRELATORS(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], the_sigmas_corr[the_start_nt:], the_rs_scheme, the_nt_ticks, 't', r'$\log\,(\lambda_{i}(t))$', 'o',  NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %s$)'%str(the_t0), yscale='log')
                corr_fig.savefig(the_location + 'Eigenvalues_' + irrep  + '_%s_log'%str(bb) + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
                
                ### Plotting the histograms of these eigenvalues
                print('Eigenvalue = %s histogram in progress...'%str(bb))
                tt = int(len(the_nt)/2)+1
                the_gauss_fig = plt.figure()
                the_nt_mean = the_mean_corr[tt]
                the_rs = np.array(the_matrix_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Resampled'%the_t0])[bb].transpose()[tt]
                the_nr_samples = np.array(the_matrix_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Resampled'%the_t0]).shape[1]

                the_mean_rs = np.mean(the_rs)
                the_means_dif = np.abs(the_nt_mean - the_mean_rs)
                the_stat_error = the_sigmas_corr[tt]
                
                vf.PLOT_HISTOGRAMS(the_rs, r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error), the_mean_rs, r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r' $%s$'%the_mean_rs, the_nt_mean, r'$ \bar{C}(t) = $ %s'%the_nt_mean, NameIrrepPlot + ' (%s) '%MomentumIrrep +  r'$\to\; \lambda_{%s}$'%str(bb) + r'($t_{0} = %s$)'%str(the_t0) +  ' [t = %s]'%(tt+the_nt[0]), the_nr_bins, r'Eigenvalue ($\lambda_{%s}$)'%str(bb))
                # plt.show()
                the_gauss_fig.savefig(the_location + 'Histogram_Eigenvalues_' + irrep + '_%s'%str(bb) + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
            
            print("..................................................\n")
            
            print('ALL Eigenvalues Log-plot in progress...')
            corr_fig = plt.figure()           
            for bb in range(len(the_data)):
                
                the_mean_corr = the_data[bb]
                the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                
                plt.errorbar(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], yerr = the_sigmas_corr[the_start_nt:], marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=3.5, label = r'$\lambda_{%s}$'%str(bb))
            plt.xlabel('t')
            plt.ylabel(r'$\log\,(\lambda_{i}(t))$')
            plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;\lambda_{i} (t_{0} = %s$)'%str(the_t0))
            plt.yscale('log')
            plt.tight_layout()
            plt.legend()
            plt.xticks(the_nt_ticks)
            # plt.show()
            corr_fig.savefig(the_location + 'ALLEigenvalues_' + irrep  + '_log' + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
        
        ### If the Operator Analysis was performed, then the modifiedf eigenvalues are also going to be plotted
        if 'Operators_Analysis' in list(the_matrix_correlator_data[irrep].keys()):
            
            if any('Ops_chosen_' in the_keys for the_keys in the_matrix_correlator_data[irrep+'/Operators_Analysis'].keys()):
                
                ### Extracting all the elements that contain information about a reduced dataset
                the_list_of_chosen_ops = list(filter(lambda x: 'Ops_chosen_' in x, the_matrix_correlator_data[irrep+'/Operators_Analysis'].keys()))
                
                print("..................................................\n")
                print('Modified eigenvalues chosen operator plots in progress...')
                
                ### Loop over those elements
                for the_op_item in the_list_of_chosen_ops:
                    
                    the_data = np.array(the_matrix_correlator_data[irrep + '/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Mean'%the_t0])
                    the_data_sigmas = np.array(the_matrix_correlator_data[irrep + '/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Covariance_matrix'%the_t0])
                    
                    ### This is the information of the selected operators (caption of the plot)
                    the_chosen_ops_string = 'Ops:'
                    for ss in the_op_item[11:]:
                        ### The opertaor written in a better form
                        the_op = the_op_list[int(ss)]
                        OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                        the_chosen_ops_string+=' '+OperatorNamePlot
                        
                    ### Loop over the modified eigenvalues with the new dataset
                    for bb in range(len(the_data)):
                        ### This is the y-axis dataset
                        the_mean_corr = the_data[bb]
                        
                        ### These are the sigmas of that data
                        the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                        
                        corr_fig = plt.figure()
                        vf.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, 't', r'$\log\,(\lambda_{i}(t))$', 'o',  NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %s$)'%str(the_t0), yscale='log')
                        
                        ### This is the text in the box
                        corr_fig.text(0.5, 0.005, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                        
                        corr_fig.savefig(the_location + 'Eigenvalues_' + the_op_item + '_' + irrep  + '_%s_log'%str(bb) + '_t0_%s'%str(the_t0)+ the_rebin + '_v%s.pdf'%the_version)
                    
                    print("..................................................\n")
                    print('ALL Modified Eigenvalues Log-plot in progress...')
                    ### Here we plot all the modified eigenvalues together
                    corr_fig = plt.figure()
                    for bb in range(len(the_data)):
                        
                        the_mean_corr = the_data[bb]
                        the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                        
                        plt.errorbar(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], yerr = the_sigmas_corr[the_start_nt:], marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=3.5, label = r'$\lambda_{%s}^{mod}$'%str(bb))
                    plt.xlabel('t')
                    plt.ylabel(r'$\log\,(\lambda_{i}(t))$')
                    plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;\lambda_{i} (t_{0} = %s$)'%the_t0)
                    plt.yscale('log')
                    plt.tight_layout()
                    plt.legend()
                    plt.xticks(the_nt_ticks)
                    # plt.show()
                    corr_fig.text(0.5, 0.005, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                    corr_fig.savefig(the_location + 'ALLEigenvalues_' + the_op_item + '_' + irrep  + '_log' + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
            
            if any('Add_Op' in the_keys for the_keys in the_matrix_correlator_data[irrep+'/Operators_Analysis'].keys()):
                
                the_list_of_chosen_ops = list(filter(lambda x: 'Add_Op' in x, the_matrix_correlator_data[irrep+'/Operators_Analysis'].keys()))
                
                print("..................................................\n")
                print('Modified eigenvalues plots in progress...')
                
                for the_op_item in the_list_of_chosen_ops:
                    
                    the_data = np.array(the_matrix_correlator_data[irrep + '/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Mean'%the_t0])
                    the_data_sigmas = np.array(the_matrix_correlator_data[irrep + '/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Covariance_matrix'%the_t0])
                    
                    ### This is the information of the selected operators (caption of the plot)
                    the_chosen_ops_string = 'Ops: ' + vf.OPERATORS_MH(the_op_list[0].decode('utf-8'))
                    for ss in range(1,int(the_op_item[7:])+1):
                        ### The operator written in a better form
                        the_op = the_op_list[ss]
                        OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                        the_chosen_ops_string+=' '+OperatorNamePlot
                    
                    ### Loop over the modified eigenvalues with the new dataset
                    for bb in range(len(the_data)):
                        the_mean_corr = the_data[bb]
                        the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                        
                        corr_fig = plt.figure()
                        vf.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, 't', r'$\log\,(\lambda_{i}(t))$', 'o',  NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %s$)'%str(the_t0), yscale='log')
                        
                        corr_fig.text(0.5, 0.007, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                        corr_fig.savefig(the_location + 'Eigenvalues_' + the_op_item + '_' + irrep  + '_%s_log'%str(bb) + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
                    
                    print("..................................................\n")
                    print('ALL Modified Eigenvalues Log-plot in progress...')
                    corr_fig = plt.figure()
                    for bb in range(len(the_data)):
                        
                        the_mean_corr = the_data[bb]
                        the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                        
                        plt.errorbar(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], yerr = the_sigmas_corr[the_start_nt:], marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=3.5, label = r'$\lambda_{%s}^{mod}$'%bb)
                    plt.xlabel('t')
                    plt.ylabel(r'$\log\,(\lambda_{i}(t))$')
                    plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;\lambda_{i} (t_{0} = %s$)'%str(the_t0))
                    plt.yscale('log')
                    plt.tight_layout()
                    plt.legend()
                    plt.xticks(the_nt_ticks)
                    # plt.show()
                    corr_fig.text(0.5, 0.005, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                    corr_fig.savefig(the_location + 'ALLEigenvalues_' + the_op_item + '_' + irrep  + '_log' + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
            
            if any('Remove_Op' in the_keys for the_keys in the_matrix_correlator_data[irrep+'/Operators_Analysis'].keys()):
                
                the_list_of_chosen_ops = list(filter(lambda x: 'Remove_Op' in x, the_matrix_correlator_data[irrep+'/Operators_Analysis'].keys()))
                
                print("..................................................\n")
                print('Modified eigenvalues plots in progress...')
                
                
                for the_op_item in the_list_of_chosen_ops:
                    
                    the_data = np.array(the_matrix_correlator_data[irrep + '/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Mean'%the_t0])
                    the_data_sigmas = np.array(the_matrix_correlator_data[irrep + '/Operators_Analysis/'+the_op_item+'/t0_%s/Eigenvalues/Covariance_matrix'%the_t0])
                    
                    ### This is the information of the selected operators (caption of the plot)
                    the_chosen_ops_string = 'Ops: '
                    for ss in range(len(the_op_list)):
                        if ss==int(the_op_item[10:]): continue
                        else:
                            ### The operator written in a better form
                            the_op = the_op_list[ss]
                            OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
                            the_chosen_ops_string+=' '+OperatorNamePlot
                            
                    ### Loop over the modified eigenvalues with the new dataset
                    for bb in range(len(the_data)):
                        the_mean_corr = the_data[bb]
                        the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                        
                        corr_fig = plt.figure()
                        vf.PLOT_CORRELATORS(the_nt, the_mean_corr, the_sigmas_corr, the_rs_scheme, the_nt_ticks, 't', r'$\log\,(\lambda_{i}(t))$', 'o',  NameIrrepPlot+ ' (%s) '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %s$)'%str(the_t0), yscale='log')
                        
                        corr_fig.text(0.5, 0.007, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                        corr_fig.savefig(the_location + 'Eigenvalues_' + the_op_item + '_' + irrep  + '_%s_log'%str(bb) + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)
                        
                    print("..................................................\n")
                    print('ALL Modified Eigenvalues Log-plot in progress...')
                    
                    corr_fig = plt.figure()
                    for bb in range(len(the_data)):
                        the_mean_corr = the_data[bb]
                        the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                        
                        plt.errorbar(the_nt[the_start_nt:], the_mean_corr[the_start_nt:], yerr = the_sigmas_corr[the_start_nt:], marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=3.5, label = r'$\lambda_{%s}^{mod}$'%str(bb))
                    plt.xlabel('t')
                    plt.ylabel(r'$\log\,(\lambda_{i}(t))$')
                    plt.title( NameIrrepPlot+ ' (%s) '%MomentumIrrep + r'$\to\;\lambda_{i} (t_{0} = %s$)'%str(the_t0))
                    plt.yscale('log')
                    plt.tight_layout()
                    plt.legend()
                    plt.xticks(the_nt_ticks)
                    # plt.show()
                    corr_fig.text(0.5, 0.007, the_chosen_ops_string, ha='center', va='bottom', fontsize=10)
                    corr_fig.savefig(the_location + 'ALLEigenvalues_' + the_op_item + '_' + irrep  + '_log' + '_t0_%s'%str(the_t0) + the_rebin + '_v%s.pdf'%the_version)


def PlotRatioHadronCorrelators(the_ratio_correlator_data, the_type_rs, the_version, the_t0, the_location):
    mr_irreps = list(the_ratio_correlator_data.keys())
    the_nr_bins = 25
    for irrep in mr_irreps:
            the_op_list = list(the_ratio_correlator_data[irrep+'/Operators'])
            the_data = np.array(the_ratio_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Mean'%the_t0])
            the_data_sigmas = np.array(the_ratio_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Covariance_matrix'%the_t0])
            for bb in range(len(the_op_list)):
                the_mean_corr = the_data[bb]
                the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                the_nt = np.array(the_ratio_correlator_data[irrep + '/Time_slices'])
                the_nt_ticks = np.arange(the_nt[0], the_nt[-1], 5)

                the_op = the_op_list[bb]
                OperatorNamePlot = vf.OPERATORS_SH(the_op.decode('utf-8'))
                da_irrep = vf.IrrepInfo(irrep)
                
                MomentumIrrep = da_irrep.TotalMomPlot
                NameIrrepPlot = da_irrep.NamePlot
                NameIrrep = da_irrep.Name
                
                print('Correlator plot in process...')
                the_corr_fig = plt.figure()
                plt.errorbar(the_nt[the_t0-the_nt[0]:], the_mean_corr[the_t0-the_nt[0]:], yerr = the_sigmas_corr[the_t0-the_nt[0]:], marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5)
                plt.xlabel('t')
                plt.ylabel(r'$\mathbb{Re}\;C(t)$')
                plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0)
                plt.tight_layout()
                plt.xticks(the_nt_ticks)
                #plt.show()
                the_corr_fig.savefig(the_location + 'Ratios/Eigenvalues_ratios_' + irrep  + '_%s_'%bb + the_rebin + 'v%s.pdf'%the_version)
                
                print('Correlator Log-plot in process...')
                the_log_corr_fig = plt.figure()
                plt.errorbar(the_nt[the_t0-the_nt[0]:], the_mean_corr[the_t0-the_nt[0]:], yerr = the_sigmas_corr[the_t0-the_nt[0]:], marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5)
                plt.xlabel('t')
                plt.ylabel(r'$\log\mathbb{Re}\;C(t)$')
                plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0)
                plt.yscale('log')
                plt.tight_layout()
                plt.xticks(the_nt_ticks)
                #plt.show()
                the_log_corr_fig.savefig(the_location + 'Ratios/Eigenvalues_ratios_' + irrep + '_%s_log'%bb + the_rebin + '_v%s.pdf'%the_version)
                
                print('Correlator histogram in progress...')
                tt = int(len(the_nt)/2)+1
                the_gauss_fig = plt.figure()
                the_mean = the_mean_corr[tt]
                the_rs = np.array(mr[irrep + '/GEVP/t0_%s/Eigenvalues/Resampled'%the_t0])[bb].transpose()[tt]

                the_mean_rs = np.mean(the_rs)
                the_means_dif = np.abs(the_mean - the_mean_rs)
                the_stat_error = the_sigmas_corr[tt]
                the_nr_samples = np.array(mr[irrep + '/GEVP/t0_%s/Eigenvalues/Resampled'%the_t0])[bb].shape[0]
                
                plt.hist(the_rs, bins=the_nr_bins, label =  r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error))
                plt.vlines(the_mean_rs, 0,  the_nr_samples*(the_nr_bins/100)*.65, colors= 'red', label = r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r' $%s$'%the_mean_rs)
                plt.vlines(the_mean, 0,  the_nr_samples*(the_nr_bins/100)*.65, colors='black', label = r'$ \bar{C}(t) = $ %s'%the_mean)
                plt.title( NameIrrep + ' (%s): '%MomentumIrrep +  r'$\lambda_{%s}$'%bb +  ' t = %s'%(tt+the_nt[0]) + r' ($t_{0} = %s$)'%the_t0)
                plt.ylabel('Frequency')
                plt.xlabel(r'Eigenvalue ($\lambda_{%s}$)'%bb)
                plt.legend()
                plt.tight_layout()
                plt.ylim([0,the_nr_samples*(the_nr_bins/100)*.6])
                # plt.show()
                the_gauss_fig.savefig(the_location + 'Histogram_Eigenvalues_ratios_' + irrep + '_%s_'%bb + the_rebin + 'v%s.pdf'%the_version)



### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------




if __name__=="__main__":
    
    myEns = str(sys.argv[1]).upper()
    myWhichCorrelator = str(sys.argv[2]).lower()
    myTypeRs = str(sys.argv[3]).lower()
    myRebinOn = str(sys.argv[4]).lower()
    
    myRb = 1
    myVersion = '_test'
    myT0 = 4 
    
    myDataLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH(SAME_THAN_CORRS_SCRIPT_OUTPUT)$/%s/'%myEns)
    
    vf.INFO_PRINTING(myWhichCorrelator, myEns)
    
    if myRebinOn=='rb': 
        rb = int(myRb)
        reBin = '_bin'+str(rb)
    else:
        reBin = '' 
        
    if myTypeRs=='jk':
        myResamplingScheme='Jackknife'
    elif myTypeRs=='bt':
        myResamplingScheme='Bootstrap'    

    if myWhichCorrelator=='s':
        mySingleCorrelatorData = h5py.File(myDataLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/SingleHadrons/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotSingleHadronCorrelators(mySingleCorrelatorData, myTypeRs, myVersion, myPlotLocation, reBin)
        
        mySingleCorrelatorData.close()
    
    elif myWhichCorrelator=='m':
        myMatrixCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/Matrices/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotMultiHadronCorrelators(myMatrixCorrelatorData, myTypeRs, myVersion, myT0, myPlotLocation, reBin)
        
        myMatrixCorrelatorData.close()
        
    elif myWhichCorrelator=='mr':
        myRatioCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/Matrices_Ratios/'%myEns +  '%s/'%myResamplingScheme)
        
        myRatioCorrelatorData.close()
    
