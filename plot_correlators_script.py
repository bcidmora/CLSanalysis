import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_functions as vf
  

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
        the_nt_ticks = np.arange(the_nt[0], the_nt[-1], 5)

        the_op_list = list(the_single_correlator_data[irrep+'/Operators'])[0]
        OperatorNamePlot = vf.OPERATORS_SH(the_op_list.decode('utf-8'))
        
        da_irrep = vf.IrrepInfo(irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot
        NameIrrep = da_irrep.Name
        
        print('Correlator plot in progress...')
        the_corr_fig = plt.figure()
        plt.errorbar(the_nt, the_mean_corr, yerr = the_sigmas_corr, marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5, label = the_rs_scheme)
        plt.xlabel('t')
        plt.ylabel(r'$\mathbb{Re}\;C(t)$')
        plt.title( OperatorNamePlot + ' (%s): '%MomentumIrrep + r' $\to$ %s'%NameIrrepPlot)
        plt.xticks(the_nt_ticks)
        plt.legend()
        plt.tight_layout()
        the_corr_fig.savefig(the_location + 'Correlator_' + irrep[:4] +'_%s'%irrep[-1] + the_rebin + '_v%s.pdf'%the_version)
        
        print('Correlator Log-plot in process...')
        the_log_corr_fig = plt.figure()
        plt.errorbar(the_nt, the_mean_corr, yerr = the_sigmas_corr, marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5, label = the_rs_scheme)
        plt.xlabel('t')
        plt.ylabel(r'$\log\mathbb{Re}\;C(t)$')
        plt.title(OperatorNamePlot + ' (%s): '%MomentumIrrep + r' $\to$ %s'%NameIrrepPlot)
        plt.yscale('log')
        plt.legend()
        plt.xticks(the_nt_ticks)
        plt.tight_layout()
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
        
        plt.hist(the_rs, bins=the_nr_bins, label =  r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error))
        plt.vlines(the_mean_rs, 0, the_nr_samples*(the_nr_bins/100)*.65, colors= 'red', label = r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r'$%s$'%the_mean_rs)
        plt.vlines(the_nt_mean, 0, the_nr_samples*(the_nr_bins/100)*.65, colors='black', label = r'$ \bar{C}(t) = $ %s'%the_nt_mean)
        plt.ylim([0,the_nr_samples*(the_nr_bins/100)*.6])
        plt.title(  OperatorNamePlot+ ' (%s): '%MomentumIrrep + ' t = %s'%(tt+the_nt[0]))
        plt.ylabel('Frequency')
        plt.xlabel('Correlator')
        plt.legend()
        plt.tight_layout()
        the_gauss_fig.savefig(the_location + 'Histogram_correlators_' + irrep[:4] + '_' + irrep[-1] + the_rebin + '_v%s.pdf'%the_version)


def PlotMultiHadronCorrelators(the_matrix_correlator_data, the_type_rs, the_version, the_t0, the_location, the_rebin):
    m_irreps = list(the_matrix_correlator_data.keys())
    the_do_eigs = True
    
    if the_type_rs=='jk':
        the_rs_scheme='Jackknife'
    elif the_type_rs=='bt':
        the_rs_scheme='Bootstrap'

    for irrep in m_irreps:
        the_op_list = list(the_matrix_correlator_data[irrep+'/Operators'])
        
        the_data_corr = vf.RESHAPING_CORRELATORS(np.array(the_matrix_correlator_data[irrep + '/Correlators/Real/Mean']))
        the_data_sigmas_corr = np.array(the_matrix_correlator_data[irrep + '/Correlators/Real/Sigmas'])
        
        try:
            the_data = np.array(the_matrix_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Mean'%the_t0])
            the_data_sigmas = np.array(the_matrix_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Covariance_matrix'%the_t0])
            the_do_eigs=True
        except KeyError:
            the_do_eigs=False        
        
        for bb in range(len(the_op_list)):
            the_nt = np.array(the_matrix_correlator_data[irrep + '/Time_slices'])
            the_nt_ticks = np.arange(the_nt[0], the_nt[-1], 5)
            
            the_op = the_op_list[bb]
            OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))
            
            da_irrep = vf.IrrepInfo(irrep)
            MomentumIrrep = da_irrep.TotalMomPlot
            NameIrrepPlot = da_irrep.NamePlot 
            NameIrrep = da_irrep.Name
            
            print('Correlator plots in process...')
                
            corr_fig = plt.figure()
            plt.errorbar(the_nt[the_t0-the_nt[0]:], the_data_corr[bb][bb][the_t0-the_nt[0]:], yerr = the_data_sigmas_corr[bb][the_t0-the_nt[0]:], marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5, label = the_rs_scheme)
            plt.xlabel('t')
            plt.ylabel(r'$\mathbb{Re}\;C(t)$')
            plt.title( NameIrrepPlot+ ' ' + ' (%s): '%MomentumIrrep + r' $\to \;C_{%s}$'%(str(bb)+str(bb)) + '= ' + OperatorNamePlot)
            plt.xticks(the_nt_ticks)
            plt.tight_layout()
            plt.legend()
            # plt.show()
            corr_fig.savefig(the_location + 'DiagonalCorrelator_' + irrep + '_%s'%bb + the_rebin + '_v%s.pdf'%the_version)
            
            print('Correlator Log-plots in progress...')
            corr_fig = plt.figure()
            plt.errorbar(the_nt[the_t0-the_nt[0]:], the_data_corr[bb][bb][the_t0-the_nt[0]:], yerr = the_data_sigmas_corr[bb][the_t0-the_nt[0]:], marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5, label = the_rs_scheme)
            plt.xlabel('t')
            plt.ylabel(r'$\log\mathbb{Re}\;C(t)$')
            plt.title( NameIrrepPlot+ ' ' +  ' (%s): '%MomentumIrrep + r' $\to \;C_{%s}$'%(str(bb)+str(bb)) + '= ' + OperatorNamePlot)
            plt.yscale('log')
            plt.legend()
            plt.tight_layout()
            plt.xticks(the_nt_ticks)
            # plt.show()
            corr_fig.savefig(the_location + 'DiagonalCorrelator_' + irrep  + '_%s_log'%bb + the_rebin + '_v%s.pdf'%the_version)
            
            print('Correlator histogram in progress...')
            tt = int(len(the_nt)/2)+1
            the_gauss_fig = plt.figure()
            the_nt_mean = the_data_corr[bb][bb][tt]
            the_rs = vf.RESHAPING_CORRELATORS_RS_NT(np.array(the_matrix_correlator_data[irrep + '/Correlators/Real/Resampled']))[bb][bb][tt]
            the_nr_samples = np.array(the_matrix_correlator_data[irrep + '/Correlators/Real/Resampled']).shape[1]

            the_mean_rs = np.mean(the_rs)
            the_means_dif = np.abs(the_nt_mean - the_mean_rs)
            the_stat_error = the_data_sigmas_corr[bb][tt]
            
            plt.hist(the_rs, bins=25, label =  r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error))
            plt.vlines(the_mean_rs, 0, 200, colors= 'red', label = r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r' $%s$'%the_mean_rs)
            plt.vlines(the_nt_mean, 0, 200, colors='black', label = r'$ \bar{C}(t) = $ %s'%the_nt_mean)
            plt.title( NameIrrep + ' ' + ' (%s): '%MomentumIrrep +  r'$C_{%s}$'%(str(bb) + str(bb)) + ' [t = %s]'%(tt+the_nt[0]))
            plt.ylabel('Frequency')
            plt.xlabel(r'$Diag(C(t))_{%s}$'%(str(bb)+str(bb)))
            plt.legend()
            plt.tight_layout()
            plt.ylim([0,int(the_nr_samples*.3)])
            # plt.show()
            the_gauss_fig.savefig(the_location + 'Histogram_DiagCorrelator_' + irrep + '_%s'%bb + the_rebin + '_v%s.pdf'%the_version)
                
            if the_do_eigs:
                the_mean_corr = the_data[bb]
                the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                
                print('Eigenvalues plot in process...')
                
                corr_fig = plt.figure()
                plt.errorbar(the_nt[the_t0-the_nt[0]:], the_mean_corr[the_t0-the_nt[0]:], yerr = the_sigmas_corr[the_t0-the_nt[0]:], marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5)
                plt.xlabel('t')
                plt.ylabel(r'$\mathbb{Re}\;C(t)$')
                plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0)
                plt.xticks(the_nt_ticks)
                plt.tight_layout()
                plt.show()
                corr_fig.savefig(the_location + 'Eigenvalues_' + irrep + '_%s'%bb + the_rebin + '_v%s.pdf'%the_version)
                
                print('Eigenvalues Log-plot in progress...')
                corr_fig = plt.figure()
                plt.errorbar(the_nt[the_t0-the_nt[0]:], the_mean_corr[the_t0-the_nt[0]:], yerr = the_sigmas_corr[the_t0-the_nt[0]:], marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5)
                plt.xlabel('t')
                plt.ylabel(r'$\log\mathbb{Re}\;C(t)$')
                plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0)
                plt.yscale('log')
                plt.tight_layout()
                plt.xticks(the_nt_ticks)
                # plt.show()
                corr_fig.savefig(the_location + 'Eigenvalues_' + irrep  + '_%s_log'%bb + the_rebin + '_v%s.pdf'%the_version)
                
                print('Eigenvalues histogram in progress...')
                tt = int(len(the_nt)/2)+1
                the_gauss_fig = plt.figure()
                the_nt_mean = the_mean_corr[tt]
                the_rs = np.array(the_matrix_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Resampled'%the_t0])[bb].transpose()[tt]
                the_nr_samples = np.array(the_matrix_correlator_data[irrep + '/GEVP/t0_%s/Eigenvalues/Resampled'%the_t0]).shape[1]

                the_mean_rs = np.mean(the_rs)
                the_means_dif = np.abs(the_nt_mean - the_mean_rs)
                the_stat_error = the_sigmas_corr[tt]
                
                plt.hist(the_rs, bins=25, label =  r'$\Delta = %s$'%'{:.10e}'.format(the_means_dif) +'\n'+ r'$\sigma = %s$'%'{:.10e}'.format(the_stat_error))
                plt.vlines(the_mean_rs, 0, 200, colors= 'red', label = r'$ \bar{C}_{%s}(t) =$'%the_type_rs + r' $%s$'%the_mean_rs)
                plt.vlines(the_nt_mean, 0, 200, colors='black', label = r'$ \bar{C}(t) = $ %s'%the_nt_mean)
                plt.title( NameIrrep + ' (%s): '%MomentumIrrep +  r'$\lambda_{%s}$'%bb +  ' t = %s'%(tt+the_nt[0]) + r' ($t_{0} = %s$)'%the_t0)
                plt.ylabel('Frequency')
                plt.xlabel(r'Eigenvalue ($\lambda_{%s}$)'%bb)
                plt.legend()
                plt.tight_layout()
                plt.ylim([0,int(the_nr_samples*.3)])
                # plt.show()
                the_gauss_fig.savefig(the_location + 'Histogram_Eigenvalues_' + irrep + '_%s'%bb + the_rebin + '_v%s.pdf'%the_version)
        
        if the_do_eigs:
            the_nt = np.array(the_matrix_correlator_data[irrep + '/Time_slices'])
            the_nt_ticks = np.arange(the_nt[0], the_nt[-1], 5)
                
            corr_fig = plt.figure()
            
            the_markers_list = ['o','v','s','p','^','*','x','d','>','D', '<','8','P','h','1']
            
            for bb in range(len(the_op_list)):
                the_mean_corr = the_data[bb]
                the_sigmas_corr = np.sqrt(np.diag(the_data_sigmas[bb]))
                
                the_op = the_op_list[bb]
                OperatorNamePlot = vf.OPERATORS_MH(the_op.decode('utf-8'))

                da_irrep = vf.IrrepInfo(irrep)
                MomentumIrrep = da_irrep.TotalMomPlot
                NameIrrepPlot = da_irrep.NamePlot 
                NameIrrep = da_irrep.Name
                
                print('Eigenvalues Log-plot in progress...')
                
                plt.errorbar(the_nt[the_t0-the_nt[0]:], the_mean_corr[the_t0-the_nt[0]:], yerr = the_sigmas_corr[the_t0-the_nt[0]:], marker=the_markers_list[bb], ls='None', ms=4.5, markeredgewidth=1.1, lw=1.5, elinewidth=1.5, zorder=3, capsize=3.5, label = r'$\lambda_{%s}$'%bb)
            plt.xlabel('t')
            plt.ylabel(r'$\log\mathbb{Re}\;C(t)$')
            plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' ($t_{0} = %s$)'%the_t0)
            plt.yscale('log')
            plt.tight_layout()
            plt.legend()
            plt.xticks(the_nt_ticks)
            # plt.show()
            corr_fig.savefig(the_location + 'ALLEigenvalues_' + irrep  + '_log' + the_rebin + '_v%s.pdf'%the_version)
    

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
    
    myRb = 2
    myVersion = 'test'
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
    
