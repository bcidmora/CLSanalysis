import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_functions as vf


def PlotSingleHadronsEffectiveMasses(the_single_correlator_data, the_rs_scheme, the_version, the_location, the_rebin):
    s_irreps = list(the_single_correlator_data.keys())
    for irrep in s_irreps:
        the_mean_corr = np.array(the_single_correlator_data[irrep + '/Effective_masses/Mean'])
        the_sigmas_corr = np.array(the_single_correlator_data[irrep + '/Effective_masses/Sigmas'])
        the_nt_corr = np.array(the_single_correlator_data[irrep + '/Time_slices'])
        the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
        the_nt_ticks = np.arange(5, the_nt_corr[-1], 5)

        the_op = list(the_single_correlator_data[irrep+'/Operators'])[0]
        OperatorNamePlot = vf.OPERATORS_SH(the_op.decode('utf-8'))
        da_irrep = vf.IrrepInfo(irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot
        NameIrrep = da_irrep.Name
        
        print('Effective Mass plot in progress...')
        the_efm_fig = plt.figure()
        plt.errorbar(the_nt, the_mean_corr, yerr = the_sigmas_corr, marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5, label='%s'%the_rs_scheme)
        plt.xlabel('t')
        plt.ylabel(r'$a_{t} \;m_{eff}(t+\frac{1}{2})$')
        plt.title( OperatorNamePlot + ' (%s): '%MomentumIrrep + r' $\to$ %s'%NameIrrepPlot)
        plt.xticks(the_nt_ticks)
        plt.legend()
        #plt.ylim([min(the_mean_corr)*0.95, max(the_mean_corr[2:])*1.1])#[the_t0-2:]
        plt.tight_layout()
        #plt.show()
        the_efm_fig.savefig(the_location + 'EffectiveMass_' + irrep[:4] +'_%s'%irrep[-1] + the_rebin + '_v%s.pdf'%the_version)
        
        
def PlotMultiHadronsEffectiveMasses(the_matrix_correlator_data, the_rs_scheme, the_version, the_t0, the_location, the_rebin):
    m_irreps = list(the_matrix_correlator_data.keys())
    the_do_eigs = True
    
    for irrep in m_irreps:
        the_op_list = list(the_matrix_correlator_data[irrep+'/Operators'])
        
        the_data_corr = np.array(the_matrix_correlator_data[irrep + '/Correlators/Real/Effective_masses'])
        the_data_sigmas_corr = np.array(the_matrix_correlator_data[irrep + '/Correlators/Real/Effective_masses_sigmas'])
        
        try:
            the_data = np.array(the_matrix_correlator_data[irrep + '/GEVP/t0_%s/Effective_masses/Mean'%the_t0])
            the_data_sigmas = np.array(the_matrix_correlator_data[irrep + '/GEVP/t0_%s/Effective_masses/Sigmas'%the_t0])
            the_do_eigs = True
        except KeyError:
            the_do_eigs = False
       
        
        for bb in range(len(the_op_list)):
            the_nt_corr = np.array(the_matrix_correlator_data[irrep + '/Time_slices'])[the_t0-2:]
            
            the_mean_corr_corr = the_data_corr[bb]
            the_sigmas_corr_corr = the_data_sigmas_corr[bb]
            the_nt_corr_corr = np.array(the_matrix_correlator_data[irrep + '/Time_slices'])
            the_nt_corr_efm = np.arange(the_nt_corr_corr[0]+0.5, the_nt_corr_corr[-1]+0.5, 1)
            the_nt_ticks_corr = np.arange(3, the_nt_corr[-1], 2)

            the_op = the_op_list[bb]
            
            OperatorNamePlot = vf.OPERATORS_SH(the_op.decode('utf-8'))
            da_irrep = vf.IrrepInfo(irrep)
            MomentumIrrep = da_irrep.TotalMomPlot
            NameIrrepPlot = da_irrep.NamePlot
            NameIrrep = da_irrep.Name
            
            efm_corr_fig = plt.figure()
            plt.errorbar(the_nt_corr_efm, the_mean_corr_corr, yerr = the_sigmas_corr_corr, marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5, label = '%s'%the_rs_scheme)
            plt.xlabel('t')
            plt.ylabel(r'$a_{t} \;m_{eff}(t+\frac{1}{2})$')
            plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;Corr_{%s}$'%(str(bb)+str(bb)) )
            plt.xticks(the_nt_ticks_corr)
            plt.legend()
            plt.tight_layout()
            # plt.show()
            efm_corr_fig.savefig(the_location + 'EffectiveMass_DiagonalCorrelators_' + irrep + '_%s'%bb + the_rebin + '_v%s.pdf'%the_version)
            
            if the_do_eigs:
                the_mean_corr = the_data[bb][the_t0-2:]
                the_sigmas_corr = the_data_sigmas[bb][the_t0-2:]
                the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
                the_nt_ticks = np.arange(5, the_nt_corr[-1], 5)
                
                print('Effective Mass plot in progress...')
                efm_fig = plt.figure()
                plt.errorbar(the_nt, the_mean_corr, yerr = the_sigmas_corr, marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5, label = '%s'%the_rs_scheme)
                plt.xlabel('t')
                plt.ylabel(r'$a_{t} \;m_{eff}(t+\frac{1}{2})$')
                plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0)
                plt.xticks(the_nt_ticks)
                plt.legend()
                plt.tight_layout()
                # plt.show()
                efm_fig.savefig(the_location + 'EffectiveMass_Eigenvalues_' + irrep + '_%s'%bb + the_rebin + '_v%s.pdf'%the_version)
            
            
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
            OperatorNamePlot = vf.OPERATORS_SH(the_op.decode('utf-8'))
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
    
