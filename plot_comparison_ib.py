import numpy as np
import matplotlib.pyplot as plt
import math

### LAYOUT AND EXTRA FUNCTIONS
import set_of_layout_functions as vfl
import set_of_analysis_functions as vfa
import set_of_plot_functions as vfp

### EXTRA LIBRARIES
import sys
import os
import h5py

import ensembles as ed


def PlotEffectiveMassComparison(the_single_correlator_data_1, the_single_correlator_data_2, the_rs_scheme, the_location, the_label_1, the_label_2, **kwargs):
    
    s_irreps = list(the_single_correlator_data_1.keys())
    
    ### Loop over the SH irreps (all of them)
    for the_irrep in s_irreps:
        ### Central values of the correlators
        the_mean_corr_1 = np.asarray(the_single_correlator_data_1[f'{the_irrep}/Effective_masses/Mean'])
        the_mean_corr_2 = np.asarray(the_single_correlator_data_2[f'{the_irrep}/Effective_masses/Mean'])
        
        ### Statistical errors of the correlators
        the_sigmas_corr_1 = np.asarray(the_single_correlator_data_1[f'{the_irrep}/Effective_masses/Sigmas'])
        the_sigmas_corr_2 = np.asarray(the_single_correlator_data_2[f'{the_irrep}/Effective_masses/Sigmas'])
        
        ##3 Time extent
        the_nt_corr = np.array(the_single_correlator_data_1[f'{the_irrep}/Time_slices'])
        the_nt = np.arange(the_nt_corr[0]+0.5, the_nt_corr[-1]+0.5, 1)
        
        the_nt_ticks = np.arange(the_nt_corr[2], the_nt_corr[-1], 2)
        
        ### The SH operator that appears in the plot
        the_op = list(the_single_correlator_data_1[f'{the_irrep}/Operators'])[0]
        theOperatorNamePlot = vfp.OPERATORS_SH(the_op.decode('utf-8'))
        
        da_irrep = vfp.IrrepInfo(the_irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        nrMomentumIrrep = da_irrep.Momentum
        NameIrrepPlot = da_irrep.NamePlot 
    
        OperatorNamePlot = vfp.SH_OPERATORS_RELABEL(theOperatorNamePlot, NameIrrepPlot, nrMomentumIrrep )
        
        the_efm_fig = plt.figure()
        plt.errorbar(the_nt, the_mean_corr_1, yerr=the_sigmas_corr_1, marker=vfp.the_markers_list[0], ls='None', ms=6, markeredgewidth=1.25, lw=1.25, elinewidth=1.25, capsize=5, label=the_label_1, color=vfp.the_colors[0])
        plt.errorbar(the_nt, the_mean_corr_2, yerr=the_sigmas_corr_2, marker=vfp.the_markers_list[1], ls='None', ms=6, markeredgewidth=1.25, lw=1.25, elinewidth=1.25, capsize=5, label=the_label_2, color=vfp.the_colors[1])
        plt.xlabel(r'$t\,/\, a$', fontsize=26)
        plt.ylabel(r'$a_{t}\,m_{\mathrm{eff}}(t+\frac{1}{2})$', fontsize=26)
        plt.title(f'{OperatorNamePlot} ({the_rs_scheme})', fontsize=22)
        plt.xticks(the_nt_ticks, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlim([the_nt[0] - 1, the_nt_ticks[-1] + 1])
        plt.legend(fontsize=18, handletextpad=0.3)
        plt.tight_layout()
        plt.show()
        the_efm_fig.savefig(f'{the_location}Comparison_EffectiveMass_{the_irrep}.pdf', bbox_inches='tight')

def PlotFitsComparison(the_ensemble_fits, the_selection, the_observable_name, the_location):
    the_nr_ensembles = len(the_selection)
    the_nr_cols = 4
    the_nr_rows = math.ceil(the_nr_ensembles / the_nr_cols)
    fig, axes = plt.subplots(the_nr_rows, the_nr_cols, figsize=(3 * the_nr_cols, 3 * the_nr_rows))
    
    axes = np.atleast_1d(axes).flatten()
    for i, (ax, ensemble) in enumerate(zip(axes, the_selection)):
        the_data = the_ensemble_fits[ensemble]

        alex_val = the_data['Alex']['Value']
        alex_err = the_data['Alex']['Error']

        barb_val = the_data['Barb']['Value']
        barb_err = the_data['Barb']['Error']

        x = np.arange(2)
        ax.errorbar( x[0], barb_val, yerr=barb_err, fmt=vfp.the_markers_list[0],  ls='None', ms=6, markeredgewidth=1.25, lw=1.25, elinewidth=1.25, capsize=5, label='Barb', color=vfp.the_colors[0])
        ax.errorbar(x[1], alex_val, yerr=alex_err, fmt=vfp.the_markers_list[1], ls='None', ms=6, markeredgewidth=1.25, lw=1.25, elinewidth=1.25, capsize=5, label='Alex', color=vfp.the_colors[1])
        
        ax.set_xticks(x, fontsize=12)
        ax.set_yticks([barb_val, alex_val])
        ax.tick_params(axis='y', labelsize=12)
        ax.set_xticklabels(['Barb', 'Alex'])
        ax.tick_params(axis='x', labelsize=12)
        ax.set_xlim(-0.75, 1.75)
        
        if i % the_nr_cols == 0:
            ax.set_ylabel(the_observable_name, fontsize=20)
        ax.set_title(ensemble, fontsize=14)
        ax.grid(True, alpha=0.3)
    for i in range(the_nr_ensembles, len(axes)):
        fig.delaxes(axes[i])
    fig.suptitle(fr'Comparison $\Omega$ mass (2-exp fit, Bootstrap)', fontsize=20)
    fig.tight_layout()
    fig.savefig(f'{the_location}EnsemblesComparisonFits_OmegaMass.pdf', bbox_inches='tight')
    plt.show()
    

if __name__=="__main__":
    
    myEffMass = False
    myFitsComp = True
    
    # myEnsemble = ['A654']
    myEnsemble = ['A654','D200','D450','N200','N203','N101','N451']
    # myEnsemble = ['A654','D200','D450','N200','N101','N451']
    
    ### BOOTSTRAP
    ensembleFits = { 'A654': {'Alex': {'Value': 0.7398, 'Error': 0.0084,},
                              'Barb': {'Value': 0.7445, 'Error': 0.0065,},},
                     'N101': {'Alex': {'Value': 0.6633, 'Error': 0.0067,},
                              'Barb': {'Value': 0.663, 'Error': 0.015,},},
                     'D450': {'Alex': {'Value': 0.6152, 'Error': 0.0031,},
                              'Barb': {'Value': 0.6152, 'Error': 0.0016,},},
                     'N451': {'Alex': {'Value': 0.5983, 'Error': 0.0016,},
                              'Barb': {'Value': 0.5997, 'Error': 0.0025,},},
                     'N452': {'Alex': {'Value': 0.5736, 'Error': 0.0020,},
                              'Barb': {'Value': 0.0, 'Error': 0.0,},},
                     'D200': {'Alex': {'Value': 0.5245, 'Error': 0.0016,},
                              'Barb': {'Value': 0.5243, 'Error': 0.0014,},},
                     'N203': {'Alex': {'Value': 0.4880, 'Error': 0.0038,},
                              'Barb': {'Value': 0.4895, 'Error': 0.0051,},},
                     'N200': {'Alex': {'Value': 0.5136, 'Error': 0.0023,},
                              'Barb': {'Value': 0.5138, 'Error': 0.0015,},},}
    ### JACKKNIFE
    # ensembleFits = { 'A654': {'Alex': {'Value': 0.7398, 'Error': 0.0084,},
    #                           'Barb': {'Value': 0.7426, 'Error': 0.0076,},},
    #                  'N101': {'Alex': {'Value': 0.6633, 'Error': 0.0067,},
    #                           'Barb': {'Value': 0.6636, 'Error': 0.0074,},},
    #                  'D450': {'Alex': {'Value': 0.6152, 'Error': 0.0031,},
    #                           'Barb': {'Value': 0.6173, 'Error': 0.0018,},},
    #                  'N451': {'Alex': {'Value': 0.5983, 'Error': 0.0016,},
    #                           'Barb': {'Value': 0.6018, 'Error': 0.0028,},},
    #                  'N452': {'Alex': {'Value': 0.5736, 'Error': 0.0020,},
    #                           'Barb': {'Value': 0.0, 'Error': 0.0,},},
    #                  'D200': {'Alex': {'Value': 0.5245, 'Error': 0.0016,},
    #                           'Barb': {'Value': 0.5248, 'Error': 0.0028,},},
    #                  'N203': {'Alex': {'Value': 0.4880, 'Error': 0.0038,},
    #                           'Barb': {'Value': 0.4906, 'Error': 0.0062,},},
    #                  'N200': {'Alex': {'Value': 0.5136, 'Error': 0.0023,},
    #                           'Barb': {'Value': 0.5155, 'Error': 0.0046,},},}
                    
    
    for ensemble in myEnsemble:
        main_location = f'{ed.outputLocation}/{ensemble}/'
        
        rs_type = 'jk'
        if rs_type=='jk':
            myResamplingScheme='Jackknife'
        elif rs_type=='bt':
            myResamplingScheme='Bootstrap' 
        baseName = f'{main_location}Single_correlators_{rs_type}_{ensemble}_omega_'
        
        if myEffMass:
            barbFile = h5py.File(f'{baseName}barb.h5', 'r')
            alexFile = h5py.File(f'{baseName}alex.h5', 'r')
            
            plot_location = f'{ed.location}/Plots/{ensemble}/SingleHadrons/{myResamplingScheme}/'

            PlotEffectiveMassComparison(barbFile, alexFile, myResamplingScheme, plot_location, 'Barb', 'Alex')
            
            barbFile.close()
            alexFile.close()
            
        if myFitsComp:
            theOmegaName = r'$ am_{\Omega} $'
            print("Not done yet.")
            PlotFitsComparison(ensembleFits, myEnsemble, theOmegaName, f'{ed.location}/Plots/')
        
        
