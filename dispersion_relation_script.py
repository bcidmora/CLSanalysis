import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from collections import defaultdict
import h5py
import statistics
import numpy as np
import math
import os
from datetime import datetime
import glob

import set_of_plot_functions as vfp
import matplotlib.gridspec as gridspec
import os
import matplotlib.ticker as ticker

import set_of_analysis_functions as vfa
import set_of_plot_functions as vfp
import set_of_layout_functions as vfl

import matplotlib
matplotlib.rcParams.update({
    "mathtext.fontset": "cm",          
    "font.family": "serif",
    "font.serif": ["CMU Serif"],       
    "axes.unicode_minus": False,       
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def DISPERSION_RELATION(fit_file, the_rs, the_tmin_choice, the_lat_size, the_disp_location, the_plot_type, the_ensemble):

    print("\n                     DISPERSION RELATION \n")

    now = datetime.now()
    the_version = now.strftime("%d_%m_%H%M")
    
    the_irreps = list(fit_file.keys())
    tmin_per_irrep = dict(zip(the_irreps, the_tmin_choice))
    
    the_norm = (2 * np.pi / the_lat_size) ** 2
    
    the_irrep_dict = defaultdict(dict)
    for the_irrep in the_irreps:
        the_info = vfp.IrrepInfo(the_irrep)
        the_hadron = the_info.Hadron
        the_hadron_isolabel = the_info.HadronIsospin
        the_momentum = the_info.Momentum
        
        the_irrep_dict[the_hadron][the_momentum] = the_irrep
        
    the_selected_hadrons = vfl.HADRON_MENU(the_irreps)
    
    for this_hadron in the_selected_hadrons:
        print(f"Processing hadron: {this_hadron}")
    
        the_momentum_map = the_irrep_dict[this_hadron]
        the_momenta = sorted(the_momentum_map.keys())
        
        the_energies, the_energies_sampled = [], []

        x_ticks_dict = {i: m for i, m in enumerate(the_momenta)}
        x_pos = list(x_ticks_dict.keys())
        
        for m in the_momenta:
            match = the_momentum_map[m]
            the_tmin = tmin_per_irrep[match]

            t_min_ll = int(list(fit_file[f'{match}/1exp/Tmin/Correlated/Mean'])[0][0])

            energy_i = np.asarray(fit_file[f'{match}/1exp/Tmin/Correlated/Mean'])[2][the_tmin - t_min_ll]
            energy_i_sample = np.asarray(fit_file[f'{match}/1exp/Tmin/Correlated/Resampled'])[the_tmin - t_min_ll]

            the_energies.append(energy_i)
            the_energies_sampled.append(energy_i_sample)

        the_energies = np.asarray(the_energies)
        the_energies_sampled = np.asarray(the_energies_sampled)
        the_momenta = np.asarray(the_momenta)

        if the_plot_type == 'abs' or the_plot_type == 'both':
            energy_dispersion, error_dispersion, energy_calculated, energies_sampled_normalized = [], [], [], []
            sigmas_dispersion, sigmas_calculated = [], []

            for k in range(the_momenta.size):
                energy_dispersion.append(vfa.CONTINUUM_DISP_REL(the_momenta[k], the_energies[0], the_norm))
                energy_calculated.append(the_energies[k] / the_energies[0])

                useful = []
                useful_2 = []
                for kk in range(the_energies_sampled[k].size):
                    useful.append(the_energies_sampled[k][kk] / the_energies_sampled[0][kk])
                    useful_2.append(vfa.CONTINUUM_DISP_REL(the_momenta[k], the_energies_sampled[0][kk], the_norm))
                energies_sampled_normalized.append(useful)
                error_dispersion.append(useful_2)

                sigmas_dispersion.append(vfa.STD_DEV(error_dispersion[k], energy_dispersion[k], the_rs))
                sigmas_calculated.append(vfa.STD_DEV(energies_sampled_normalized[k], energy_calculated[k], the_rs))

            data_min, data_max = [], []
            range_size, center = [], []

            for i in range(the_momenta.size):
                mini = min(energy_calculated[i] + sigmas_calculated[i], energy_calculated[i] - sigmas_calculated[i], energy_dispersion[i] + sigmas_dispersion[i], energy_dispersion[i] - sigmas_dispersion[i])
                maxi = max(energy_calculated[i] + sigmas_calculated[i], energy_calculated[i] - sigmas_calculated[i], energy_dispersion[i] + sigmas_dispersion[i], energy_dispersion[i] - sigmas_dispersion[i])
                center_i = np.mean([mini, maxi])
                
                center.append(center_i)
                data_min.append(mini)
                data_max.append(maxi)
            
            data_min = np.asarray(data_min)
            data_max = np.asarray(data_max)
            range_size = max(list(data_max - data_min))

            if len(the_momenta) < 6:
                fig = plt.figure(figsize=(10, 7))
            else: fig = plt.figure(figsize=(21, 7))
            cols = len(the_momenta)
            
            gs = gridspec.GridSpec(2, cols, height_ratios=[1, 1])  
            ax_top = fig.add_subplot(gs[0, :])
            axes_bottom = [fig.add_subplot(gs[1, i]) for i in range(cols)]  
            
            for i, ax in enumerate(axes_bottom):
                ax.errorbar(x_pos[i] , energy_dispersion[i], yerr= sigmas_dispersion[i], fmt='^', markersize=8, capsize=4, color=vfp.the_colors[1], ecolor=vfp.greys[1])
                ax.errorbar(x_pos[i], energy_calculated[i], yerr=sigmas_calculated[i], fmt='o', markersize=8, capsize=4, color=vfp.the_colors[0], ecolor=vfp.greys[0])

                ax.set_ylim(center[i] - range_size / 2 - 0.005, center[i] + range_size / 2 + 0.005)

                ax.tick_params(axis='both', which='major', length=10, direction='inout', labelsize=16, top=True, right=True)
                ax.set_xticks([x_pos[i]])
                ax.set_xticklabels([the_momenta[i]])

            ax_top.errorbar(x_pos, energy_dispersion, yerr=sigmas_dispersion, fmt='^', markersize=8, capsize=4, color=vfp.the_colors[1], ecolor=vfp.greys[1], label=r'$D(am, d^{2})$')
            ax_top.errorbar(x_pos, energy_calculated, yerr=sigmas_calculated, fmt='o', markersize=8, capsize=4, color=vfp.the_colors[0], ecolor=vfp.greys[0], label=f'Fit results')

            ax_top.set_xlabel(r'Momentum $(\mathrm{d}^{2})$', fontsize=20)

            ax_top.set_ylabel(r'$aE_{lab}/am$', fontsize=20)
            axes_bottom[0].set_ylabel(r'$aE_{lab}/am$', fontsize=20)

            ax_top.set_title(f'Dispersion Relation', fontsize=20)
            ax_top.legend(fontsize=14)
            ax_top.tick_params(axis='both', which='major', length=10, direction='inout', labelsize=16, top=True, right=True)
            ax_top.set_xticks(x_pos, [x_ticks_dict[i] for i in x_pos])
            fig.suptitle(rf'$I =$ {vfp.OPERATORS_SH_ISOSPIN(this_hadron)}, ' + vfp.GET_RESAMPLING_BINNING(fit_file)[0] + ', ' + vfp.GET_RESAMPLING_BINNING(fit_file)[1] + ', ' + the_ensemble, fontsize=16)

            fig.tight_layout()
            fig.savefig(f'{the_disp_location}{this_hadron}_{the_version}_absolute_dispersion_relation.pdf', format='pdf')


        if the_plot_type == 'rel' or the_plot_type == 'both':

            energies_relative, energies_sampled_relative, sigmas_relative = [], [], []

            for i in range(the_momenta.size):
                energies_relative.append(vfa.RELATIVE_DISP_REL(the_momenta[i], the_energies[0], the_energies[i], the_norm))
                useful = []
                for k in range(the_energies_sampled[i].size):
                    useful.append(vfa.RELATIVE_DISP_REL(the_momenta[i], the_energies_sampled[0][k],the_energies_sampled[i][k], the_norm))
                energies_sampled_relative.append(useful)
                sigmas_relative.append(vfa.STD_DEV(energies_sampled_relative[i], energies_relative[i], the_rs))

            fig, ax = plt.subplots(figsize=(8, 4))
            ax.errorbar(x_pos, energies_relative, yerr=sigmas_relative, fmt='o', markersize=6, capsize=3, color=vfp.colors[0], ecolor=vfp.greys[0])
            ax.tick_params(axis='both', which='major', length=10, direction='inout', labelsize=12, top=True, right=True)
            ax.set_xlabel(r'Momentum $d^{2}$', fontsize=16)
            ax.axhline(1, 0, 10, linestyle=':', color=vfp.greys[2])
            ax.set_ylabel(r'$aE/ D(am, d^{2})$', fontsize=16)
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.set_xticks(x_pos, [x_ticks_dict[i] for i in x_pos])
            ax.set_title(f'Dispersion Relation', fontsize=20)
            fig.suptitle(rf'$I =$ {vfp.OPERATORS_SH_ISOSPIN(this_hadron)}, ' + vfp.GET_RESAMPLING_BINNING(fit_file)[0] + ', ' + vfp.GET_RESAMPLING_BINNING(fit_file)[1] + ', ' + the_ensemble, fontsize=16)
            fig.tight_layout()
            fig.savefig(f'{the_disp_location}{this_hadron}_v{the_version}_relative_dispersion_relation.pdf', format='pdf')
