import numpy as np
import re
import os
import h5py
from collections import defaultdict
from datetime import datetime
import set_of_analysis_functions as vfa
import set_of_plot_functions as vfp
import set_of_layout_functions as vfl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
# from adjustText import adjust_text
from itertools import cycle

def CORRELATED_DIFFERENCES(non_int_list, t_mins_range_s, t_mins_range_m, singles_files, multi_file, save_file_dir, t_mins_shift, the_lat_size, the_nr_irreps, the_nr_op, the_T0, the_td, the_rs, the_plot, the_threshold_masses, the_ensemble, the_iso_label, two_body, three_body, fmToMev, the_bin):

    print("\n                     CORRELATED DIFFERENCES COMPUTATION + PLOTTING\n")

    irreps = list(multi_file.keys())
    now = datetime.now()
    version = now.strftime("%d_%m_%H%M")

    save_file = h5py.File(save_file_dir + f"Correlated_differences_{the_rs}{the_bin}_fits_v{version}.h5", "w")

    try:
        test = int(the_td)
        suffix = '_pv'
    except ValueError:
        suffix = '_run'
    except TypeError:
        suffix = ''

    norm = (2 * np.pi / the_lat_size) ** 2

    if the_plot:

        thresholds = {}

        # --- Two-body thresholds
        for particles in two_body:
            if all(p in the_threshold_masses for p in particles):
                ch_name = "+".join(particles)
                #label_name = "+".join(labels)

                thresholds[ch_name] = {
                    "value": sum(the_threshold_masses[p] for p in particles),
                    #"label": label_name
                }

        # --- Three-body thresholds
        for particles in three_body:
            if all(p in the_threshold_masses for p in particles):
                ch_name = "+".join(particles)
                #label_name = "+".join(labels)

                thresholds[ch_name] = {
                    "value": sum(the_threshold_masses[p] for p in particles),
                    #"label": label_name
                }

        labels_x = [vfp.GET_IRREP_LOGO(irreps[j]) + ' ' + vfp.IRREP_TO_INDEX(irreps[j].split('_')[1]) for j in
                    range(len(irreps))]
        ### Plot Overview
        fig_overview, ax_overview = plt.subplots(figsize=(14, 8))
        ax_overview.set_xticks(range(len(labels_x)))
        ax_overview.set_xticklabels(labels_x, rotation=90, va='bottom', ha='center')
        ax_overview.set_xlim(-0.8, len(irreps) - 1 + 0.6)
        ax_overview.set_ylabel(r'$a E_{\mathrm{com}}$', fontsize=20)
        ax_overview.tick_params(axis='both', which='major',
                       length=10, direction='inout', top=True, right=True, labelsize=20)
        offset_o = mtransforms.ScaledTranslation(0, -1.6, fig_overview.dpi_scale_trans)

        ### Plot Truncation
        fig_truncated, ax_truncated = plt.subplots(figsize=(14, 8))
        ax_truncated.set_xticks(range(len(labels_x)))
        ax_truncated.set_xticklabels(labels_x, rotation=90, va='bottom', ha='center')
        ax_truncated.set_xlim(-0.8, len(irreps) - 1 + 0.6)
        ax_truncated.set_ylabel(r'$E_{\mathrm{com}}$ [MeV]', fontsize=20)
        ax_truncated.tick_params(axis='both', which='major',
                                length=10, direction='inout', top=True, right=True, labelsize=20)
        offset_t = mtransforms.ScaledTranslation(0, -1.6, fig_truncated.dpi_scale_trans)

        for (ch, data), ls in zip(thresholds.items(), cycle(vfp.line_style)):
            ax_overview.axhline(y=data["value"], color=vfp.greys[0], linestyle=ls, label=vfp.CONVERT_LABEL_TO_LATEX(ch, vfp.particel_dictionary), alpha=0.3)

        for label in ax_overview.get_xticklabels():
            label.set_transform(label.get_transform() + offset_o)

        for label in ax_truncated.get_xticklabels():
            label.set_transform(label.get_transform() + offset_t)

        min_name, min_data = min(thresholds.items(), key=lambda x: x[1]["value"])

        min_th = min_data["value"]
        label_t = vfp.CONVERT_LABEL_TO_LATEX(min_name, vfp.particel_dictionary)

        ax_truncated.axhline(
            vfa.UNIT_CONVERSION(min_th, the_lat_size, fmToMev), color=vfp.greys[0], linestyle='--', label=label_t, alpha=0.3)

    # getting the reference energy levels for the plot and the energies of the single mesons
    non_int_energies_reference, non_int_energies_reference_rs, non_int_energies_mean, non_int_energies_rs = vfa.NON_INTERACTING_LEVELS(
        non_int_list, t_mins_range_s, singles_files, t_mins_shift, norm)

    for index, j in enumerate(the_nr_irreps):

        irrep  = irreps[j]
        frame_momentum = int(irrep[3])
        if irrep in save_file: del save_file[irrep]
        new_irr_group = save_file.create_group(irrep)
        op = list(multi_file[irrep].keys())

        for i in the_nr_op[index]:
            if "/".join([irrep, op[i]]) in save_file: del save_file[irrep][op[i]]
            new_op_group = new_irr_group.create_group(op[i])
            n_basis = len(list(multi_file[irrep][op[i]].get('Operators')))

            mean_diff = []
            mean_sigmas = []
            mean_for_plot = []
            upper_limit = []
            lower_limit = []

            for n in range(n_basis):

                multi_hads_mean = list(multi_file[irrep][op[i]].get(f'1exp/t0_{the_T0}{suffix}/Tmin/Correlated/Mean/lambda_{n}'))[2][t_mins_range_m[j][i][n]- t_mins_shift]
                multi_hads_rs = multi_file[irrep][op[i]].get(f'1exp/t0_{the_T0}{suffix}/Tmin/Correlated/Resampled/lambda_{n}')[()][t_mins_range_m[j][i][n]- t_mins_shift]

                index_closest = np.argmin(np.abs(non_int_energies_mean[j] - multi_hads_mean))
                mean_diff_one_ev = multi_hads_mean - non_int_energies_mean[j][index_closest]
                mean_diff.append(mean_diff_one_ev)

                rs_diff_one_ev = multi_hads_rs - np.array(non_int_energies_rs[j][index_closest], dtype=float)

                # To plot the differences, the reusults from the dispersion relation are added up
                mean_diff_plus_disp_one_ev = mean_diff_one_ev + non_int_energies_reference[j][index_closest]


                mean_sigma_one_ev = vfa.STD_DEV(rs_diff_one_ev, mean_diff_one_ev, "bt")
                mean_sigmas.append(mean_sigma_one_ev)

                mean_com_one_ev = vfa.CENTER_OF_MASS(frame_momentum,
                                            mean_diff_plus_disp_one_ev, norm)
                mean_for_plot.append(mean_com_one_ev)

                # Take differences for plot later, add non-interacting levels and boost to com frame
                upper_limit_one_ev = vfa.CENTER_OF_MASS(frame_momentum, mean_diff_plus_disp_one_ev + mean_sigma_one_ev, norm)
                lower_limit_one_ev = vfa.CENTER_OF_MASS(frame_momentum, mean_diff_plus_disp_one_ev - mean_sigma_one_ev, norm)

                upper_limit.append(upper_limit_one_ev)
                lower_limit.append(lower_limit_one_ev)

            cd = new_op_group.create_group('Correlated_Differences')
            data = cd.create_group('Data')
            plot = cd.create_group('Plot_Data_Plus_Non_Interacting_COM')
            data.create_dataset('Mean_Difference', data=mean_diff)
            data.create_dataset('Mean_Sigma', data=mean_sigmas)
            plot.create_dataset('Mean', data=mean_for_plot)
            plot.create_dataset('Upper_Limit', data=upper_limit)
            plot.create_dataset('Lower_Limit', data=lower_limit)


            if the_plot:

                ### Plot Overview
                text_o = []
                for k in range(len(non_int_energies_reference[j])):
                    text_o.append(ax_overview.text(j - 0.35, vfa.CENTER_OF_MASS(frame_momentum,non_int_energies_reference[j][k], norm),
                                         vfp.CONVERT_LABEL_TO_LATEX(non_int_list[j][k], vfp.particel_dictionary), va='center', ha='right'))
                    label = r'Non-Inerating Dispersion $O_{(p^2)}$' if k == 0 else None
                    ax_overview.plot([j - 0.2, j + 0.2], [vfa.CENTER_OF_MASS(frame_momentum, non_int_energies_reference[j][k],norm), vfa.CENTER_OF_MASS( frame_momentum, non_int_energies_reference[j][k], norm)], linestyle='--',
                            color=vfp.all_colors[2], label=label)
                for l in range(len(mean_for_plot)):
                    ax_overview.fill_between([j - 0.2, j + 0.2], y1=lower_limit[l], y2=upper_limit[l], color=vfp.greys[1],
                                    alpha=0.3)
                    label = 'Measured Mean' if l == 0 else None
                    ax_overview.plot([j - 0.2, j + 0.2], [mean_for_plot[l], mean_for_plot[l]], color=vfp.all_colors[0], label=label)
                if j == 0: ax_overview.legend(fontsize=16, loc='upper right', ncols=2)
                adjust_text(text_o)

                ### Plot Truncation
                text_t = []
                first_non_int = True
                for k in range(len(non_int_energies_reference[j])):
                    if all(non_int_energies_reference[j][k] < th["value"] for th in thresholds.values()):
                        text_t.append(ax_truncated.text(j - 0.35, vfa.UNIT_CONVERSION(vfa.CENTER_OF_MASS(frame_momentum, non_int_energies_reference[j][k], norm), the_lat_size, fmToMev),
                                            vfp.CONVERT_LABEL_TO_LATEX(non_int_list[j][k], vfp.particel_dictionary),
                                            va='center',
                                            ha='right'))
                        label = 'No Interaction' if first_non_int else None
                        ax_truncated.plot([j - 0.2, j + 0.2], [vfa.UNIT_CONVERSION(vfa.CENTER_OF_MASS(frame_momentum, non_int_energies_reference[j][k], norm), the_lat_size, fmToMev),
                                                     vfa.UNIT_CONVERSION(vfa.CENTER_OF_MASS(frame_momentum, non_int_energies_reference[j][k],  norm), the_lat_size, fmToMev)],
                                linestyle='--',
                                color=vfp.all_colors[2],
                                label=label)
                        first_non_int = False
                    else:
                        continue
                first_int = True
                for l in range(len(mean_for_plot)):
                    if all(mean_for_plot[l] < th["value"] for th in thresholds.values()):
                        ax_truncated.fill_between([j - 0.2, j + 0.2],
                                        y1=vfa.UNIT_CONVERSION(lower_limit[l], the_lat_size, fmToMev),
                                        y2=vfa.UNIT_CONVERSION(upper_limit[l], the_lat_size, fmToMev),
                                        color=vfp.greys[1], alpha=0.3)
                        label_i = 'Measured Mean' if first_int else None
                        ax_truncated.plot([j - 0.2, j + 0.2], [vfa.UNIT_CONVERSION(mean_for_plot[l], the_lat_size, fmToMev),
                                                     vfa.UNIT_CONVERSION(mean_for_plot[l], the_lat_size, fmToMev)],
                                color=vfp.all_colors[0],
                                label=label_i)
                        first_int = False
                adjust_text(text_t)
                if j == 0: ax_truncated.legend(fontsize=16, loc='upper right')

    if the_plot:
        fig_overview.suptitle(f'Overview Fit Results {the_ensemble}, $I$ = {the_iso_label}, {vfp.GET_RESAMPLING_BINNING(multi_file)[1]}, {vfp.GET_RESAMPLING_BINNING(multi_file)[0]}', fontsize=24)
        fig_overview.tight_layout()
        fig_overview.savefig(save_file_dir + f'v_{version}_overview.pdf', format='pdf')

        fig_truncated.suptitle(f'Overview Truncated Fit Results {the_ensemble}, $I$ = {the_iso_label}, {vfp.GET_RESAMPLING_BINNING(multi_file)[1]}, {vfp.GET_RESAMPLING_BINNING(multi_file)[0]}', fontsize=24)
        fig_truncated.tight_layout()
        fig_truncated.savefig(save_file_dir + f'v_{version}_truncated.pdf', format='pdf')
