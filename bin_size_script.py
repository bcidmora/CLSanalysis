import set_of_analysis_functions as vfa
import fitting_script as fs
import correlators_script as cs
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import h5py
import numpy as np
import glob
import os
from pathlib import Path
import set_of_plot_functions as vfp
from datetime import datetime
import re
now = datetime.now()

plt.rcParams["font.family"] = "sans"
plt.rcParams["mathtext.fontset"] = "dejavusans"

def BinSizeAnalysis(data_file, bin_folder, rs_scheme, nr_irreps, the_irreps, t_min, max_bin_size, chosen_fit_range, chosen_bin_size, version, the_weight, **kwargs):

    # get the arguments that need to be passed to FitSingleCorrelators
    only_one_tmin = kwargs.get('one_tmin')
    the_type_fit = kwargs.get('type_fit')
    type_correlated_fit = kwargs.get('type_correlation')
    nr_kbt = kwargs.get('kbt')
    nr_cfgs = kwargs.get('number_cfgs')
    kbt_samples = kwargs.get('own_kbt_list')
    isospin = kwargs.get('isospin_label')
    ensemble = kwargs.get('ensemble')
    plots_only = kwargs.get('plots_only')

    # create directory for correlator files
    correlator_folder = vfa.DIRECTORY_EXISTS(f'{bin_folder}Single_correlators_{rs_scheme}/')
    fits_folder = vfa.DIRECTORY_EXISTS(f'{bin_folder}Fits_{rs_scheme}_v{version}_tmax_{chosen_fit_range[1]}/')
    images_folder = vfa.DIRECTORY_EXISTS(f'{bin_folder}Plots/')


    if not plots_only:
        ###Correlator part
        # create list of time slices for the fit -> fit routine needs a list
        nr = len(the_irreps)
        list_tmax = [chosen_fit_range[1]] * nr
        
        # check if correlators were already calculated -> if folder is empty, calculate them
        if not os.listdir(correlator_folder):
            for bin_i in range(1, max_bin_size + 1):
                cs.SingleCorrelatorAnalysis(data_file, correlator_folder, version, rs_scheme, the_irreps, the_weight, nr_irreps, rebin_on='rb', rb=bin_i, kbt= nr_kbt, number_cfgs=nr_cfgs, own_kbt_list=kbt_samples)

        elif os.listdir(correlator_folder):
            bin_values = [
                int(m.group(1))
                for f in Path(correlator_folder).iterdir()
                if f.is_file() and (m := re.search(r'bin(\d+)', f.name))
            ]

            if not bin_values:
                raise RuntimeError('No bin files found!')

            max_bin = max(bin_values)
            missing = list(range(max_bin + 1, max_bin_size + 1))
            if missing:
                for bin_i in missing:
                    cs.SingleCorrelatorAnalysis(data_file, correlator_folder, version, rs_scheme, the_irreps, the_weight, rebin_on='rb', rb=bin_i, kbt=nr_kbt, number_cfgs=nr_cfgs, own_kbt_list=kbt_samples)
        bin_size= int()

        ### Fits part
        # perform all the fits for the correlator files for the chosen maximal fit range if the version does not already exist
        if not os.listdir(fits_folder):

            for file in Path(correlator_folder).iterdir():

                # get bin size from file name
                parts = file.name.split('bin')
                if len(parts) > 1:
                    number = ''
                    for char in parts[1]:
                        if char.isdigit():
                            number += char
                        else:
                            break
                    if number:
                        bin_size = int(number)

                fit_correlator = h5py.File(
                    fits_folder + f'Single_correlators_{rs_scheme}_bin{bin_size}_fits.h5', 'a')

                fit_file = h5py.File(file, 'r')
                fs.FitSingleCorrelators(fit_file, fit_correlator, rs_scheme, nr_irreps, t_min, list_tmax, one_tmin=only_one_tmin,
                                    type_fit=the_type_fit, type_correlation=type_correlated_fit)
                fit_file.close()
                fit_correlator.close()

        elif os.listdir(fits_folder):

            bin_values = [
                int(m.group(1))
                for f in Path(fits_folder).iterdir()
                if f.is_file() and (m := re.search(r'bin(\d+)', f.name))
            ]

            if not bin_values:
                raise RuntimeError('No bin files found!')

            max_bin = max(bin_values)
            missing = list(range(max_bin + 1, max_bin_size + 1))
            if missing:
                for bin_i in missing:
                    fit_correlator = h5py.File(
                        fits_folder +  f'Single_correlators_{rs_scheme}_bin{bin_i}_fits.h5', 'a')

                    fit_file = h5py.File(correlator_folder + f'Single_correlators_{rs_scheme}_bin{bin_i}_v{version}.h5', 'r')
                    fs.FitSingleCorrelators(fit_file, fit_correlator, rs_scheme, nr_irreps, t_min, list_tmax,
                                            one_tmin=only_one_tmin,
                                            type_fit=the_type_fit, type_correlation=type_correlated_fit)
                    fit_file.close()
                    fit_correlator.close()

    # test, if the bin sizes are consecutive in the fit file, start at one and only occur once before we plot any nonsense
    bin_values = [
        int(m.group(1))
        for f in Path(fits_folder).iterdir()
        if f.is_file() and (m := re.search(r'bin(\d+)', f.name))
    ]

    bin_values_sorted = sorted(bin_values)
    print(bin_values_sorted)

    # check if fits_folder is empty
    if not bin_values_sorted:
        raise ValueError('No fit files found in fits folder!')

    # check for duplicates
    if len(bin_values_sorted) != len(set(bin_values_sorted)):
        raise ValueError('Duplicate bin sizes detected in fits folder!')

    # check starting at 1
    if bin_values_sorted[0] != 1:
        raise ValueError(f'Bin numbering does not start at 1 (starts at {bin_values_sorted[0]})')

    # check consecutiveness
    expected = list(range(1, bin_values_sorted[-1] + 1))

    if bin_values_sorted != expected:
        missing = sorted(set(expected) - set(bin_values_sorted))
        raise ValueError(f'Missing bin sizes detected: {missing}')

    ### Plot part
    # find the file with bin size 1 as norm
    t_min_shift = int()
    files = [f for f in Path(fits_folder).glob('*bin1_*') if f.is_file()]
    if len(files) == 1:
        fit_norm = files[0]
        t_min_shift = int(list(h5py.File(fit_norm, 'r').get(f'{the_irreps[nr_irreps[0]]}/1exp/Tmin/Correlated/Mean'))[0][0])

    else: raise ValueError('There are too many files with bin=1 in the Single Correlators folder. Every bin size should occur only once!')

    # plot the actual analysis results

    sigmas_list, chis_list, n_bin_list = [], [], []
    t_min_index = chosen_fit_range[0] - t_min_shift
    if t_min_index < 0:
        raise ValueError('Chosen t_min in myBinSizeFitRange is smaller than range that is covered by the FitSingleCorrelators routine')

    sigmas_1 = list(h5py.File(fit_norm, 'r').get(f'{the_irreps[nr_irreps[0]]}/1exp/Tmin/Correlated/Mean'))[3][t_min_index] ** 2
    chi_1 = list(h5py.File(fit_norm, 'r').get(f'{the_irreps[nr_irreps[0]]}/1exp/Tmin/Correlated/Mean'))[4][t_min_index]



    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 5))
    ax1, ax2 = axes.flatten()
    for file in Path(fits_folder).glob('*'):
        if file.is_file():
            # extract bin number
            match = re.search(r'bin(\d+)', file.name)
            if match:
                bin_size = int(match.group(1))
                # only proceed if bin_size <= max_bin_size
                if bin_size <= max_bin_size:
                    data = h5py.File(file, 'r')
                    sigma = list(data.get(f'{the_irreps[nr_irreps[0]]}/1exp/Tmin/Correlated/Mean'))[3][t_min_index] ** 2
                    chi = list(data.get(f'{the_irreps[nr_irreps[0]]}/1exp/Tmin/Correlated/Mean'))[4][t_min_index]
                    n_bin = int(vfp.GET_RESAMPLING_BINNING(data)[2])
                    sigmas_list.append(sigma / sigmas_1)
                    chis_list.append(chi / chi_1)
                    n_bin_list.append(n_bin)
    n_bin_list = np.array(n_bin_list)
    indices = np.argsort(n_bin_list)
    n_bin_sorted = n_bin_list[indices]
    chis_sorted = np.array(chis_list)[indices]
    sigmas_sorted = np.array(sigmas_list)[indices]
    ax1.scatter(n_bin_sorted, sigmas_sorted, color= vfp.colors[0], s=80)
    ax2.scatter(n_bin_sorted, chis_sorted, color= vfp.colors[0], s=80)
    ax1.scatter(n_bin_sorted[chosen_bin_size-1], sigmas_sorted[chosen_bin_size-1], color= vfp.colors[5],
                label=r'chosen $N_{\mathrm{bin}}$', s=90, marker='s')
    ax2.scatter(n_bin_sorted[chosen_bin_size-1], chis_sorted[chosen_bin_size-1], color= vfp.colors[5],
                label=r'chosen $N_{\mathrm{bin}}$', s=90, marker='s')
    ax1.tick_params(axis='both', which='major',
                    length=10, direction='inout', top=True, right=True, labelsize=20)
    ax2.tick_params(axis='both', which='major',
                    length=10, direction='inout', top=True, right=True, labelsize=20)

    if max_bin_size > 10:
        ax1.xaxis.set_major_locator(MultipleLocator(2))
        ax2.xaxis.set_major_locator(MultipleLocator(2))

    else:
        ax1.xaxis.set_major_locator(MultipleLocator(1))
        ax2.xaxis.set_major_locator(MultipleLocator(1))


    ax1.set_xlabel(r'$N_{\mathrm{bin}}$', fontsize=24)
    ax2.set_xlabel(r'$N_{\mathrm{bin}}$', fontsize=24)
    ax1.tick_params(axis='both', labelsize=20)
    ax2.tick_params(axis='both', labelsize=20)
    ax1.set_ylabel(r'$\sigma^2_{N_{\mathrm{bin}}}$ / $\sigma^2_{N_1}$', fontsize=24)
    ax2.set_ylabel(r'$\chi ^2_{N_{\mathrm{bin}}}$ / $\chi^2_{N_1}$ ', fontsize=24)
    fig.suptitle(rf'$I$ = {isospin} {vfp.GET_IRREP_LOGO(the_irreps[nr_irreps[0]])} {vfp.IRREP_TO_INDEX(the_irreps[nr_irreps[0]][5:])} for $t_{{\mathrm{{min}}}}$ = {chosen_fit_range[0]}, $t_{{\mathrm{{max}}}}$ = {chosen_fit_range[1]}, {ensemble}', fontsize=24)
    ax1.legend(fontsize=20)
    ax2.legend(fontsize=20)
    fig.tight_layout()
    fig.savefig(images_folder + f'bin-size_analysis_v{version}_{now.strftime("%d_%m_%H%M")}.pdf', format='pdf')

    print('.............................................................................')
    print('                   Bin-Size Analysis Completed                               ')
    print('.............................................................................')
