import numpy as np
import h5py
import time
import sys
import os
import set_of_analysis_functions as vfa
import set_of_layout_functions as vfl
from iminuit import Minuit
import warnings
warnings.filterwarnings('ignore')


def SingleCorrelatorEffectiveMass(the_single_correlator_data, the_type_rs,**kwargs):   
    
    print("                  EFFECTIVE MASSES COMPUTATION IB\n")
    
    ### The irreps
    the_list_name_irreps = list(the_single_correlator_data.keys())
    
    begin_time = time.time()     
    ### Extracting data from file
    this_data = the_single_correlator_data[the_list_name_irreps[0]]
    
    ### Mean values of the real part of the correlator to get the effective masses
    the_mean_corr_real = np.asarray(this_data['Correlators/Real/Mean'])
    
    ### The real part of the resampled data
    the_rs_real = np.transpose(np.asarray(this_data['Correlators/Real/Resampled']), (1,0))
    
    ### Effective Mass computation (central values)
    the_em_rs_f =  vfa.EFF_MASS(the_mean_corr_real, 1)
    
    ### Loop over the resampled data            
    the_em_rs = np.transpose(np.array([vfa.EFF_MASS(the_rs_real.real[l],1) for l in range(len(the_rs_real))]), (1,0))
    
    ### Mean resampled data
    the_mrs_f_real_rs = np.mean(the_em_rs, axis=1)
    
    ### Sigma values for the resamples
    the_sigma_eff_mass = vfa.STD_DEV_MEAN(the_em_rs, the_mrs_f_real_rs, the_type_rs)
    
    ### CORRECTIONS
    this_data_mass = the_single_correlator_data[the_list_name_irreps[1]]
    this_data_qed = the_single_correlator_data[the_list_name_irreps[2]]
    
    ### Central values 
    the_mean_corr_real_mass = np.asarray(this_data_mass['Correlators/Real/Mean'])
    the_mean_corr_real_qed = np.asarray(this_data_qed['Correlators/Real/Mean'])
            
    ### Resamples
    the_rs_real_mass = np.transpose(np.asarray(this_data_mass['Correlators/Real/Resampled']), (1,0))
    the_rs_real_qed = np.transpose(np.asarray(this_data_qed['Correlators/Real/Resampled']), (1,0))
    
    ### Effective Mass computation Corrections
    the_em_rs_f_mass =  vfa.EFF_MASS_CORRECTIONS(the_mean_corr_real_mass, the_mean_corr_real)
    the_em_rs_f_qed =  vfa.EFF_MASS_CORRECTIONS(the_mean_corr_real_qed, the_mean_corr_real)

    ### Loop over the resampled data            
    the_em_rs_mass = np.transpose(np.asarray([vfa.EFF_MASS_CORRECTIONS(the_rs_real_mass.real[l],the_rs_real[l]) for l in range(len(the_rs_real))]), (1,0))
    
    the_em_rs_qed = np.transpose(np.asarray([vfa.EFF_MASS_CORRECTIONS(the_rs_real_qed.real[l],the_rs_real[l]) for l in range(len(the_rs_real))]), (1,0))
    
    ### Mean values resamples corrections
    the_mrs_f_real_rs_mass = np.mean(the_em_rs_mass, axis=1)
    the_mrs_f_real_rs_qed = np.mean(the_em_rs_qed, axis=1)
    
    ### Standard deviation
    the_sigma_eff_mass_ms = vfa.STD_DEV_MEAN(the_em_rs_mass, the_mrs_f_real_rs_mass, the_type_rs)
    the_sigma_eff_mass_qed = vfa.STD_DEV_MEAN(the_em_rs_qed, the_mrs_f_real_rs_qed, the_type_rs)
    
    ### Ratio of central values
    the_ratio_mean_mass = the_mean_corr_real_mass / the_mean_corr_real
    the_ratio_mean_qed = the_mean_corr_real_qed / the_mean_corr_real
    
    the_single_correlator_data[f'{the_list_name_irreps[1]}/Correlators/Real/Mean'][...] = the_ratio_mean_mass
    the_single_correlator_data[f'{the_list_name_irreps[2]}/Correlators/Real/Mean'][...] = the_ratio_mean_qed
    
    ### Ratios of resamples
    the_ratio_rs_mass = np.transpose(the_rs_real_mass / the_rs_real, (1,0))
    the_ratio_rs_qed = np.transpose(the_rs_real_qed / the_rs_real, (1,0))
    
    the_single_correlator_data[f'{the_list_name_irreps[1]}/Correlators/Real/Resampled'][...] = the_ratio_rs_mass
    the_single_correlator_data[f'{the_list_name_irreps[2]}/Correlators/Real/Resampled'][...] = the_ratio_rs_qed
    
    ### New covariance matrix
    the_ratio_cov_matrix_mass = np.asarray(vfa.COV_MATRIX(the_ratio_rs_mass, vfa.MEAN(the_ratio_rs_mass), the_type_rs))
    the_ratio_cov_matrix_qed = np.asarray(vfa.COV_MATRIX(the_ratio_rs_qed, vfa.MEAN(the_ratio_rs_qed), the_type_rs))
    
    the_single_correlator_data[f'{the_list_name_irreps[1]}/Correlators/Real/Covariance_matrix'][...] = the_ratio_cov_matrix_mass
    the_single_correlator_data[f'{the_list_name_irreps[2]}/Correlators/Real/Covariance_matrix'][...] = the_ratio_cov_matrix_qed
    
    the_effs_masses = [the_em_rs_f, the_em_rs_f_mass, the_em_rs_f_qed]
    the_sigmas = [the_sigma_eff_mass, the_sigma_eff_mass_ms, the_sigma_eff_mass_qed]
    
    for j in range(len(the_list_name_irreps)):
        this_data = the_single_correlator_data[the_list_name_irreps[j]]
        
        ### If the Effective Masses were already computed, this part gets deleted and created a new branch.
        if 'Effective_masses' in this_data.keys(): 
            del the_single_correlator_data[f'{the_list_name_irreps[j]}/Effective_masses']
        
        group_em = this_data.create_group('Effective_masses')
        
        group_em.create_dataset('Mean', data = the_effs_masses[j])
        group_em.create_dataset('Sigmas', data = the_sigmas[j])
        
        print(f'Irrep nr.: {j+1} out of {len(the_list_name_irreps)}')
    end_time = time.time()
    print(f'TOTAL TIME TAKEN: {(end_time-begin_time)/60} mins')




def FitSingleCorrelators(the_data, the_fit_data, the_type_rs, the_list_tmaxs, the_irreps, the_fit_key, type_correlated_fit, **kwargs): 
    
    the_fit_map_ib = {'1': {"iso": {"n_params": 2, "model": vfa.SINGLE_EXPONENTIAL, "params": ('a0', 'e0'),},
                            "ib": {"n_params": 2, "model": vfa.SINGLE_EXP_CORRECTIONS_IB, "params": ('a0', 'e0')},
                            "label": "Single Exponential Fit",},
                      '2': {"iso": {"n_params": 4, "model": vfa.DOUBLE_EXPONENTIAL, "params": ('a0', 'e0', 'b', 'dm'),},
                            "ib": {"n_params": 6, "model": vfa.DOUBLE_EXP_CORRECTIONS_IB, "params": ('b0', 'a0', 'b', 'dm', 'e0', 'e1')},
                            "label": "Double Exponential Fit"},}    
                      
    ### Which correlator to fit: isoQCD or ib corrections
    the_corr_choice = kwargs.get('iso_or_ib')
    
    ### The type of fit
    the_type_fit = the_fit_map_ib[the_fit_key]
    
    print(f"                     FITTING: IB {the_type_fit.get("label", "")}\n")
    
    ### The name of the irreps
    the_s_irreps = list(the_data.keys())
    
    ### Storing all the correlators before calculating anything
    the_corr_temp = {key: np.asarray(the_data[f'{key}/Correlators/Real/Mean'], dtype=np.float64) for key in the_s_irreps}
    the_corr_rs_temp = {key: vfa.NT_TO_NCFGS(np.asarray(the_data[f'{key}/Correlators/Real/Resampled'], dtype=np.float64)) for key in the_s_irreps}
    
    ### The isoQCD irrep, because we need these correlators for all of the fits
    the_iso_irrep = the_s_irreps[0]
    
    ### The isoQCD correlator
    the_iso_corr = the_corr_temp[the_iso_irrep]
    
    ### The isoQCD resamples
    the_iso_corr_rs = the_corr_rs_temp[the_iso_irrep]
    
    ### Now separting the fit forms into the isoQCD and the IB corrections
    if the_corr_choice=='isoQCD' or the_corr_choice=='both':
        
        ### isoQCD part
        the_iso = the_type_fit['iso'] 
        
        ### Choosing the configuration for the isoQCD part
        the_n_params_iso, da_minimization_iso, the_fit_params_iso = the_iso["n_params"], the_iso["model"], the_iso["params"]
        
        ### List of operators of this irrep
        the_op_list = list(the_data[f'{the_iso_irrep}/Operators'])
        
        vfl.PRINT_IB_INFO(0, the_iso_irrep, the_op_list)
        
        dis_irrep = the_fit_data.require_group(the_iso_irrep)
        
        if 'Operators' not in dis_irrep:
            dis_irrep.create_dataset('Operators', data = the_op_list)
        
        fit_group = dis_irrep.require_group(f'{the_fit_key}exp')
        the_tmin_data = fit_group.require_group('Tmin')

        the_cov_matrix = np.asarray(the_data[f'{the_iso_irrep}/Correlators/Real/Covariance_matrix'], dtype=np.float64)
        the_nt = np.asarray(the_data[f'{the_iso_irrep}/Time_slices'])
        the_eff_energy_hint = np.asarray(the_data[f'{the_iso_irrep}/Effective_masses/Mean'])
        
        ### Checking the tmin and tmax variables
        the_ul = int(the_list_tmaxs[0]) - the_nt[0]
        the_ll = np.arange(2, int(the_ul * 0.85))
    
        if type_correlated_fit == 'Correlated':
            the_cov_matrix_fit = the_cov_matrix
            if 'Correlated' in the_tmin_data:
                del the_tmin_data['Correlated']
            the_fit_data_group = the_tmin_data.create_group('Correlated')
        elif type_correlated_fit == 'Uncorrelated':
            the_cov_matrix_fit = np.diag(np.diag(the_cov_matrix))
            if 'Uncorrelated' in the_tmin_data:
                del the_tmin_data['Uncorrelated']
            the_fit_data_group = the_tmin_data.create_group('Uncorrelated')
        else:
            raise ValueError("Invalid fit correlation type")
        
        the_results = {'the_energies': [], 'the_sigmas': [], 'the_chi_vals': [], 'the_sigmas_chi': [], 'the_resampled': [], 'the_amplitudes': [],}
        
        begin_time = time.time()
        for the_yy in the_ll:
            print(f'Tmin = {the_yy + the_nt[0]} || TMax = {the_ul + the_nt[0]}')
            
            the_nt_slice = the_nt[the_yy:the_ul]
            the_corr_fit_slice = the_iso_corr[the_yy:the_ul]
            
            the_inverse_cov_m = np.linalg.inv(vfa.SHRINK_MATRIX(the_cov_matrix_fit, the_yy, the_ul))
            
            da_hint = vfa.BEST_GUESS(the_corr_fit_slice, the_nt_slice, the_fit_key)
            if np.any(np.isnan(da_hint)):
                the_dof = np.zeros(the_n_params_iso)
                the_dof[0] = 0.1
                the_dof[1] = the_eff_energy_hint[the_yy]
            else:
                the_dof = da_hint.copy()

            the_fit_choice = vfa.My_Fits(da_minimization_iso, the_nt_slice, the_corr_fit_slice, the_inverse_cov_m, the_dof, np.float64(0.))
            the_fit = Minuit(the_fit_choice, the_dof, name = the_fit_params_iso)

            the_fit.errordef, the_fit.tol = 1e-8, 1e-10
            the_fit.scan()
            the_fit.migrad(iterate=10, ncall=5000)
            
            e0 = np.double(the_fit.values['e0'])
            a0 = np.double(the_fit.values['a0'])
            
            the_results['the_energies'].append(e0)
            the_results['the_amplitudes'].append(a0)
            the_results['the_chi_vals'].append(np.double(the_fit.fval))            
            
            the_dof_rs = the_dof.copy()
            the_dof_rs[0], the_dof_rs[1] = a0, e0
            
            the_n_rs = len(the_iso_corr_rs)
            the_chi_vals_rs = np.empty(the_n_rs)
            the_resampled_vals = np.empty(the_n_rs + 1)
            the_resampled_vals[0] = e0

            for zz in range(the_n_rs):                
                corr_rs_slice = the_iso_corr_rs[zz, the_yy:the_ul]
                
                the_fit_choice_rs = vfa.My_Fits(da_minimization_iso, the_nt_slice, corr_rs_slice, the_inverse_cov_m, the_dof_rs, np.float64(0.))
                the_fit_rs = Minuit(the_fit_choice_rs, the_dof_rs, name = the_fit_params_iso)
    
                the_fit_rs.errordef, the_fit_rs.tol = 1e-8, 1e-7
                the_fit_rs.scan()
                the_fit_rs.migrad(iterate=10, ncall=5000)
                
                e0_rs = np.double(the_fit_rs.values['e0'])
                
                the_chi_vals_rs[zz] = the_fit_rs.fval
                the_resampled_vals[zz + 1] = e0_rs

            rs_vals = the_resampled_vals[1:]
            
            the_results['the_sigmas'].append(vfa.STD_DEV(rs_vals, np.mean(rs_vals), the_type_rs))
            the_results['the_sigmas_chi'].append(vfa.STD_DEV(the_chi_vals_rs, np.mean(the_chi_vals_rs), the_type_rs))
            the_results['the_resampled'].append(the_resampled_vals)
        
        the_fit_data_group.create_dataset('Resampled', data = np.asarray(the_results['the_resampled']))
        the_fit_data_group.create_dataset('Mean', data = np.array([the_ll + the_nt[0], [the_ul + the_nt[0]]*(len(the_ll)), the_results['the_energies'], the_results['the_sigmas'], the_results['the_chi_vals'], the_results['the_sigmas_chi'], the_results['the_amplitudes']]))
        
        print(f'Minimization {the_type_fit.get("label", "")}: DONE!')
        end_time = time.time()        
        print(f'Time taken: {round((end_time-begin_time)/60,2)} min')
        
    elif the_corr_choice=='ib' or the_corr_choice=='both':
        
        ### ib correctons
        the_ib = the_type_fit['ib'] 
        
        ### Choosing the configuration for the ib corrections part
        the_n_params_ib, da_minimization_ib, the_fit_params_ib = the_ib["n_params"], the_ib["model"], the_ib["params"]
        
        for jj in range(1,3):
            the_irrep = the_s_irreps[jj]
            
            ### List of operators of this irrep
            the_op_list = list(the_data[f'{the_irrep}/Operators'])
            vfl.PRINT_IB_INFO(0, the_irrep, the_op_list)
            
            dis_irrep = the_fit_data.require_group(the_irrep)
            
            if 'Operators' not in dis_irrep:
                dis_irrep.create_dataset('Operators', data = the_op_list)
            
            fit_group = dis_irrep.require_group(f'{the_fit_key}exp')
            the_tmin_data = fit_group.require_group('Tmin')
            
            ### Now let's look at the correlators 
            the_corr_fit = the_corr_temp[the_irrep]
            the_corr_fit_rs = the_corr_rs_temp[the_irrep]
            
            the_cov_matrix = np.asarray(the_data[f'{the_irrep}/Correlators/Real/Covariance_matrix'], dtype=np.float64)
            the_nt = np.asarray(the_data[f'{the_irrep}/Time_slices'])
            the_eff_energy_hint = np.asarray(the_data[f'{the_irrep}/Effective_masses/Mean'])
            
            ### Checking the tmin and tmax variables
            the_ul = int(the_list_tmaxs[jj]) - the_nt[0]
            the_ll = np.arange(2, int(the_ul * 0.85))
        
            if type_correlated_fit == 'Correlated':
                the_cov_matrix_fit = the_cov_matrix
                if 'Correlated' in the_tmin_data:
                    del the_tmin_data['Correlated']
                the_fit_data_group = the_tmin_data.create_group('Correlated')
            elif type_correlated_fit == 'Uncorrelated':
                the_cov_matrix_fit = np.diag(np.diag(the_cov_matrix))
                if 'Uncorrelated' in the_tmin_data:
                    del the_tmin_data['Uncorrelated']
                the_fit_data_group = the_tmin_data.create_group('Uncorrelated')
            else:
                raise ValueError("Invalid fit type")
            
            the_results = {'the_energies': [], 'the_sigmas': [], 'the_chi_vals': [], 'the_sigmas_chi': [], 'the_resampled': [], 'the_amplitudes': [],}
            
            begin_time = time.time()
            for the_yy in the_ll:
                print(f'Tmin = {the_yy + the_nt[0]} || TMax = {the_ul + the_nt[0]}')
                
                the_nt_slice = the_nt[the_yy:the_ul]
                the_corr_fit_slice = the_corr_fit[the_yy:the_ul]
                
                the_inverse_cov_m = np.linalg.inv(vfa.SHRINK_MATRIX(the_cov_matrix_fit, the_yy, the_ul))
                the_dof = np.zeros(the_n_params_ib)
                the_dof[0], the_dof[1] = 0.1, the_eff_energy_hint[the_yy]

                the_fit_choice = vfa.My_Fits(da_minimization_ib, the_nt_slice, the_corr_fit_slice, the_inverse_cov_m, the_dof, np.float64(0.))                
                the_fit = Minuit(the_fit_choice, the_dof, name = the_fit_params_ib)

                the_fit.errordef, the_fit.tol = 1e-8, 1e-10
                the_fit.scan()
                the_fit.migrad(iterate=10, ncall=5000)
                
                e0 = np.double(the_fit.values['e0'])
                a0 = np.double(the_fit.values['a0'])
                
                the_results['the_energies'].append(e0)
                the_results['the_chi_vals'].append(np.double(the_fit.fval)) 
                
                the_dof_rs = the_dof.copy()
                the_dof_rs[0], the_dof_rs[1] = a0, e0
                
                the_n_rs = len(the_iso_corr_rs)
                the_chi_vals_rs = np.empty(the_n_rs)
                the_resampled_vals = np.empty(the_n_rs + 1)
                the_resampled_vals[0] = e0

                for zz in range(the_n_rs):
                    corr_rs_slice = the_corr_fit_rs[zz, the_yy:the_ul]

                    the_fit_choice_rs = vfa.My_Fits(da_minimization_ib, the_nt_slice, corr_rs_slice, the_inverse_cov_m, the_dof_rs, np.float64(0.))
                    the_fit_rs = Minuit(the_fit_choice_rs, the_dof_rs, name = the_fit_params_ib)

                    the_fit_rs.errordef, the_fit_rs.tol = 1e-8,  1e-7
                    the_fit_rs.scan()
                    the_fit_rs.migrad(iterate=10, ncall=5000)
                    
                    the_chi_vals_rs[zz] = the_fit_rs.fval
                    e0_rs = np.double(the_fit_rs.values['e0'])

                    the_resampled_vals[zz + 1] = e0_rs

                rs_vals = the_resampled_vals[1:]
                
                the_results['the_sigmas'].append(vfa.STD_DEV(rs_vals, np.mean(rs_vals), the_type_rs))
                the_results['the_sigmas_chi'].append(vfa.STD_DEV(the_chi_vals_rs, np.mean(the_chi_vals_rs), the_type_rs))
                the_results['the_resampled'].append(the_resampled_vals)
            
            the_fit_data_group.create_dataset('Resampled', data = np.asarray(the_results['the_resampled']))
            
            the_fit_data_group.create_dataset('Mean', data = np.array([the_ll + the_nt[0], [the_ul + the_nt[0]]*(len(the_ll)), the_results['the_energies'], the_results['the_sigmas'], the_results['the_chi_vals'], the_results['the_sigmas_chi']]))
            
        print(f'Minimization {the_type_fit.get("label", "")}: DONE!')
        end_time = time.time()        
        print(f'Time taken: {round((end_time-begin_time)/60,2)} min')
    
    else:
        sys.exit("No valid choice.")
