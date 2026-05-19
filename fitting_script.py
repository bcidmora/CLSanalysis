import numpy as np
import h5py
import time
from iminuit import Minuit
import sys
import os
import set_of_analysis_functions as vfa
import set_of_layout_functions as vfl

import warnings
warnings.filterwarnings('ignore')

the_fit_map = {
        '1': (2, vfa.SINGLE_EXPONENTIAL, ('a0', 'e0'), "Single Exponential Fit"),
        '2': (4, vfa.DOUBLE_EXPONENTIAL, ('a0', 'e0', 'a1','e1'), "Double Exponential Fit"),
        'g': (4, vfa.GEOMETRIC_FORM, ('a0', 'e0', 'b','m'), "Geometric series Fit" ),}


def FitSingleCorrelators(the_data, the_fit_data, the_type_rs, the_list_tmaxs, the_irreps, the_type_fit, type_correlated_fit, **kwargs): 
    
    ### If only one tmin needs to be done for the fitting
    the_only_one_tmin = kwargs.get('one_tmin')
    
    ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(the_irreps)
    
    ### The name of the irreps
    the_s_irreps = the_irreps[the_first_irrep:the_last_irrep]
    
    try:
        the_n_params, da_minimization, the_fit_params, the_string_fit = the_fit_map[the_type_fit]
        the_dof_template = np.zeros(the_n_params)
    except KeyError:
        raise ValueError(f"Invalid fit type: {the_type_fit}")
    
    print(f"                     FITTING: {the_string_fit}\n")
    
    
    np_nan, np_double, np_empty, np_mean = np.isnan, np.double, np.empty, np.mean

    BEST_GUESS = vfa.BEST_GUESS
    SHRINK_MATRIX = vfa.SHRINK_MATRIX
    My_Fits = vfa.My_Fits
    STD_DEV = vfa.STD_DEV
    MinuitCls = Minuit
    
    irrep_index_map = {r: i for i, r in enumerate(the_irreps)}
    
    ### Loop over the irreps
    for the_irrep in the_s_irreps:
        
        j = irrep_index_map[the_irrep]
        
        ### List of operators of this irrep
        the_op_list = list(the_data[f'{the_irrep}/Operators'])
        print('----------------------------------------------------------------------------------------')
        print(f'IRREP ({j+1}/{len(the_irreps)}): {the_irrep}\n   --->>   Operators list: ')
        
        for item in the_op_list: print(f'           {item.decode("utf-8")}')    
        print('----------------------------------------------------------------------------------------')
        
        ### Check if this part already exists
        dis_irrep = the_fit_data.require_group(the_irrep)
        
        if 'Operators' not in dis_irrep:
            dis_irrep.create_dataset('Operators', data=the_op_list)
        
        fit_key = f'{the_type_fit}exp'
        if fit_key not in dis_irrep:
            the_one_exp_fit = dis_irrep.create_group(fit_key)
        else:
            the_one_exp_fit = dis_irrep[fit_key]
        
        the_tmin_data = the_one_exp_fit.require_group('Tmin')
        
        the_corr = the_data[f'{the_irrep}/Correlators']
        the_corr_fit = np.asarray(the_corr['Real/Mean'], dtype=np.float64)
        the_corr_fit_rs = vfa.NT_TO_NCFGS(np.asarray(the_corr['Real/Resampled'], dtype=np.float64))
        
        the_cov_matrix = np.asarray(the_corr['Real/Covariance_matrix'], dtype=np.float64)
        the_eff_energy_hint = np.asarray(the_data[f'{the_irrep}/Effective_masses/Mean'])
        the_nt = np.asarray(the_data[f'{the_irrep}/Time_slices'])
        
        if the_only_one_tmin and kwargs.get('the_tmin') is not None:
            the_ul = int(the_list_tmaxs[j]) - the_nt[0]
            the_ll = [the_tmin - the_nt[0]]
        else:
            the_ul = int(the_list_tmaxs[j]) - the_nt[0]
            the_ll = np.arange(2, int(the_ul * 0.85))
        
        if type_correlated_fit == 'Correlated':
            the_fit_data_group = the_tmin_data.require_group('Correlated')
            the_cov_matrix_fit = the_cov_matrix
        elif type_correlated_fit == 'Uncorrelated':
            the_fit_data_group = the_tmin_data.require_group('Uncorrelated')
            the_cov_matrix_fit = np.diag(np.diag(the_cov_matrix))
        else:
            raise ValueError("Invalid fit type")

        the_results = {'the_energies': [], 'the_sigmas': [], 'the_chi_vals': [], 'the_sigmas_chi': [], 'the_resampled': []}
        
        begin_time = time.time()
        for the_yy in the_ll:
            print(f'Tmin = {the_yy + the_nt[0]} || TMax = {the_ul + the_nt[0]}')
            
            the_nt_slice = the_nt[the_yy:the_ul]
            the_corr_fit_slice = the_corr_fit[the_yy:the_ul]
            
            da_hint = vfa.BEST_GUESS(the_corr_fit_slice, the_nt_slice, the_type_fit)
            
            if not np.any(np.isnan(da_hint)):
                the_dof = da_hint
            else:
                the_dof = the_dof_template.copy()
                the_dof[0] = 0.1
                the_dof[1] = the_eff_energy_hint[the_yy]
                
            the_inverse_cov_m = np.linalg.inv(SHRINK_MATRIX(the_cov_matrix_fit, the_yy, the_ul))
            
            the_fit_choice = My_Fits(da_minimization, the_nt_slice, the_corr_fit_slice, the_inverse_cov_m, the_dof, np.float64(0.))   
            
            the_fit = MinuitCls(the_fit_choice, the_dof, name = the_fit_params)
            the_fit.errordef, the_fit.tol = 1e-8, 1e-10
            the_fit.scan()
            the_fit.migrad(iterate=10,ncall=5000)
            
            e0 = np_double(the_fit.values['e0'])
            
            the_results['the_energies'].append(e0)
            the_results['the_chi_vals'].append(the_fit.fval)
            
            the_dof_rs = the_dof.copy()
            the_dof_rs[0], the_dof_rs[1] = np_double(the_fit.values['a0']), e0          
            
            the_n_rs = len(the_corr_fit_rs)
            the_chi_vals_rs = np.empty(the_n_rs)
            the_resampled_vals = np.empty(the_n_rs + 1)
            the_resampled_vals[0] = e0

            for zz in range(the_n_rs):
                the_corr_rs_slice = the_corr_fit_rs[zz,the_yy:the_ul]
                
                the_fit_choice_rs = My_Fits(da_minimization, the_nt_slice, the_corr_rs_slice, the_inverse_cov_m, the_dof_rs, np_double(0.))
                
                the_fit_rs = MinuitCls(the_fit_choice_rs, the_dof, name = the_fit_params)
                
                the_fit_rs.errordef, the_fit_rs.tol = 1e-8, 1e-7
                the_fit_rs.scan()
                the_fit_rs.migrad(iterate=10, ncall=5000)
                
                e0_rs = np_double(the_fit_rs.values['e0'])
                the_chi_vals_rs[zz] = the_fit_rs.fval
                the_resampled_vals[zz + 1] = e0_rs

            the_rs_vals = the_resampled_vals[1:]
            
            the_results['the_sigmas'].append(vfa.STD_DEV(the_rs_vals, np.mean(the_rs_vals), the_type_rs))
            the_results['the_sigmas_chi'].append(vfa.STD_DEV(the_chi_vals_rs, np.mean(the_chi_vals_rs), the_type_rs))
            the_results['the_resampled'].append(the_resampled_vals)
            
        the_results['the_resampled'] = np.array(the_results['the_resampled'])
        
        vfl.REPLACE_DATASET(the_fit_data_group, 'Resampled', the_results['the_resampled'])
        vfl.REPLACE_DATASET(the_fit_data_group, 'Mean', np.array([the_ll + the_nt[0], [the_ul + the_nt[0]]*(len(the_ll)), the_results['the_energies'], the_results['the_sigmas'], the_results['the_chi_vals'], the_results['the_sigmas_chi']]))

        print(f'Minimization {the_string_fit}: DONE!')
        end_time = time.time()        
        print(f'Time taken: {round((end_time-begin_time)/60,2)} min')
     
    
    
def FitMultiCorrelators(the_data, the_fit_data, the_type_rs, the_list_tmaxs, the_irreps, the_type_fit, the_type_correlated_fit, **kwargs):
    
    ### Fits only one minimum time slice
    the_only_one_tmin = kwargs.get('one_tmin')
    
    ### Fits only one t0 chosen
    the_only_one_t0 = kwargs.get('one_t0')
    
    ### For the single-pivot diagonalization
    the_td = kwargs.get('the_td')
    
    ### The names of the irreps in this ensemble
    the_m_irreps = list(the_data.keys())
    
    ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(the_irreps)
    
    ### The name of the irreps
    the_m_irreps = the_irreps[the_first_irrep:the_last_irrep]
    
    ### If only one t0 was chosen, here it takes it
    if the_only_one_t0:
        the_t0_s = [int(kwargs.get('chosen_t0'))]
    else:
        the_t0_s = sorted([int(item[3:]) for item in list(the_data[f'{the_m_irreps[0]}/GEVP'])])
    
    try:
        the_n_params, da_minimization, the_fit_params, the_string_fit = the_fit_map[the_type_fit]
        the_dof = np.zeros(the_n_params)
    except KeyError:
        raise ValueError(f"Invalid fit type: {the_type_fit}")
    
    print(f"                     FITTING: {the_string_fit}\n")
    
    ### Loop over the irreps
    for the_irrep in the_m_irreps:           
         
        ### Getting the time slices of this dataset
        the_nt = np.asarray(the_data[f'{the_irrep}/Time_slices'])
        
        ### These are the operators used for the full matrix
        the_op_list = list(the_data[f'{the_irrep}/Operators'])
        
        print('----------------------------------------------------------------------------------------')
        print(f'IRREP ({the_irreps.index(the_irrep)+1}/{len(the_irreps)}): {the_irrep}\n   --->>   Operators list: ')
        for item in the_op_list: 
            print(f'           {item.decode("utf-8")}')
        print('----------------------------------------------------------------------------------------')
        
        ### This is creating a path for the irrep
        if the_irrep not in the_fit_data.keys(): 
            dis_irrep = the_fit_data.create_group(the_irrep)
        else: 
            dis_irrep = the_fit_data[the_irrep]
        
        ### It creates a path to store the operators
        if 'Operators' not in dis_irrep.keys():
            dis_irrep.create_dataset('Operators', data = the_op_list)
        
        ### Searching if the path for the data exists
        if f'{the_type_fit}exp' not in dis_irrep.keys(): 
            the_exp_fit = dis_irrep.create_group(f'{the_type_fit}exp')
        else:
            the_exp_fit = the_fit_data[f'{the_irrep}/{the_type_fit}exp']
            
        ### Checking the GEVP was done before
        if 'GEVP' in list(the_data[the_irrep].keys()):
            
            begin_time_tmin = time.time()
            for the_t0 in the_t0_s:
                print(f'T0 equal: {the_t0}')
                
                ### Checking if the path inside this file exists
                if f't0_{the_t0}' not in the_exp_fit.keys():   
                    t0_group = the_exp_fit.create_group(f't0_{the_t0}')
                    tmin_data = t0_group.create_group('Tmin')
                else:
                    t0_group = the_exp_fit[f't0_{the_t0}']
                    tmin_data = t0_group['Tmin']
                
                ### Retrieving data for the fits (eigenvalues and resamples)
                the_corr = the_data[f'{the_irrep}/GEVP/t0_{the_t0}']
                vfa.DOING_THE_FITTING(the_corr, the_nt, the_type_rs, the_irreps, the_irrep, tmin_data, the_type_correlated_fit, the_type_fit, the_only_one_tmin, the_t0, the_list_tmaxs, da_minimization, the_fit_params)
            end_time_tmin = time.time()
            print(f'Time taken: {round((end_time_tmin-begin_time_tmin)/60,2)} mins')
            print(f'Minimization {the_type_fit}: E vs Tmin DONE!')
        else: raise ValueError('GEVP needs to be calculated first :D')
            
                
def FitRatioMultiCorrelators(the_data, the_fit_data, the_type_rs, the_list_tmaxs, the_irreps, the_type_fit, the_type_correlated_fit, **kwargs):
    
    print("                     FITTING \n")
    
    ### Fits only one minimum time slice
    the_only_one_tmin = kwargs.get('one_tmin')
    
    ### Fits only one t0 chosen
    the_only_one_t0 = kwargs.get('one_t0')
    
    ### The names of the irreps in this ensemble
    the_m_irreps = list(the_data.keys())
    
    if the_only_one_t0:
        the_t0_s = [int(kwargs.get('chosen_t0'))]
    else:
        the_t0_s = sorted([int(item[3:]) for item in list(the_data[f'{the_m_irreps[0]}/GEVP'])])
    
    ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(the_m_irreps)
    
    ### The name of the irreps
    the_m_irreps = the_irreps[the_first_irrep:the_last_irrep]
            
    try:
        the_n_params, da_minimization, the_fit_params, the_string_fit = the_fit_map[the_type_fit]
        the_dof = np.zeros(the_n_params)
    except KeyError:
        raise ValueError(f"Invalid fit type: {the_type_fit}")
    
    ### Loop over the irreps
    for the_irrep in the_m_irreps:           
         
        ### Getting the time slices of this dataset
        the_nt = np.asarray(the_data[f'{the_irrep}/Time_slices'])
        
        ### These are the operators used for the full matrix
        the_op_list = list(the_data[f'{the_irrep}/Operators'])
        
        ### These are the non-interacting levels
        the_non_int = list(the_data[f'{the_irrep}/Single_hadron_corrs'])
        
        print('----------------------------------------------------------------------------------------')
        print(f'IRREP ({the_irreps.index(the_irrep)+1}/{len(the_irreps)}): {the_irrep}\n   --->>   Operators list: ')
        
        for item in the_op_list: 
            print(f'           {item.decode("utf-8")}')
        print('----------------------------------------------------------------------------------------')
        
        ### This is creating a path for the irrep
        if the_irrep not in the_fit_data.keys(): 
            dis_irrep = the_fit_data.create_group(the_irrep)
        else: 
            dis_irrep = the_fit_data[the_irrep]
        
       ### It creates a path to store the operators
        if 'Operators' not in dis_irrep.keys():
            dis_irrep.create_dataset('Operators', data = the_op_list)
        
        if 'Single_hadron_corrs' not in dis_irrep.keys():
            dis_irrep.create_dataset('Single_hadron_corrs', data = the_non_int)
        
        ### Searching if the path for the data exists
        if f'{the_type_fit}exp' not in dis_irrep.keys(): 
            the_exp_fit = dis_irrep.create_group(f'{the_type_fit}exp')
        else:
            the_exp_fit = the_fit_data[f'{the_irrep}/{the_type_fit}exp']
            
        begin_time_tmin = time.time()
        for the_t0 in the_t0_s:
            print(f'T0 equal: {the_t0}')
            
             ### Checking if the path inside this file exists
            if f't0_{the_t0}' not in the_exp_fit.keys():   
                t0_group = the_exp_fit.create_group(f't0_{the_t0}')
                tmin_data = t0_group.create_group('Tmin')
            else:
                t0_group = the_exp_fit[f't0_{the_t0}']
                tmin_data = t0_group['Tmin']
            
            ### Retrieving data for the fits (eigenvalues and resamples)
            the_corr = the_data[f'{the_irrep}/GEVP/t0_{the_t0}']
            
             ### Loop over all the t0 chosen
            the_corr_fit = np.asarray(the_corr['Eigenvalues/Mean']).real 
            the_corr_fit_rs = np.asarray(the_corr['Eigenvalues/Resampled']).real 
            
            ### Retrieving the covariance matrix for the fits
            the_cov_matrix = np.asarray(the_corr['Eigenvalues/Covariance_matrix']).real 
            
            ### Effective Masses are used as a prior
            the_eff_energy_hint = np.asarray(the_corr['Effective_masses/Mean'])
            
            ### Creating a folder for correlated or uncorrelated fits
            if the_type_correlated_fit=='Correlated':
                if 'Correlated' in tmin_data.keys(): del tmin_data['Correlated']
                the_fit_data_corr = tmin_data.create_group('Correlated')
            elif the_type_correlated_fit=='Uncorrelated':
                if 'Uncorrelated' in tmin_data.keys(): del tmin_data['Uncorrelated']
                the_fit_data_corr = tmin_data.create_group('Uncorrelated')
            
            ### Creating a folder for the fits of the central values
            the_mean_data = the_fit_data_corr.create_group('Mean')
            
            ### Creating a folder for the resampled data fits
            the_rs_data = the_fit_data_corr.create_group('Resampled')
            
            ### Loop over eigenvalues
            for ls in range(len(the_corr_fit)):

                ### Now also loop over each non-interacting combination
                for nn in range(the_corr_fit.shape[1]):
                    
                    ### This is the modified time range to go from item 0 until the end without caring about labeling
                    nt_mod = np.arange(0,len(the_nt))
                    if the_only_one_tmin: 
                        
                        ### Upper limit for the fit
                        the_ul = [x-the_nt[0] for x in the_list_tmaxs[the_irreps.index(the_irrep)]]
                        
                        ### Lower limit for the fit
                        the_ll = [nt_mod[5]]                    
                    else: 
                        ### Upper limir for the fit
                        the_ul = [x-the_nt[0] for x in the_list_tmaxs[the_irreps.index(the_irrep)]]
                        
                        ### Lower limit for the fit depends on the upper limit
                        the_ll = np.arange(nt_mod[0]+1, int(nt_mod[-1]*.75))
                    
                    ### Choosing the covariance matrix depending on the type of correlated fit
                    if the_type_correlated_fit=='Correlated':
                        the_cov_matrix_fit = the_cov_matrix[ls,nn,:,:]
                    elif the_type_correlated_fit=='Uncorrelated':
                        the_cov_matrix_fit = np.diag(np.diag(the_cov_matrix[ls,nn,:,:]))          
                    
                    the_results = {'the_energies': [], 'the_sigmas': [], 'the_chi_vals': [], 'the_sigmas_chi': [], 'the_resampled': []}
                    another_list = []
                    
                    the_corr_fit_slice = the_corr_fit[ls,nn]
                    the_corr_fit_rs_slice = the_corr_fit_rs[ls,nn]                    
                    
                    ### Loop over all the tmins
                    for yy in the_ll:
                        print(f'Tmin = {yy+the_nt[0]} || TMax = {the_ul[ls]+the_nt[0]}')
                        another_useful_list = []
                        
                        ### This is finding a good guess to make the fit converge easier.
                        the_dof[0] = np.float64(0.1)
                        the_dof[1] = np.float64(the_eff_energy_hint[ls,nn,yy])
                        
                        ### Reduces the covariance matrix to the size of the time range chosen and takes the inverse of it                        
                        the_inverse_cov_m = np.linalg.inv(vfa.SHRINK_MATRIX(the_cov_matrix_fit, yy, the_ul[ls]))
                        
                        ### This chooses the fit function to use 
                        the_fit_choice = vfa.My_Fits(da_minimization, the_nt[yy:the_ul[ls]], the_corr_fit_slice[yy:the_ul[ls]], the_inverse_cov_m, the_dof, np.float64(the_t0))
                        
                        ### Fitting started
                        the_fit = Minuit(the_fit_choice, the_dof, name = the_fit_params)
                        
                        the_fit.errordef, the_fit.tol = 1e-8, 1e-10
                        the_fit.scan()
                        the_fit.migrad(iterate=10,ncall=5000)
                        
                        ### Energy values
                        e0 = np.float128(the_fit.values['e0'])
                        
                        ### The fitted energy results from the central values are used as an initial guess for the resamples
                        the_dof_rs = the_dof
                        the_dof_rs[0], the_dof_rs[1] = np.float64(the_fit.values['a0']), e0
                        
                        another_useful_list.append(e0)
                        the_results['the_energies'].append(e0)
                        the_results['the_chi_vals'].append(np.float64(the_fit.fval))
                        
                        chi_vals_rs_list = []
                        ### Loop over the resamples
                        for zz in range(the_corr_fit_rs.shape[2]):
                            my_fit_choice_rs = vfa.My_Fits(da_minimization, the_nt[yy:the_ul[ls]], the_corr_fit_rs_slice[zz, yy:the_ul[ls]], the_inverse_cov_m, the_dof_rs, np.float64(the_t0))
                            the_fit_rs = Minuit(my_fit_choice_rs, the_dof_rs, name=the_fit_params)
                            
                            the_fit_rs.errordef, the_fit_rs.tol = 1e-8, 1e-7
                            the_fit_rs.scan()
                            the_fit_rs.migrad(iterate=10, ncall=5000)
                            
                            e0_rs = np.float64(the_fit_rs.values['e0']); 
                            chi_vals_rs_list.append(the_fit_rs.fval); 
                            another_useful_list.append(e0_rs)
                        
                        ### This is the sigma for the fittings
                        sigma_fit_rs = vfa.STD_DEV(another_useful_list[1:], np.mean(another_useful_list[1:]), the_type_rs)
                        sigma_chi_rs = vfa.STD_DEV(chi_vals_rs_list, np.mean(chi_vals_rs_list), the_type_rs)
                        
                        the_results['the_sigmas'].append(sigma_fit_rs)
                        the_results['the_sigmas_chi'].append(sigma_chi_rs)
                        another_list.append(np.array(another_useful_list))
                    print(f'E = {ls} READY')    

                    lambda_group_name = f"lambda_{ls}_nonint_{nn}"

                    the_rs_data.create_dataset(lambda_group_name, data = np.asarray(another_list))
                    the_mean_data.create_dataset(lambda_group_name, data = np.asarray([the_ll + the_nt[0], [the_ul[ls]+the_nt[0]]*len(the_ll), the_results['the_energies'], the_results['the_sigmas'], the_results['the_chi_vals'], the_results['the_sigmas_chi']]))                    
                    
                end_time_tmin = time.time()
                print(f'Time taken: {round((end_time_tmin-begin_time_tmin)/60,2)} min')
                print(f'Minimization {the_type_fit} exp: E vs Tmin DONE!')
            

if __name__== "__main__":
    print("Nothing to do here.")
