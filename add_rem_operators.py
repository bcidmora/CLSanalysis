import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_analysis_functions as vfa

def OperatorsAnalysis(the_matrix_correlator_data, the_type_rs, the_irreps, the_t0_min, the_t0_max, the_location, the_chosen_op_list, **kwargs):
    
    print("                     OPERATORS ANALYSIS PROCESS \n")
    ### The list of total irreps
    the_m_irreps =  list(the_matrix_correlator_data.keys())
    
    ### Resampling scheme
    if the_type_rs=='jk':
        the_resampling_scheme = 'Jackknife'
    elif the_type_rs=='bt':
        the_resampling_scheme = 'Bootstrap'
    
    ### Getting the t0 min and t0 max to do the GEVP
    if (the_t0_min is None or the_t0_max is None) and kwargs.get('the_td') is None:
        sys.exit('Error: T0 min or T0 max not valid. \nQuitting.')
    
    ### What type of sorting of the eigenstates
    the_sorting = kwargs.get('sorting')
    the_sorting_map = {
        None : (vfa.SORTING_EIGENVALUES, "Sorting states based on Eigenvalues."),
        'eigenvals' : (vfa. SORTING_EIGENVALUES, "Sorting states based on Eigenvalues."),
        'vecs_fix' : (vfa.SORTING_EIGENVECTORS, "Sorting states by Eigenvectors with a fixed reference time slice."),
        'vecs_fix_norm' : (vfa.SORTING_EIGENVECTORS_NORMALIZED, "Sorting states by normalized Eigenvectors with a fixed reference time slice."),
        'vecs_var' : (vfa.SORTING_EIGENVECTORS_CHANGING_TSLICE, "Sorting states by Eigenvectors with a varying reference time slice."),
        'vecs_var_norm' : (vfa.SORTING_EIGENVECTORS_NORMALIZED_CHANGING_TSLICE, "Sorting states by normalized Eigenvectors with a varying reference time slice." ), 
        'vecs_var_rs_mean' : (vfa.SORTING_EIGENVECTORS_RS_MEAN, "Sorting states by Eigenvectors with a varying reference time slice. The resamples are sorted in the same way as the mean."),}
    
    the_sorting_process, the_msg = the_sorting_map.get(the_sorting, the_sorting_map[None])
    print(the_msg)
    
    the_rs_sorting = kwargs.get('rs_sorting')
    if the_rs_sorting is None:
        the_rs_sorting_process = vfa.SORTING_EIGENVALUES_NEW
    else:
        the_rs_sorting_process = the_sorting_process
        
    the_matrix_correlator_data_reduced = h5py.File(the_location,'w')    
    
    begin_time = time.time()
    for jj, the_irrep in enumerate(the_m_irreps):
        
        print('\n----------------------------------------------')
        print(f'     IRREP ({jj+1}/{len(the_irreps)}): {the_irrep}')
        
        ### The data to analyse. 
        this_data = the_matrix_correlator_data[the_irrep]
        
        ### The operators list and the time slices
        the_op_list, the_nt = list(this_data['Operators']), np.asarray(this_data['Time_slices'])
        
        ### This is the original correlator
        this_data_corr = this_data['Correlators/Real/']
        
        if len(the_chosen_op_list[jj])>=1:
        
            the_mod_op_list = [the_op_list[the_item] for the_item in the_chosen_op_list[jj]]

            the_mod_data = this_data_corr['Mean'][:, the_chosen_op_list[jj]][:, :,  the_chosen_op_list[jj]]            
            
            the_mod_data_rs = this_data_corr['Resampled'][:, :, the_chosen_op_list[jj]][:, :, :, the_chosen_op_list[jj]]
            the_mod_data_sigmas = this_data_corr['Sigmas'][the_chosen_op_list[jj], :]
            
            print(f"Size of the Correlation matrix: {the_mod_data.shape[-1]}x{the_mod_data.shape[-1]}\nTime slices: {the_nt[0]} - {the_nt[-1]}\nResampling data: {the_resampling_scheme} {the_mod_data_rs.shape[1]}\n----------------------------------------------")
            
            print('      OPERATORS LIST \n----------------------------------------------')
            for i in range(len(the_mod_op_list)):
                print(f'       {the_mod_op_list[i].decode('utf-8')}')
            
            g_i = the_matrix_correlator_data_reduced.create_group(the_irrep)
            g_i.create_dataset('Correlators/Real/Mean', data = the_mod_data)
            g_i.create_dataset('Correlators/Real/Resampled', data = the_mod_data_rs)
            g_i.create_dataset('Correlators/Real/Sigmas', data = the_mod_data_sigmas)
            g_i.create_dataset('Time_slices', data = the_nt)
            g_i.create_dataset('Operators', data = the_mod_op_list) 
            
            group_gevp = g_i.create_group('GEVP')
            
            vfa.DOING_THE_GEVP([the_t0_min, the_t0_max], the_nt, the_mod_data, the_mod_data_rs, the_type_rs, the_sorting, the_sorting_process, the_rs_sorting_process, group_gevp)
            
        elif len(the_chosen_op_list[jj])<1:
            print("OBS: No reduction in the operator's set.")
            
            print(f"Size of the Correlation matrix: {this_data_corr['Mean'].shape[-1]}x{this_data_corr['Mean'].shape[-1]}\nTime slices: {the_nt[0]} - {the_nt[-1]}\nResampling data: {the_resampling_scheme} {this_data_corr['Resampled'].shape[1]}\n----------------------------------------------")
            
            print('      OPERATORS LIST \n----------------------------------------------')
            for i in range(len(the_op_list)):
                print(f'       {the_op_list[i].decode('utf-8')}')
            
            g_i = the_matrix_correlator_data_reduced.create_group(the_irrep)
            g_i.create_dataset('Correlators/Real/Mean', data = this_data_corr['Mean'])
            g_i.create_dataset('Correlators/Real/Resampled', data = this_data_corr['Resampled'])
            g_i.create_dataset('Correlators/Real/Sigmas', data = this_data_corr['Sigmas'])
            g_i.create_dataset('Time_slices', data = the_nt)
            g_i.create_dataset('Operators', data = the_op_list) 
            
            group_gevp = g_i.create_group('GEVP')
            vfa.DOING_THE_GEVP([the_t0_min, the_t0_max], the_nt, this_data_corr['Mean'], this_data_corr['Resampled'], the_type_rs, the_sorting, the_sorting_process, the_rs_sorting_process, group_gevp)
        
    end_time = time.time()
    print(f'TIME TAKEN: {(end_time-begin_time)/60} mins')
            


### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------


if __name__== "__main__":
    print("Nothing to do here.")
    
    
