import numpy as np
from scipy.linalg import eigh
from scipy.linalg import fractional_matrix_power
import h5py
import time
import sys
import os
import set_of_analysis_functions as vfa

def EigenvaluesExtraction(the_matrix_correlator_data, the_type_rs, the_irreps, the_t0_min, the_t0_max, **kwargs):   
    
    print("                     SOLVING GEVP \n")
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
        the_rs_sorting_process = vfa.SORTING_EIGENVALUES
    else:
        the_rs_sorting_process = the_sorting_process

    begin_time = time.time()       
    for the_irrep in the_m_irreps:
        
        ### The data to analise is extracted here
        this_data = the_matrix_correlator_data[the_irrep]
        
        ### The list of operators of the correlation matrix and the time slices.
        the_op_list, the_nt = list(this_data['Operators']), np.asarray(this_data['Time_slices'])
        
        ### The size of each correlation matrix
        the_size_matrix = len(the_op_list)
        
        ### The resampled correlators
        the_rs_real = np.asarray(this_data['Correlators/Real/Resampled'])
        
        ### The central values of the original correlators
        the_mean_corr = np.asarray(this_data['Correlators/Real/Mean'], dtype = np.float64)    
        
        print('\n----------------------------------------------')
        print(f'     IRREP ({the_irreps.index(the_irrep)+1}/{len(the_irreps)}): {the_irrep}')
        print(f'Size of the Correlation matrix: {the_size_matrix}x{the_size_matrix}\nTime slices: {the_nt[0]} - {the_nt[-1]}\nResampling data ({the_resampling_scheme}): {the_rs_real.shape[1]}\n----------------------------------------------')
        print('      OPERATORS LIST \n----------------------------------------------')
        
        for i in range(the_size_matrix):
            print(f'       {the_op_list[i].decode('utf-8')}')
        
        if 'GEVP' in this_data.keys(): 
            del the_matrix_correlator_data[f'{the_irrep}/GEVP']
        group_gevp = this_data.create_group('GEVP')
        
        if kwargs.get('the_td')==None:
            ### This is a loop over the t0s
            vfa.DOING_THE_GEVP([the_t0_min, the_t0_max], the_nt, the_mean_corr, the_rs_real, the_type_rs, the_sorting, the_sorting_process, the_rs_sorting_process, group_gevp)
        else:
            vfa.DOING_THE_GEVP_SINGLE_PIVOT([the_t0_min, the_t0_max], the_nt, the_mean_corr, the_rs_real, the_type_rs, the_sorting, the_sorting_process, the_rs_sorting_process, group_gevp, t_diag=int(kwargs.get('the_td')))
    end_time = time.time()
    print(f'TIME TAKEN: {(end_time-begin_time)/60} mins')


### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------


if __name__== "__main__":
    print("Nothing to do here.")
