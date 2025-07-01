import numpy as np
from scipy.linalg import eigh
from scipy.linalg import fractional_matrix_power
import h5py
import time
import sys
import os
import set_of_functions as vf


def EigenvaluesExtraction(the_matrix_correlator_data, the_type_rs, the_irreps, **kwargs):   
    
    ### The list of total irreps
    the_m_irreps =  list(the_matrix_correlator_data.keys())
    
    ### Resampling scheme
    if the_type_rs=='jk':
        the_resampling_scheme = 'Jackknife'
    elif the_type_rs=='bt':
        the_resampling_scheme = 'Bootstrap'
    
    ### Getting the t0 min and t0 max to do the GEVP
    if kwargs.get('t0_min')==None or kwargs.get('t0_max')==None:
        print('Error: T0 min or T0 max not valid.')
        sys.exit('Quitting.')
    else:
        the_t0_min, the_t0_max = int(kwargs.get('t0_min')), int(kwargs.get('t0_max'))
    
    ### What type of sorting of the eigenstates
    the_sorting = kwargs.get('sorting')
    if the_sorting==None or the_sorting=='eigenvals':
        print("Sorting states based on Eigenvalues.")
        ### This function returns the eigenvalues sorted from the largest to the smallest.
        the_sorting_process = vf.SORTING_EIGENVALUES
        
    elif the_sorting=='vecs_fix':
        print("Sorting states by Eigenvectors with a fixed reference time slice.")
        ### This function returns the eigenvalues sorted based on the orthogonality of the eigenvectors based on a reference time slice where the eigenstated are already sorted.
        the_sorting_process = vf.SORTING_EIGENVECTORS        
        
    elif the_sorting=='vecs_fix_norm':
        ### This function returns the eigenvalues sorted based on the orthogonality of the normalized eigenvectors based on a reference time slice where the eigenstated are already sorted.
        print("Sorting states by normalized Eigenvectors with a fixed reference time slice.")
        the_sorting_process = vf.SORTING_EIGENVECTORS_NORMALIZED
        
    elif the_sorting=='vecs_var':
        print("Sorting states by Eigenvectors with a varying reference time slice.")
        ### This function returns the eigenvalues sorted based on the orthogonality of the eigenvectors based on the previous reference time slice.
        the_sorting_process = vf.SORTING_EIGENVECTORS_CHANGING_TSLICE
        
    elif the_sorting=='vecs_var_norm':
        print("Sorting states by normalized Eigenvectors with a varying reference time slice.")
        ### This function returns the eigenvalues sorted based on the orthogonality of the normalized eigenvectors based on the previous reference time slice.
        the_sorting_process = vf.SORTING_EIGENVECTORS_NORMALIZED_CHANGING_TSLICE

    begin_time = time.time()       
    for the_irrep in the_m_irreps:
        
        ### The data to analise is extracted here
        this_data = the_matrix_correlator_data[the_irrep]
        
        ### The list of operators of the correlation matrix and the time slices.
        the_op_list, the_nt = list(this_data.get('Operators')), np.array(this_data.get('Time_slices'))
        
        ### The size of each correlation matrix
        the_size_matrix = len(the_op_list)
        
        ### The resampled correlators
        the_rs_real = np.array(this_data.get('Correlators/Real/Resampled'))
        
        ### The central values of the original correlators
        the_mean_corr_real = np.array(this_data.get('Correlators/Real/Mean'))    
        the_mean_corr = np.array(the_mean_corr_real, dtype=np.float64)
        
        print('\n----------------------------------------------')
        print('     IRREP (%s/'%str(the_irreps.index(the_irrep)+1) + str(len(the_irreps)) +'): ', the_irrep)
        print('Size of the Correlation matrix: ' + str(the_size_matrix)+ 'x' + str(the_size_matrix) +  '\nTime slices: '+str(the_nt[0])+' - '+str(the_nt[-1]) + '\nResampling data (%s): '%the_resampling_scheme+ str(the_rs_real.shape[1]) +  '\n----------------------------------------------')
        print('      OPERATORS LIST \n----------------------------------------------')
        for i in range(the_size_matrix):
            print('       '+str(the_op_list[i].decode('utf-8')))
        
        if 'GEVP' in this_data.keys(): del the_matrix_correlator_data[the_irrep+'/GEVP']
        group_gevp = this_data.create_group('GEVP')
        
        ### This is a loop over the t0s
        vf.DOING_THE_GEVP([the_t0_min, the_t0_max], the_nt, the_mean_corr, the_rs_real, the_type_rs, the_sorting, the_sorting_process, group_gevp)
    end_time = time.time()
    print('TIME TAKEN: ' + str((end_time-begin_time)/60) +' mins')




### ------------------------------- END FUNCTIONS ----------------------------------------------------


### --------------------------------------------------------------------------------------------------



### ------------------------------- START EXECUTING --------------------------------------------------


if __name__=="__main__":
    
    myEns = str(sys.argv[1]).upper()
    myTypeRs = str(sys.argv[2]).lower()
    myRebinOn = str(sys.argv[3]).lower()
    myRb = sys.argv[4]
    
    mySorting = 'eigenvals' #'eigenvals' # 'vecs_fix' # 'vecs_fix_norm' # 'vecs_var' # 'vecs_var_norm'
    myVersion = 'test'
    
    myT0Min = int(input('T0 min: '))
    myT0Max = int(input('T0 max: ')) 
    
    myLocation = os.path.expanduser('~')+'$YOUR_OUTPUT_PATH(SAME_THAN_CORRS_SCRIPT)$/%s/'%myEns
    
    if myRebinOn=='rb': 
        rb = int(myRb)
        reBin = '_bin'+str(rb)
    else:
        reBin = ''  
    
    myArchivo = myLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion
    
    myMatrixCorrelatorData = h5py.File(myArchivo,'r+')
    
    EigenvaluesExtraction(myMatrixCorrelatorData, myTypeRs, t0_min = myT0Min, t0_max = myT0Max, sorting=mySorting )
    myMatrixCorrelatorData.close()   
    
    print('-'*(len(myArchivo)+1))
    print('Saved as: \n' + myArchivo)
    print('_'*(len(myArchivo)+1))
