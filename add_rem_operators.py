import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_functions as vf

def OperatorsAnalysis(the_matrix_correlator_data, the_type_rs, the_operator_analysis_method, the_irreps, **kwargs):
    
    ### The list of total irreps
    the_m_irreps =  list(the_matrix_correlator_data.keys())
    
    ### Resampling scheme
    if the_type_rs=='jk':
        the_resampling_scheme = 'Jackknife'
    elif the_type_rs=='bt':
        the_resampling_scheme = 'Bootstrap'
        
    ### If not all irreps are analysed, this number can be changed.
    if kwargs.get('nr_irreps')!=None:
        the_nr_irreps = int(kwargs.get('nr_irreps'))
    else:
        the_nr_irreps = len(the_m_irreps)    
    
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
        the_sorting_process = vf.SORTING_EIGENVALUES
    elif the_sorting=='vecs_fix':
        print("Sorting states by Eigenvectors with a fixed reference time slice.")
        the_sorting_process = vf.SORTING_EIGENVECTORS        
    elif the_sorting=='vecs_fix_norm':
        print("Sorting states by normalized Eigenvectors with a fixed reference time slice.")
        the_sorting_process = vf.SORTING_EIGENVECTORS_NORMALIZED
    elif the_sorting=='vecs_var':
        print("Sorting states by Eigenvectors with a varying reference time slice.")
        the_sorting_process = vf.SORTING_EIGENVECTORS_CHANGING_TSLICE
    elif the_sorting=='vecs_var_norm':
        print("Sorting states by normalized Eigenvectors with a varying reference time slice.")
        the_sorting_process = vf.SORTING_EIGENVECTORS_NORMALIZED_CHANGING_TSLICE
        
    ### If not all irreps are analysed, this number can be changed.
    if the_operator_analysis_method=='adding':
        the_op_method = vf.ADD_ROWS_COLS
        the_start_op = 1
        the_last_op = -1
        the_name_ops_analysis = 'Add_'
        print("OPERATOR ANALYSIS: Adding operators to the correlation matrix one by one starting with a 2x2 matrix.")
    elif the_operator_analysis_method=='removing':
        the_op_method = vf.REMOVE_ROWS_COLS
        the_start_op = 0
        the_last_op = 0
        the_name_ops_analysis = 'Remove_'
        print("OPERATOR ANALYSIS: Removing operators from the correlation matrix one at the time.")
    elif the_operator_analysis_method=='from_list':
        the_op_method = vf.CHOOSE_OPS
        the_chosen_op_list = kwargs.get('ops_analysis_list') # This list contains the operators that stay
        if the_chosen_op_list==None or len(the_chosen_op_list)==0:
           sys.exit('Exit Error: No list of operators chosen, please choose one.')
        
        
    begin_time = time.time()
    for the_irrep in the_m_irreps:
        
        ### The data to analyse. 
        this_data = the_matrix_correlator_data[the_irrep]
        
        ### The operators list and the time slices
        the_op_list, the_nt = list(this_data.get('Operators')), np.array(this_data.get('Time_slices'))
        
        if 'Operators_Analysis' not in this_data.keys(): 
            the_group_rows_cols = this_data.create_group('Operators_Analysis')
        else:
            the_group_rows_cols = this_data.get('Operators_Analysis')
        
         ### This is just reshaping the data to make it easier to analyse
        the_mod_data = vf.RESHAPING_CORRELATORS(np.array(this_data.get('Correlators/Real/Mean')))
        the_mod_data_rs = vf.RESHAPING_CORRELATORS_RS_NT(np.array(this_data.get('Correlators/Real/Resampled')))
        
        print('\n----------------------------------------------')
        print('     IRREP (%s/'%str(the_irreps.index(the_irrep)+1) + str(len(the_irreps)) +'): ', the_irrep)
        
        ### If you want to choose from a list, then you have to provide the list.
        if the_operator_analysis_method=='from_list':
            
            ### Here is checking if the list of operators chosen is the same than the original to not repeat the process
            if len(the_chosen_op_list[the_irreps.index(the_irrep)])<len(the_op_list) and len(the_chosen_op_list[the_irreps.index(the_irrep)])!=0:
                
                the_chosen_op_list_j = the_chosen_op_list[the_irreps.index(the_irrep)]
                
                ### The string has the positions of the chosen operators from the original full list
                the_ops_chosen_string = 'Ops_chosen_'
                for ss in the_chosen_op_list_j:
                    the_ops_chosen_string+=str(ss)
                
                if the_ops_chosen_string in this_data['Operators_Analysis'].keys(): del the_matrix_correlator_data[the_irrep+'/Operators_Analysis/'+the_ops_chosen_string]
                
                group_i = the_group_rows_cols.create_group(the_ops_chosen_string)
                
                ### This is getting the new correlation matrix with the chosen operators
                the_new_corrs = the_op_method(the_mod_data,the_mod_data_rs, the_chosen_op_list_j)
                
                ### The mean values of this new correlation matrix
                the_mean_corr = np.array(the_new_corrs[0], dtype=np.float64)
                
                ### The resampled values of this new correlation matrix
                the_rs_real = np.array(the_new_corrs[1], dtype=np.float64)
                
                ### Reshaping the datasets to do the GEVP
                the_mean_corr = vf.RESHAPING_EIGENVALS_MEAN(the_mean_corr)
                the_rs_real = vf.RESHAPING_EIGENVALS_RS(the_rs_real)
                
                print('Size of the Correlation matrix: ' + str(the_mean_corr.shape[-1])+ 'x' + str(the_mean_corr.shape[-1]) +  '\nTime slices: '+str(the_nt[0])+' - '+str(the_nt[-1]) + '\nResampling data: %s '%the_resampling_scheme + str(the_rs_real.shape[1]) + '\n----------------------------------------------')
                
                vf.DOING_THE_GEVP([the_t0_min, the_t0_max], the_nt, the_mean_corr, the_rs_real, the_type_rs, the_sorting, the_sorting_process, group_i)
                
        else:
            
            the_chosen_op_list_j = np.arange(the_start_op,len(the_op_list)+the_last_op)
             
            ### Loop over the operators
            for ii in the_chosen_op_list_j:
                ### Checking that the branch for this analysis was created before
                if (the_name_ops_analysis + 'Op_%s'%str(ii)) in this_data['Operators_Analysis'].keys(): 
                    del the_matrix_correlator_data[the_irrep+'/Operators_Analysis/' + the_name_ops_analysis + 'Op_%s'%str(ii)]
            
                ### It creates now the branch for this specific operator
                group_i = the_group_rows_cols.create_group(the_name_ops_analysis + 'Op_%s'%str(ii))
                
                ### Removing one operator (column and row)
                the_new_corrs = the_op_method(the_mod_data,the_mod_data_rs,ii)
                
                ### The mean values of this new correlation matrix
                the_mean_corr = np.array(the_new_corrs[0], dtype=np.float64)
                
                ### The resampled values of this new correlation matrix
                the_rs_real = np.array(the_new_corrs[1], dtype=np.float64)
                
                 ### Reshaping the datasets
                the_mean_corr = vf.RESHAPING_EIGENVALS_MEAN(the_mean_corr)
                the_rs_real = vf.RESHAPING_EIGENVALS_RS(the_rs_real)
                
                print('Size of the Correlation matrix: ' + str(the_mean_corr.shape[-1])+ 'x' + str(the_mean_corr.shape[-1]) +  '\nTime slices: '+str(the_nt[0])+' - '+str(the_nt[-1]) + '\nResampling data (%s): '%the_resampling_scheme + str(the_rs_real.shape[1]) + '\n----------------------------------------------')
        
                vf.DOING_THE_GEVP([the_t0_min, the_t0_max], the_nt, the_mean_corr, the_rs_real, the_type_rs, the_sorting, the_sorting_process, group_i)
    end_time = time.time()
            

                


if __name__=="__main__":
    
    myEns = str(sys.argv[1]).upper()
    myWhichCorrelator = "m"
    myTypeRs = str(sys.argv[2]).lower()
    myRebinOn = str(sys.argv[3]).lower()
    myRb = 1
    myVersion = 'test'
    myKbt = 500
    myNrIrreps=1
    
    if myRebinOn=='rb': 
        rb = int(myRb)
        reBin = '_bin'+str(rb)
    else:
        reBin = '' 
    
    if myEns == 'N451': from files_n451 import *
    elif myEns == 'N201': from files_n201 import * 
    elif myEns == 'D200': from files_d200 import *
    elif myEns == 'X451': from files_x451 import *
    
    myWeight = weight
    myLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'$YOUR_OUTPUT_PATH(TOTAL_SAME_THAN_CORRS_SCRIPT)$/%s/'%myEns)
    
    myCnfgs = ncfgs
    
    vf.INFO_PRINTING(myWhichCorrelator, myEns)
    
    myArchivo = myLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion
    
    myMatrixCorrelatorData = h5py.File(myArchivo,'r+')
    
    OperatorsAnalysis(myMatrixCorrelatorData, myTypeRs, myOperatorMethod, ops_analysis_list = myListOperators)#, nr_irreps=myNrIrreps)

    
    myArchivo.close()
    # print('-'*(len(savedLocation)+1))
    # print('Saved as: \n' + savedLocation)
    # print('_'*(len(savedLocation)+1))
    
    
