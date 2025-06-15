import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_functions as vf


def REMOVING_COLS_ROWS(the_matrix_correlator_data, the_type_rs, **kwargs):
    
    ### The list of total irreps
    the_list_name_irreps =  list(the_matrix_correlator_data.keys())
    
    ### If not all irreps are analysed, this number can be changed.
    if kwargs.get('nr_irreps')!=None:
        the_nr_irreps = int(kwargs.get('nr_irreps'))
    else:
        the_nr_irreps = len(the_list_name_irreps)    
    
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
        
    begin_time = time.time()
    for j in range(the_nr_irreps):
        
        ### The data to analyse. 
        this_data = the_matrix_correlator_data[the_list_name_irreps[j]]
        
        ### The operators list and the time slices
        the_op_list, the_nt = list(this_data.get('Operators')), np.array(this_data.get('Time_slices'))
        
        ### The size of the correlation matrix 
        the_size_matrix = len(the_op_list)
        
        ### If this analysis was done before, then that branch gets deleted and a new one is created.
        if 'Operators_Analysis' in this_data.keys(): del the_matrix_correlator_data[the_list_name_irreps[j]+'/Operators_Analysis']
        
        the_group_rows_cols = this_data.create_group('Operators_Analysis')
        
        ### This is just reshaping the data to make it easier to analyse
        the_mod_data = vf.RESHAPING_CORRELATORS(np.array(this_data.get('Correlators/Real/Mean')))
        the_mod_data_rs = vf.RESHAPING_CORRELATORS_RS_NT(np.array(this_data.get('Correlators/Real/Resampled')))
        
        print('\n----------------------------------------------')
        print('     IRREP (%s/'%str(j+1) + str(len(the_list_name_irreps)) +'): ', the_list_name_irreps[j])
            
        for ii in range(the_size_matrix):
            group_i = the_group_rows_cols.create_group('Op_%s'%str(ii))
            
            ### Removing one operator (column and row)
            the_new_corrs = vf.REMOVE_ROWS_COLS(the_mod_data,the_mod_data_rs,ii)
            
            ### The mean values of this new correlation matrix
            the_mean_corr = np.array(the_new_corrs[0], dtype=np.float128)
            
            ### The resampled values of this new correlation matrix
            the_rs_real = np.array(the_new_corrs[1], dtype=np.float128)
            
            print(the_mean_corr.shape,the_rs_real)
            
            the_mean_corr = vf.RESHAPING_EIGENVALS_MEAN(the_mean_corr,the_size_matrix-1)
            the_rs_real = vf.RESHAPING_EIGENVALS_RS(the_rs_real,the_size_matrix-1)
            
            print('Size of the Correlation matrix: ' + str(the_size_matrix-1)+ 'x' + str(the_size_matrix-1) +  '\nTime slices: '+str(the_nt[0])+' - '+str(the_nt[-1]) + '\nResampling data: %s '%the_type_rs+ '\n----------------------------------------------')
            print('      OPERATOR REMOVED \n----------------------------------------------')
            print('       '+str(the_op_list[ii].decode('utf-8')))
            
            the_t0_init=0
            for the_t0_init in range(np.abs(the_t0_min - the_nt[0]), (the_t0_max - the_nt[0]) + 1):                   
                the_ct0_mean = np.array(the_mean_corr[the_t0_init])
                
                the_evals_mean, the_evecs_mean, the_evecs_mean_ct0 = [], [], []
                the_evalues_rs, the_evectors_rs, the_evectors_rs_ct0 =[], [], []
                
                ttt=0
                for ttt in range(len(the_mean_corr)):
                    try:
                        the_evs_mean, the_evc_mean = eigh(the_mean_corr[ttt], b=the_ct0_mean, eigvals_only=False) #GEVP
                        the_evs_mean_nongevp, the_evec_mean_nongevp = eigh(the_mean_corr[ttt], eigvals_only=False) #NON-GEVP
                        the_evals_mean.append(the_evs_mean)
                        the_evecs_mean_ct0.append(the_evc_mean)
                        the_evecs_mean.append(the_evec_mean_nongevp)
                    except np.linalg.LinAlgError:
                        print("WARNING: Matrix isn't positive definite anymore. Skipping T0 = %s"%str(the_t0_init + the_nt[0]))
                        break
                the_evals_mean, the_evecs_mean = vf.SORTING_EIGENVALUES(the_t0_init, the_evals_mean, the_evecs_mean)
                if the_sorting!=None or the_sorting!='eigenvals':
                    the_evals_mean, the_evecs_mean = the_sorting_process(the_t0_init, the_evals_mean, the_evecs_mean)
                    
                ttt=0
                for ttt in range(len(the_mean_corr)):
                    the_evalues_rs_raw, the_evectors_rs_raw, the_evectors_rs_ct0_raw = [], [], []
                    xyz = 0
                    try:
                        for xyz in range(the_rs_real.shape[1]):
                            the_ct0 = np.array(the_rs_real[the_t0_init][xyz])
                            dis_resample = the_rs_real[ttt][xyz]
                            the_ew_rs, the_ev_rw = eigh(dis_resample, b=the_ct0, eigvals_only=False) #GEVP 
                            the_ew_rs_nongevp, the_ev_rw_nongevp = eigh(dis_resample, eigvals_only=False) #NON-GEVP
                            
                            the_evalues_rs_raw.append(np.array(the_ew_rs))
                            the_evectors_rs_ct0_raw.append(np.array(the_ev_rw, dtype=np.float128))
                            the_evectors_rs_raw.append(np.array(the_ev_rw_nongevp,dtype=np.float128))
                        the_evalues_rs.append(np.array(the_evalues_rs_raw))
                        the_evectors_rs_ct0.append(np.array(the_evectors_rs_ct0_raw))
                        the_evectors_rs.append(np.array(the_evectors_rs_raw))
                        
                    except np.linalg.LinAlgError: break
                    
                if len(evalues_rs)>0:
                    the_mod_evectors_rs = vf.RESHAPING_EIGEN_FOR_SORTING(np.array(the_evectors_rs))
                    the_mod_evals_rs = vf.RESHAPING_EIGEN_FOR_SORTING(np.array(the_evalues_rs))
                    for xyz in range(len(the_mod_evals_rs)):
                        the_mod_evals_rs[xyz], the_mod_evectors_rs[xyz] = vf.SORTING_EIGENVALUES(the_t0_init, the_mod_evals_rs[xyz], the_mod_evectors_rs[xyz])
                    the_evectors_rs = vf.RESHAPING_EIGEN_FOR_SORTING_REVERSE(the_mod_evectors_rs)
                    the_evalues_rs = vf.RESHAPING_EIGEN_FOR_SORTING_REVERSE(the_mod_evals_rs)
                    
                    group_t0 = group_i.create_group('t0_%s'%(the_t0_init+the_nt[0]))
                    
                    the_evecs_mean = np.array(the_evecs_mean)

                    the_eigevals_final_mean = vf.NT_TO_NCFGS(the_evals_mean)
                    
                    the_evals_fits_rs = np.array(vf.RESHAPING_EIGENVALS_FOR_FITS(np.array(evalues_rs), the_size_matrix-1), dtype=np.float128)
                    
                    group_eigvecs = group_t0.create_group('Eigenvectors')
                    group_eigvecs.create_dataset('Resampled', data=evectors_rs)
                    group_eigvecs.create_dataset('Mean', data=the_evecs_mean)
                    
                    the_l, the_sigma_2 = 0, []
                    for the_l in range(the_size_matrix):
                        dis_eign = vf.NCFGS_TO_NT(the_evals_fits_rs[the_l])
                        evals_fits_rs_mean = vf.MEAN(dis_eign)
                        the_sigma_2.append(vf.COV_MATRIX(dis_eign, evals_fits_rs_mean, the_type_rs))

                    group_eigns = group_t0.create_group('Eigenvalues')
                    group_eigns.create_dataset('Resampled', data = the_evals_fits_rs)
                    group_eigns.create_dataset('Mean', data = the_eigevals_final_mean)
                    group_eigns.create_dataset('Covariance_matrix', data = np.array(the_sigma_2))
                    print('T0 = %s'%str(the_t0_init + the_nt[0]) + '... DONE')
        j+=1
    end_time = time.time()
            


def ADDING_COLS_ROWS(the_matrix_correlator_data, the_type_rs, **kwargs):
    
    ### The list of total irreps
    the_list_name_irreps =  list(the_matrix_correlator_data.keys())
    
    ### Resampling scheme
    if the_type_rs=='jk':
        the_resampling_scheme = 'Jackknife'
    elif the_type_rs=='bt':
        the_resampling_scheme = 'Bootstrap'
    
    ### If not all irreps, then this can be changed
    if kwargs.get('nr_irreps')!=None:
        the_nr_irreps = int(kwargs.get('nr_irreps'))
    else:
        the_nr_irreps = len(the_list_name_irreps)    
    
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
        
    
    begin_time = time.time()
    for j in range(the_nr_irreps):
        this_data = the_matrix_correlator_data[the_list_name_irreps[j]]
        
        the_op_list, the_nt = list(this_data.get('Operators')), np.array(this_data.get('Time_slices'))
        the_size_matrix = len(the_op_list)
        
        if 'Operators_Analysis' in this_data.keys(): del the_matrix_correlator_data[the_list_name_irreps[j]+'/Operators_Analysis']
        
        the_group_rows_cols = this_data.create_group('Operators_Analysis')
        
        the_mod_data = vf.RESHAPING_CORRELATORS(np.array(this_data.get('Correlators/Real/Mean')))
        the_mod_data_rs = vf.RESHAPING_CORRELATORS_RS_NT(np.array(this_data.get('Correlators/Real/Resampled')))
        
        print('\n----------------------------------------------')
        print('     IRREP (%s/'%str(j+1) + str(len(the_list_name_irreps)) +'): ', the_list_name_irreps[j])
            
        for ii in range(the_size_matrix):
            group_i = the_group_rows_cols.create_group('Op_%s'%str(ii))
            
            the_new_corrs = vf.ADD_ROWS_COLS(the_mod_data,the_mod_data_rs,ii+1)
            the_mean_corr = np.array(the_new_corrs[0], dtype=np.float128)
            the_rs_real = np.array(the_new_corrs[1], dtype=np.float128)
            
            the_mean_corr = vf.RESHAPING_EIGENVALS_MEAN(the_mean_corr,len(the_mean_corr))
            the_rs_real = vf.RESHAPING_EIGENVALS_RS(the_rs_real,len(the_rs_real))
            
            print('Size of the Correlation matrix: ' + str(len(the_mean_corr))+ 'x' + str(len(the_mean_corr)) +  '\nTime slices: '+str(the_nt[0])+' - '+str(the_nt[-1]) + '\nResampling data (%s) '%the_resampling_scheme + str(the_rs_real.shape[-1]) + '\n----------------------------------------------')
            print('      OPERATOR ADDED \n----------------------------------------------')
            print('       '+str(the_op_list[ii].decode('utf-8')))

            the_t0_init=0
            for the_t0_init in range(np.abs(the_t0_min - the_nt[0]), (the_t0_max - the_nt[0]) + 1):                   
                the_ct0_mean = np.array(the_mean_corr[the_t0_init])
                
                the_evals_mean, the_evecs_mean, the_evecs_mean_ct0 = [], [], []
                the_evalues_rs, the_evectors_rs, the_evectors_rs_ct0 =[], [], []
                
                ttt=0
                for ttt in range(len(the_mean_corr)):
                    try:
                        the_evs_mean, the_evc_mean = eigh(the_mean_corr[ttt], b=the_ct0_mean, eigvals_only=False) # GEVP   
                        the_evs_mean_nongevp, the_evec_mean_nongevp = eigh(the_mean_corr[ttt], eigvals_only=False) #NON-GEVP
                        the_evals_mean.append(the_evs_mean)
                        the_evecs_mean_ct0.append(the_evc_mean)
                        the_evecs_mean.append(the_evec_mean_nongevp)
                    except np.linalg.LinAlgError:
                        print("WARNING: Matrix isn't positive definite anymore. Skipping T0 = %s"%str(the_t0_init + the_nt[0]))
                        break
                the_evals_mean, the_evecs_mean = vf.SORTING_EIGENVALUES(the_t0_init, the_evals_mean, the_evecs_mean)
                if the_sorting!=None or the_sorting!='eigenvals':
                    the_evals_mean, the_evecs_mean = the_sorting_process(the_t0_init, the_evals_mean, the_evecs_mean)
                
                ttt=0
                for ttt in range(len(the_mean_corr)):
                    the_evalues_rs_raw, the_evectors_rs_raw, the_evectors_rs_ct0_raw = [], [], []
                    xyz = 0
                    try:
                        for xyz in range(the_rs_real.shape[1]):
                            the_ct0 = np.array(the_rs_real[the_t0_init][xyz])
                            dis_resample = the_rs_real[ttt][xyz]
                            the_ew_rs, the_ev_rw = eigh(dis_resample, b=the_ct0, eigvals_only=False) #GEVP 
                            the_ew_rs_nongevp, the_ev_rw_nongevp = eigh(dis_resample, eigvals_only=False) #NON-GEVP
                            
                            the_evalues_rs_raw.append(np.array(the_ew_rs))
                            the_evectors_rs_ct0_raw.append(np.array(the_ev_rw, dtype=np.float128))
                            the_evectors_rs_raw.append(np.array(the_ev_rw_nongevp,dtype=np.float128))
                        the_evalues_rs.append(np.array(the_evalues_rs_raw))
                        the_evectors_rs_ct0.append(np.array(the_evectors_rs_ct0_raw))
                        the_evectors_rs.append(np.array(the_evectors_rs_raw))
                        
                    except np.linalg.LinAlgError: break
                
                if len(the_evalues_rs)>0:
                    
                    the_mod_evectors_rs = vf.RESHAPING_EIGEN_FOR_SORTING(np.array(the_evectors_rs))
                    the_mod_evals_rs = vf.RESHAPING_EIGEN_FOR_SORTING(np.array(the_evalues_rs))
                    
                    for xyz in range(len(the_mod_evals_rs)):
                        the_mod_evals_rs[xyz], the_mod_evectors_rs[xyz] = vf.SORTING_EIGENVALUES(the_t0_init, the_mod_evals_rs[xyz], the_mod_evectors_rs[xyz])
                    
                    the_evectors_rs = vf.RESHAPING_EIGEN_FOR_SORTING_REVERSE(the_mod_evectors_rs)
                    the_evalues_rs = vf.RESHAPING_EIGEN_FOR_SORTING_REVERSE(the_mod_evals_rs)
                    
                    group_t0 = group_i.create_group('t0_%s'%(the_t0_init+the_nt[0]))

                    the_eigevals_final_mean = vf.NT_TO_NCFGS(the_evals_mean)
                    the_evals_fits_rs = np.array(vf.RESHAPING_EIGENVALS_FOR_FITS(np.array(the_evalues_rs), the_evalues_rs.shape[-1]), dtype=np.float128)
                    
                    l, the_sigma_2 = 0, []
                    for l in range(the_evals_fits_rs.shape[0]):
                        dis_eign = vf.NCFGS_TO_NT(the_evals_fits_rs[l])
                        the_evals_fits_rs_mean = vf.MEAN(dis_eign)
                        the_sigma_2.append(vf.COV_MATRIX(dis_eign, the_evals_fits_rs_mean, the_type_rs))

                    the_evecs_mean_ct0 = np.array(the_evecs_mean_ct0)
                    # the_evecs_mean = np.array(the_evecs_mean)
                    
                    group_eigvecs = group_t0.create_group('Eigenvectors')
                    group_eigvecs.create_dataset('Mean', data=the_evecs_mean_ct0)       
                    # group_eigvecs.create_dataset('Mean', data=the_evecs_mean)
                    group_eigvecs.create_dataset('Resampled', data=the_evectors_rs_ct0)
                    # group_eigvecs.create_dataset('Resampled', data=the_evectors_rs)
                    
                    group_eigns = group_t0.create_group('Eigenvalues')
                    group_eigns.create_dataset('Resampled', data = the_evals_fits_rs)
                    group_eigns.create_dataset('Mean', data = the_eigevals_final_mean)
                    group_eigns.create_dataset('Covariance_matrix', data = np.array(the_sigma_2))
                    print('T0 = %s'%str(the_t0_init + the_nt[0]) + '... DONE')
        j+=1
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
    
    # REMOVING_COLS_ROWS(myMatrixCorrelatorData, myTypeRs)#, nr_irreps=myNrIrreps)
    ADDING_COLS_ROWS(myMatrixCorrelatorData, myTypeRs)#, nr_irreps=myNrIrreps)

    
    myArchivo.close()
    # print('-'*(len(savedLocation)+1))
    # print('Saved as: \n' + savedLocation)
    # print('_'*(len(savedLocation)+1))
    
    
