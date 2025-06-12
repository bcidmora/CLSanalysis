import numpy as np
from scipy.linalg import eigh
import h5py
import time
import sys
import os
import set_of_functions as vf



def EigenvaluesExtraction(the_matrix_correlator_data, the_type_rs, **kwargs):   
    
    ### The list of total irreps
    the_list_name_irreps =  list(the_matrix_correlator_data.keys())
    
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
    for j in range(len(the_list_name_irreps)):        
        this_data = the_matrix_correlator_data[the_list_name_irreps[j]]
        the_op_list, nt = list(this_data.get('Operators')), np.array(this_data.get('Time_slices'))
        the_size_matrix = len(the_op_list)
        
        the_rs_real = np.array(this_data.get('Correlators/Real/Resampled'))
        the_mean_corr_real = np.array(this_data.get('Correlators/Real/Mean'))    
        the_mean_corr = np.array(the_mean_corr_real, dtype=np.float128)
        
        print('\n----------------------------------------------')
        print('     IRREP (%s/'%str(j+1) + str(len(the_list_name_irreps)) +'): ', the_list_name_irreps[j])
        print('Size of the Correlation matrix: ' + str(the_size_matrix)+ 'x' + str(the_size_matrix) +  '\nTime slices: '+str(nt[0])+' - '+str(nt[-1]) + '\nResampling data (%s): '%the_resampling_scheme+ str(the_rs_real.shape[1]) +  '\n----------------------------------------------')
        print('      OPERATORS LIST \n----------------------------------------------')
        for i in range(the_size_matrix):
            print('       '+str(the_op_list[i].decode('utf-8')))
        
        if 'GEVP' in this_data.keys(): del the_matrix_correlator_data[the_list_name_irreps[j]+'/GEVP']
        group_gevp = this_data.create_group('GEVP')
        
        the_t0_init=0
        for the_t0_init in range(np.abs(the_t0_min - nt[0]), (the_t0_max - nt[0]) + 1):               
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
                    print("WARNING: Matrix isn't positive definite anymore. Skipping T0 = %s"%str(the_t0_init + nt[0]))
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
                
                group_t0 = group_gevp.create_group('t0_%s'%(the_t0_init+nt[0]))
                
                the_eigevals_final_mean = vf.NT_TO_NCFGS(the_evals_mean)
                the_evals_fits_rs = np.array(vf.RESHAPING_EIGENVALS_FOR_FITS(np.array(the_evalues_rs), the_size_matrix), dtype=np.float128)
                
                group_eigvecs = group_t0.create_group('Eigenvectors')
                group_eigvecs.create_dataset('Resampled', data=the_evectors_rs_ct0)
                
                l, the_sigma_2 = 0, []
                for l in range(the_size_matrix):
                    dis_eign = vf.NCFGS_TO_NT(the_evals_fits_rs[l])
                    evals_fits_rs_mean = vf.MEAN(dis_eign)
                    the_sigma_2.append(vf.COV_MATRIX(dis_eign, evals_fits_rs_mean, the_type_rs))
                
                the_evecs_mean_ct0 = np.array(the_evecs_mean_ct0)
                
                group_eigvecs.create_dataset('Mean', data=the_evecs_mean_ct0)
                group_eigns = group_t0.create_group('Eigenvalues')
                group_eigns.create_dataset('Resampled', data = the_evals_fits_rs)
                group_eigns.create_dataset('Mean', data = the_eigevals_final_mean)
                group_eigns.create_dataset('Covariance_matrix', data = np.array(the_sigma_2))
                print(".......................")
                print('T0 = %s'%str(the_t0_init + nt[0]) + '... DONE')
        j+=1
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
    
    EigenvaluesExtraction(myMatrixCorrelatorData, myTypeRs, t0_min = myT0Min, t0_max = myT0Max)
    myMatrixCorrelatorData.close()   
    
    print('-'*(len(myArchivo)+1))
    print('Saved as: \n' + myArchivo)
    print('_'*(len(myArchivo)+1))
