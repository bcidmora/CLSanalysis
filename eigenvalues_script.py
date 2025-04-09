import numpy as np
from scipy.linalg import eigh
import h5py
import time
import sys
import os
import set_of_functions as vf


def EigenvaluesExtraction(the_matrix_correlator_data, the_type_rs, **kwargs):   
    
    the_list_name_irreps =  list(the_matrix_correlator_data.keys())
    
    if kwargs.get('t0_min')==None or kwargs.get('t0_max')==None:
        print('Error: T0 min or T0 max not valid.')
        sys.exit('Quitting.')
    else:
        t0_min, t0_max = int(kwargs.get('t0_min')), int(kwargs.get('t0_max'))

    begin_time = time.time()
    for j in range(len(the_list_name_irreps)):        
        this_data = the_matrix_correlator_data[the_list_name_irreps[j]]
        the_op_list, nt = list(this_data.get('Operators')), np.array(this_data.get('Time_slices'))
        the_size_matrix = len(the_op_list)
        
        rs_real = np.array(this_data.get('Correlators/Real/Resampled'))
        mean_corr_real = np.array(this_data.get('Correlators/Real/Mean'))    
        mean_corr = np.array(mean_corr_real, dtype=np.float64)
        
        print('\n----------------------------------------------')
        print('     IRREP (%s/'%str(j+1) + str(len(the_list_name_irreps)) +'): ', the_list_name_irreps[j])
        print('Size of the Correlation matrix: ' + str(the_size_matrix)+ 'x' + str(the_size_matrix) +  '\nTime slices: '+str(nt[0])+' - '+str(nt[-1]) + '\nResampling data: %s '%the_type_rs+ '\n----------------------------------------------')
        print('      OPERATORS LIST \n----------------------------------------------')
        for i in range(the_size_matrix):
            print('       '+str(the_op_list[i].decode('utf-8')))
        
        if 'GEVP' in this_data.keys(): del the_matrix_correlator_data[the_list_name_irreps[j]+'/GEVP']
        group_gevp = this_data.create_group('GEVP')
        
        t0_init=0
        for t0_init in range(np.abs(t0_min - nt[0]), (t0_max - nt[0]) + 1):               
            my_ct0_mean = np.array(mean_corr[t0_init])
            
            evals_mean, evecs_mean, evecs_mean_ct0 = [], [], []
            evalues_rs, evectors_rs, evectors_rs_ct0 =[], [], []
            
            ttt=0
            for ttt in range(len(mean_corr)):
                try:
                    evs_mean, evc_mean = eigh(mean_corr[ttt], b=my_ct0_mean, eigvals_only=False, type=1, driver='gv')   
                    if ttt>t0_init:
                        evals_mean.append(np.flip(np.array(evs_mean)))
                        evec = np.flip(evc_mean); ppp=0
                        for ppp in range(len(evec)):
                            evec[ppp] = np.array(np.flip(evec[ppp]))
                        evc_mean = evec
                    else:
                        evals_mean.append(np.array(evs_mean))
                        evc_mean = np.array(evc_mean)
                    evecs_mean.append(evc_mean)
                except np.linalg.LinAlgError:
                    print("WARNING: Matrix isn't positive definite anymore. Skipping T0 = %s"%str(t0_init + nt[0]))
                    break

            ttt=0
            for ttt in range(len(mean_corr)):
                evalues_rs_raw, evectors_rs_raw, evectors_rs_ct0_raw = [], [], []
                xyz = 0
                try:
                    for xyz in range(rs_real.shape[1]):
                        my_ct0 = np.array(rs_real[t0_init][xyz])
                        dis_resample = rs_real[ttt][xyz]
                        ew_rs, ev_rw = eigh(dis_resample, b= my_ct0, eigvals_only=False,type=1,driver='gv')
                        if ttt>t0_init: 
                            evalues_rs_raw.append(np.flip(np.array(ew_rs)))
                            evecs = np.flip(ev_rw); ppp=0
                            for ppp in range(len(evecs)):
                                evecs[ppp] = np.array(np.flip(evecs[ppp]), dtype=np.float64)
                            evecs_flip = evecs
                        else: 
                            evalues_rs_raw.append(np.array(ew_rs))
                            evecs_flip = np.array(ev_rw, dtype=np.float64)
                        evectors_rs_raw.append(np.array(evecs_flip))
                    evalues_rs.append(np.array(evalues_rs_raw)); 
                    evectors_rs.append(np.array(evectors_rs_raw))
                except np.linalg.LinAlgError: break
            
            if len(evalues_rs)>0:
                group_t0 = group_gevp.create_group('t0_%s'%(t0_init+nt[0]))
                
                eigevals_final_mean = vf.NT_TO_NCFGS(evals_mean)
                evals_fits_rs = np.array(vf.RESHAPING_EIGENVALS_FOR_FITS(np.array(evalues_rs), the_size_matrix), dtype=np.float64)
                
                group_eigvecs = group_t0.create_group('Eigenvectors')
                group_eigvecs.create_dataset('Resampled', data=evectors_rs)
                
                l=0; sigma_2 = []
                for l in range(the_size_matrix):
                    dis_eign = vf.NCFGS_TO_NT(evals_fits_rs[l])
                    evals_fits_rs_mean = vf.MEAN(dis_eign)
                    sigma_2.append(vf.COV_MATRIX(dis_eign, evals_fits_rs_mean, the_type_rs))
                evecs_mean = np.array(evecs_mean)
                group_eigvecs.create_dataset('Mean', data=evecs_mean)
                group_eigns = group_t0.create_group('Eigenvalues')
                group_eigns.create_dataset('Resampled', data = evals_fits_rs)
                group_eigns.create_dataset('Mean', data = eigevals_final_mean)
                group_eigns.create_dataset('Covariance_matrix', data = np.array(sigma_2))
                print('T0 = %s'%str(t0_init + nt[0]) + '... DONE')
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
