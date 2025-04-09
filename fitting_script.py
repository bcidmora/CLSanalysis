import numpy as np
import h5py
import time
from iminuit import Minuit
import sys
import os
import set_of_functions as vf

import warnings
warnings.filterwarnings('ignore')


def FitSingleCorrelators(the_data, the_fit_data, the_type_rs, the_list_tmaxs, **kwargs): 
    only_one_tmin = kwargs.get('one_tmin')
    the_type_fit = kwargs.get('type_fit')
    type_correlated_fit = kwargs.get('type_correlation')
    the_names_irreps = list(the_data.keys())
    
    if the_type_fit=='1':
        dof = np.zeros((1,2))
        da_minimization = vf.SINGLE_EXPONENTIAL
        fit_params = ('a0', 'e0')
    elif the_type_fit=='2':
        dof = np.zeros((1,4))
        da_minimization = vf.DOUBLE_EXPONENTIAL
        fit_params = ('a0', 'e0', 'a1','e1')

    for j in range(len(the_names_irreps)):
        the_op_list = list(the_data[the_names_irreps[j]+'/Operators'])
        print('----------------------------------------------------------------------------------------')
        print('IRREP: ' + str(the_names_irreps[j]) + '\n   --->>   Operators list: ')
        for item in the_op_list: print('           ' + str(item.decode("utf-8")))    
        print('----------------------------------------------------------------------------------------')
        
        if the_names_irreps[j] not in the_fit_data.keys(): dis_irrep = the_fit_data.create_group(the_names_irreps[j])
        else: dis_irrep = the_fit_data[the_names_irreps[j]]
        
        if '%sexp'%the_type_fit not in dis_irrep.keys(): 
            one_exp_fit = dis_irrep.create_group('%sexp'%the_type_fit)
            tmin_data = one_exp_fit.create_group('Tmin')
        else:
            one_exp_fit = the_fit_data[the_names_irreps[j]+'/%sexp'%the_type_fit]
            tmin_data = one_exp_fit['Tmin']
        
        corr_info = the_data[the_names_irreps[j]].get('Correlators')
        data_corr_fit = np.array(corr_info.get('Real/Mean'), dtype=np.float64)
        data_corr_fit_rs_raw = np.array(corr_info.get('Real/Resampled'), dtype=np.float64)
        data_corr_fit_rs = vf.NT_TO_NCFGS(data_corr_fit_rs_raw)

        data_cov_matrix = np.array(corr_info.get('Real/Covariance_matrix'), dtype=np.float64)
        eff_energy_hint = np.array(the_data[the_names_irreps[j]].get('Effective_masses/Mean'))
        nt = np.array(the_data[the_names_irreps[j]].get('Time_slices'))  
        
        nt_mod = np.arange(0,len(nt))
        if only_one_tmin: 
            shift_nt = nt[0]-nt_mod[0] #Shift of the real nt versus the nt starting from 0
            ul = int(the_list_tmaxs[j]) - shift_nt # This is the max nt given by us, shifted.
            ll = [nt_mod[5]] # lower limit (min nt for the fit)
        else: 
            shift_nt = nt[0]-nt_mod[0]
            ul = int(the_list_tmaxs[j]) - shift_nt
            ll = np.arange(nt_mod[1],int(ul*0.7)) #it skips the very first one, because in many cases, the correlator here is zero.
        
        if type_correlated_fit=='Correlated':
            if 'Correlated' in tmin_data.keys(): del tmin_data['Correlated']
            fit_data = tmin_data.create_group('Correlated')
            cov_matrix_fit = data_cov_matrix
        elif type_correlated_fit=='Uncorrelated':
            if 'Uncorrelated' in tmin_data.keys(): del tmin_data['Uncorrelated']
            fit_data = tmin_data.create_group('Uncorrelated') 
            cov_matrix_fit = np.zeros((len(data_cov_matrix), len(data_cov_matrix)))
            np.fill_diagonal(cov_matrix_fit, np.diag(data_cov_matrix))

        energies_list, sigmas_list, chi_vals_list, sigmas_chi_list  = [], [], [], []
        another_list = []
        
        begin_time = time.time()
        for yy in ll:
            print('Tmin = ' + str(yy+shift_nt) + '|| TMax = %s'%(ul+shift_nt))
            
            another_useful_list = []
            
            da_hint = vf.BEST_GUESS(data_corr_fit[yy:ul], nt[yy:ul], the_type_fit)
            if False in np.isnan(da_hint):
                dof = da_hint
            else: 
                dof = np.zeros((1,len(da_hint)))
                dof = dof[0]
                dof[0], dof[1] = np.float64(0.1), np.float64(eff_energy_hint[yy])
                
            small_cov = vf.SHRINK_MATRIX(cov_matrix_fit, yy, ul)   
            
            sigma_matrix = np.array(small_cov, dtype=np.float64)
            inverse_cov_m = np.linalg.inv(sigma_matrix) 
            
            the_fit_choice = vf.My_Fits(da_minimization, nt[yy:ul], data_corr_fit[yy:ul], inverse_cov_m, dof, np.float64(0.))            
            the_fit = Minuit(the_fit_choice, dof, name=fit_params)
            
            the_fit.errordef, the_fit.tol = 1e-8, 1e-10
            the_fit.scan()
            the_fit.migrad(iterate=10,ncall=5000)
            
            e0 = np.float64(the_fit.values['e0'])
            
            dof_rs = dof
            dof_rs[0], dof_rs[1] = np.float64(the_fit.values['a0']), e0
            
            energies_list.append(e0); chi_vals_list.append(the_fit.fval); another_useful_list.append(e0)
            
            chi_vals_rs_list = []
            zz=0; 
            for zz in range(len(data_corr_fit_rs)):
                the_fit_choice_rs = vf.My_Fits(da_minimization, nt[yy:ul], data_corr_fit_rs[zz][yy:ul], inverse_cov_m, dof, np.float64(0.))
                the_fit_rs = Minuit(the_fit_choice_rs, dof, name=fit_params)
                
                the_fit_rs.errordef, the_fit_rs.tol = 1e-8, 1e-7
                the_fit_rs.scan()
                the_fit_rs.migrad(iterate=10, ncall=5000)
                
                e0_rs = np.float64(the_fit_rs.values['e0']); chi_vals_rs_list.append(the_fit_rs.fval); another_useful_list.append(e0_rs)

            sigma_fit_rs = vf.STD_DEV(another_useful_list[1:], np.mean(another_useful_list[1:]), the_type_rs)
            sigma_chi_rs = vf.STD_DEV(chi_vals_rs_list, np.mean(chi_vals_rs_list), the_type_rs)
            
            sigmas_list.append(sigma_fit_rs); sigmas_chi_list.append(sigma_chi_rs); another_list.append(np.array(another_useful_list))
          
        fit_data.create_dataset('Resampled', data = np.array(another_list))
        fit_data.create_dataset('Mean', data = np.array([ll + shift_nt, [ul + shift_nt]*(len(ll)), energies_list, sigmas_list, chi_vals_list, sigmas_chi_list]))
        
        print('Minimization %s exp: DONE!'%the_type_fit)
        print('Irrep nr. ' + str(j+1) + ' out of '+str(len(the_data)))
        end_time = time.time()        
        print('Time taken: ' + str(round((end_time-begin_time)/60,2))+' min')
    
    
def FitMultiCorrelators(the_data, the_fit_data, the_type_rs, the_list_tmaxs, **kwargs):
    only_one_tmin = kwargs.get('one_tmin')
    only_one_t0 = kwargs.get('one_t0')
    the_type_fit = kwargs.get('type_fit')
    type_correlated_fit = kwargs.get('type_correlation')
    the_names_irreps = list(the_data.keys())
    
    if only_one_t0:
        t0_s = [int(kwargs.get('chosen_t0'))]
    else:
        t0_s = sorted([int(item[3:]) for item in list(the_data[the_names_irreps[0]+'/GEVP'])])
    
    if the_type_fit=='1':
        dof = np.zeros((1,2))
        da_minimization = vf.SINGLE_EXPONENTIAL
        fit_params = ('a0', 'e0')
    elif the_type_fit=='2':
        dof = np.zeros((1,4))
        da_minimization = vf.DOUBLE_EXPONENTIAL
        fit_params = ('a0', 'e0', 'a1','e1')

    for j in range(len(the_data)):                
        the_op_list = list(the_data[the_names_irreps[j]+'/Operators'])
        
        print('----------------------------------------------------------------------------------------')
        print('IRREP: '+ str(the_names_irreps[j]) + '\n   --->>   Operators list: ')
        for item in the_op_list: print('           ' + str(item.decode("utf-8")))
        print('----------------------------------------------------------------------------------------')
        
        if the_names_irreps[j] not in the_fit_data.keys(): dis_irrep = the_fit_data.create_group(the_names_irreps[j])
        else: dis_irrep = the_fit_data[the_names_irreps[j]]
        
        if 'Operators' not in dis_irrep.keys():
            dis_irrep.create_dataset('Operators', data=the_op_list)
        if kwargs.get('ratio_on')!=None:
            if 'Single_hadron_corrs' in dis_irrep.keys(): del the_fit_data[the_names_irreps[j]+'/Single_hadron_corrs']
            dis_irrep.create_dataset('Single_hadron_corrs', data= list(the_data[the_names_irreps[j]+'/Single_hadron_corrs']))
        
        if '%sexp'%the_type_fit not in dis_irrep.keys(): 
            one_exp_fit = dis_irrep.create_group('%sexp'%the_type_fit)
        else:
            one_exp_fit = the_fit_data[the_names_irreps[j]+'/%sexp'%the_type_fit]
            
        begin_time_tmin = time.time()
        for t0 in t0_s:
            print('T0 equal: ', t0)
            if 't0_%s'%t0 not in one_exp_fit.keys():   
                t0_group = one_exp_fit.create_group('t0_%s'%t0)
                tmin_data = t0_group.create_group('Tmin')
            else:
                t0_group = one_exp_fit['t0_%s'%t0]
                tmin_data = t0_group['Tmin']
            
            corr_info = the_data[the_names_irreps[j]].get('GEVP/t0_%s'%t0)
            data_corr_fit = np.array(corr_info.get('Eigenvalues/Mean')).real 
            data_corr_fit_rs = np.array(corr_info.get('Eigenvalues/Resampled')).real 
            
            data_cov_matrix = np.array(corr_info.get('Eigenvalues/Covariance_matrix')).real 
            eff_energy_hint = np.array(corr_info.get('Effective_masses/Mean'))
            nt = np.array(the_data[the_names_irreps[j]].get('Time_slices'))

            
            if type_correlated_fit=='Correlated':
                if 'Correlated' in tmin_data.keys(): del tmin_data['Correlated']
                fit_data = tmin_data.create_group('Correlated')
            elif type_correlated_fit=='Uncorrelated':
                if 'Uncorrelated' in tmin_data.keys(): del tmin_data['Uncorrelated']
                fit_data = tmin_data.create_group('Uncorrelated') 
            mean_data = fit_data.create_group('Mean')
            rs_data = fit_data.create_group('Resampled')
            
            for ls in range(len(data_corr_fit)):
                nt_mod = np.arange(0,len(nt))
                if only_one_tmin: 
                    shift_nt = nt[0]-nt_mod[0]
                    ul = [x-shift_nt for x in the_list_tmaxs[j]]
                    ll = [nt_mod[5]]                    
                else: 
                    shift_nt = nt[0]-nt_mod[0]
                    ul = [x-shift_nt for x in the_list_tmaxs[j]]
                    ll = np.arange(nt_mod[0]+1, int(ul[ls]*(.8 - (0.025*ls))))
                    
                if type_correlated_fit=='Correlated':
                    cov_matrix_fit = data_cov_matrix[ls]
                elif type_correlated_fit=='Uncorrelated':
                    cov_matrix_fit = np.zeros((len(data_cov_matrix[ls]), len(data_cov_matrix[ls])))
                    np.fill_diagonal(cov_matrix_fit, np.diag(data_cov_matrix[ls]))  
                
                energies_list, sigmas_list, chi_vals_list, sigmas_chi_list = [], [], [], []
                another_list = []
                for yy in ll:
                    print('Tmin = ' + str(yy+shift_nt) + '|| TMax = %s'%(ul[ls]+shift_nt))
                    another_useful_list = []
                    da_hint = vf.BEST_GUESS(data_corr_fit[ls][yy:ul[ls]], nt[yy:ul[ls]], the_type_fit) 
                    if False in np.isnan(da_hint):
                        dof = da_hint
                    else: 
                        dof = np.zeros((1,len(da_hint)));
                        dof = dof[0]
                        dof[0] = np.float64(0.1)
                        dof[1] = np.float64(eff_energy_hint[ls][yy])
                
                    small_cov = vf.SHRINK_MATRIX(cov_matrix_fit, yy, ul[ls])
                    
                    sigma_matrix = np.array(small_cov, dtype=np.float64)
                    inverse_cov_m = np.linalg.inv(sigma_matrix)
                    
                    the_fit_choice = vf.My_Fits(da_minimization, nt[yy:ul[ls]], data_corr_fit[ls][yy:ul[ls]], inverse_cov_m, dof, np.float64(t0))
                    the_fit = Minuit(the_fit_choice, dof, name=fit_params)
                    
                    the_fit.errordef, the_fit.tol = 1e-8, 1e-10
                    the_fit.scan()
                    the_fit.migrad(iterate=10,ncall=5000)
                    
                    e0 = np.float64(the_fit.values['e0'])
                    
                    dof_rs = dof
                    dof_rs[0], dof_rs[1] = np.float64(the_fit.values['a0']), e0
                    
                    energies_list.append(e0); chi_vals_list.append(np.float64(the_fit.fval)); another_useful_list.append(e0)
                    
                    zz=0
                    chi_vals_rs_list = []
                    for zz in range(data_corr_fit_rs.shape[1]):
                        my_fit_choice_rs = vf.My_Fits(da_minimization, nt[yy:ul[ls]], data_corr_fit_rs[ls][zz][yy:ul[ls]], inverse_cov_m, dof, np.float64(t0))
                        the_fit_rs = Minuit(my_fit_choice_rs, dof, name=fit_params)
                        
                        the_fit_rs.errordef, the_fit_rs.tol = 1e-8, 1e-7
                        the_fit_rs.scan()
                        the_fit_rs.migrad(iterate=10, ncall=5000)
                        
                        e0_rs = np.float64(the_fit_rs.values['e0']); chi_vals_rs_list.append(the_fit_rs.fval); another_useful_list.append(e0_rs)
                
                    sigma_fit_rs = vf.STD_DEV(another_useful_list[1:], np.mean(another_useful_list[1:]), the_type_rs)
                    sigma_chi_rs = vf.STD_DEV(chi_vals_rs_list, np.mean(chi_vals_rs_list), the_type_rs)
                    
                    sigmas_list.append(sigma_fit_rs); sigmas_chi_list.append(sigma_chi_rs); another_list.append(np.array(another_useful_list))
                print('E = %s READY'%ls)    
                
                rs_data.create_dataset('lambda_%s'%ls, data=np.array(another_list))
                mean_data.create_dataset('lambda_%s'%ls, data =np.array([ll + shift_nt, [ul[ls]+shift_nt]*len(ll), energies_list, sigmas_list, chi_vals_list, sigmas_chi_list]))              
        end_time_tmin=time.time()
        print('Time taken: ' + str(round((end_time_tmin-begin_time_tmin)/60,2))+' min')
        print('Minimization %s exp: E vs Tmin DONE!'%the_type_fit)
        j+=1
        print('Irrep nr. ' + str(j) + ' out of '+str(len(the_data)))
        
        

if __name__=="__main__":
    
    myEns = str(sys.argv[1]).upper()
    myWhichCorrelator = str(sys.argv[2]).lower()[0]
    myTypeRs = str(sys.argv[3]).lower()
    myRebinOn = str(sys.argv[4])
    
    myRb = 2
    myVersion = 'test'
    
    myTypeFit = '1'
    myTypeCorrelation = 'Correlated'
    myOneTMin = True 
    myOneT0 =  True
    myT0 = 4
    
    if myEns == 'N451': from files_n451 import listTMaxSingleHads, listTMaxMultiHads
    elif myEns == 'N201': from files_n201 import listTMaxSingleHads, listTMaxMultiHads 
    elif myEns == 'D200': from files_d200 import listTMaxSingleHads, listTMaxMultiHads
    elif myEns == 'X451': from files_x451 import listTMaxSingleHads, listTMaxMultiHads
    
    myLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'$YOUR_OUTPUT_PATH(SAME_THAN_CORRS_SCRIPT_OUTPUT)$/%s/'%myEns)
    
    if myRebinOn=='rb':
        reBin ='_bin'+str(myRb)
    else:
        reBin=''  
    
    vf.INFO_PRINTING(myWhichCorrelator, myEns)
    
    if myWhichCorrelator=='s':
        myData = h5py.File(myLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion, 'r') 
        
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_SingleHadrons/')
        myFitData = h5py.File(myFitsLocation + 'Single_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        FitSingleCorrelators(myData, myFitData, myTypeRs, listTMaxSingleHads, one_tmin=myOneTMin, type_fit=myTypeFit, type_correlation=myTypeCorrelation)
        
    elif myWhichCorrelator=='m':
        myData = h5py.File(myLocation + '/Matrix_correlators_' +  myTypeRs + reBin + '_v%s.h5'%myVersion,'r')
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_Matrices/')
        myFitData = h5py.File(myFitsLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        FitMultiCorrelators(myData, myFitData, myTypeRs, listTMaxMultiHads, type_fit = myTypeFit, type_correlation = myTypeCorrelation, one_tmin = myOneTMin, one_t0 = myOneT0, chosen_t0 = myT0)
     
    elif myWhichCorrelator=='mr':
        myData = h5py.File(myLocation + '/Matrix_correlators_ratios_' + ratioStr + myTypeRs + reBin + '_v%s.h5'%myVersion,'r')
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_Ratios/')
        myFitData = h5py.File(myFitsLocation + 'Matrix_correlators_ratios_' + ratioStr + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        FitMultiCorrelators(myData, myFitData, myTypeRs, listTMaxMultiHads, type_fit = myTypeFit, type_correlation = myTypeCorrelation, one_tmin = myOneTMin, one_t0 = myOneT0, chosen_t0 = myT0, ratio_on='yes')

    myFitData.close()
    myData.close()
