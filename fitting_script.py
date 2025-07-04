import numpy as np
import h5py
import time
from iminuit import Minuit
import sys
import os
import set_of_functions as vf

import warnings
warnings.filterwarnings('ignore')


def FitSingleCorrelators(the_data, the_fit_data, the_type_rs, the_list_tmaxs, the_irreps, **kwargs): 
    
    print("                     FITTING \n")
    
    ### If only one tmin needs to be done for the fitting
    the_only_one_tmin = kwargs.get('one_tmin')
    
    ### The type of fit
    the_type_fit = kwargs.get('type_fit')
    
    ### Correlated or uncorrelated
    type_correlated_fit = kwargs.get('type_correlation')
    
    ### The name of the irreps
    the_s_irreps = list(the_data.keys())
    
    ### Here we have three types of fit functions (1-exp, 2-exp and geometric)
    if the_type_fit=='1':
        the_dof = np.zeros((1,2))
        da_minimization = vf.SINGLE_EXPONENTIAL
        the_fit_params = ('a0', 'e0')
    elif the_type_fit=='2':
        the_dof = np.zeros((1,4))
        da_minimization = vf.DOUBLE_EXPONENTIAL
        the_fit_params = ('a0', 'e0', 'a1','e1')
    elif the_type_fit=='g':
        the_dof = np.zeros((1,4))
        da_minimization = vf.GEOMETRIC_FORM
        the_fit_params = ('a0', 'e0', 'b','m')
    
    ### Loop over the irreps
    for the_irrep in the_s_irreps:
        
        ### List of operators of this irrep
        the_op_list = list(the_data[the_irrep+'/Operators'])
        print('----------------------------------------------------------------------------------------')
        print('IRREP (%s/'%str(the_irreps.index(the_irrep)+1) + '%s): '%str(len(the_irreps))  + str(the_irrep) + '\n   --->>   Operators list: ')
        for item in the_op_list: print('           ' + str(item.decode("utf-8")))    
        print('----------------------------------------------------------------------------------------')
        
        ### Check if this part already exists
        if the_irrep not in the_fit_data.keys(): dis_irrep = the_fit_data.create_group(the_irrep)
        else: dis_irrep = the_fit_data[the_irrep]
        
        if 'Operators' not in dis_irrep.keys(): 
            dis_irrep.create_dataset('Operators', data=the_op_list)
        
        if '%sexp'%the_type_fit not in dis_irrep.keys(): 
            the_one_exp_fit = dis_irrep.create_group('%sexp'%the_type_fit)
            the_tmin_data = the_one_exp_fit.create_group('Tmin')
        else:
            the_one_exp_fit = the_fit_data[the_irrep+'/%sexp'%the_type_fit]
            the_tmin_data = the_one_exp_fit['Tmin']
        
        the_corr = the_data[the_irrep].get('Correlators')
        the_corr_fit = np.array(the_corr.get('Real/Mean'), dtype=np.float64)
        the_corr_fit_rs_raw = np.array(the_corr.get('Real/Resampled'), dtype=np.float64)
        the_corr_fit_rs = vf.NT_TO_NCFGS(the_corr_fit_rs_raw)

        the_cov_matrix = np.array(the_corr.get('Real/Covariance_matrix'), dtype=np.float64)
        the_eff_energy_hint = np.array(the_data[the_irrep].get('Effective_masses/Mean'))
        the_nt = np.array(the_data[the_irrep].get('Time_slices'))  
        
        nt_mod = np.arange(0,len(the_nt))
        if the_only_one_tmin: 
            the_ul = int(the_list_tmaxs[the_irreps.index(the_irrep)]) - the_nt[0] # This is the max nt given by us, shifted.
            the_ll = [nt_mod[5]] # lower limit (min nt for the fit)
        else: 
            the_ul = int(the_list_tmaxs[the_irreps.index(the_irrep)]) - the_nt[0]
            the_ll = np.arange(nt_mod[1],int(the_ul*0.7)) #it skips the very first one, because in many cases, the correlator here is zero.
        
        if type_correlated_fit=='Correlated':
            if 'Correlated' in the_tmin_data.keys(): del the_tmin_data['Correlated']
            the_fit_data = the_tmin_data.create_group('Correlated')
            the_cov_matrix_fit = the_cov_matrix
        elif type_correlated_fit=='Uncorrelated':
            if 'Uncorrelated' in the_tmin_data.keys(): del the_tmin_data['Uncorrelated']
            the_fit_data = the_tmin_data.create_group('Uncorrelated') 
            the_cov_matrix_fit = np.zeros((len(the_cov_matrix), len(the_cov_matrix)))
            np.fill_diagonal(the_cov_matrix_fit, np.diag(the_cov_matrix))

        the_energies_list, the_sigmas_list, the_chi_vals_list, the_sigmas_chi_list  = [], [], [], []
        another_list = []
        
        begin_time = time.time()
        for the_yy in the_ll:
            print('Tmin = ' + str(the_yy+the_nt[0]) + '|| TMax = %s'%(the_ul+the_nt[0]))
            
            another_useful_list = []
            
            da_hint = vf.BEST_GUESS(the_corr_fit[the_yy:the_ul], the_nt[the_yy:the_ul], the_type_fit)
            if False in np.isnan(da_hint):
                the_dof = da_hint
            else: 
                the_dof = np.zeros((1,len(da_hint)))
                the_dof = the_dof[0]
                the_dof[0], the_dof[1] = np.float64(0.1), np.float64(the_eff_energy_hint[the_yy])
                
            the_small_cov = vf.SHRINK_MATRIX(the_cov_matrix_fit, the_yy, the_ul)   
            
            the_sigma_matrix = np.array(the_small_cov, dtype=np.float64)
            the_inverse_cov_m = np.linalg.inv(the_sigma_matrix) 
            
            the_fit_choice = vf.My_Fits(da_minimization, the_nt[the_yy:the_ul], the_corr_fit[the_yy:the_ul], the_inverse_cov_m, the_dof, np.float64(0.))            
            the_fit = Minuit(the_fit_choice, the_dof, name=the_fit_params)
            
            the_fit.errordef, the_fit.tol = 1e-8, 1e-10
            the_fit.scan()
            the_fit.migrad(iterate=10,ncall=5000)
            
            e0 = np.float128(the_fit.values['e0'])
            
            the_dof_rs = the_dof
            the_dof_rs[0], the_dof_rs[1] = np.float128(the_fit.values['a0']), e0
            
            the_energies_list.append(e0); the_chi_vals_list.append(the_fit.fval); another_useful_list.append(e0)
            
            chi_vals_rs_list = []
            zz=0; 
            for zz in range(len(the_corr_fit_rs)):
                the_fit_choice_rs = vf.My_Fits(da_minimization, the_nt[the_yy:the_ul], the_corr_fit_rs[zz][the_yy:the_ul], the_inverse_cov_m, the_dof_rs, np.float128(0.))
                the_fit_rs = Minuit(the_fit_choice_rs, the_dof, name=the_fit_params)
                
                the_fit_rs.errordef, the_fit_rs.tol = 1e-8, 1e-7
                the_fit_rs.scan()
                the_fit_rs.migrad(iterate=10, ncall=5000)
                
                e0_rs = np.float128(the_fit_rs.values['e0']); chi_vals_rs_list.append(the_fit_rs.fval); another_useful_list.append(e0_rs)

            the_sigma_fit_rs = vf.STD_DEV(another_useful_list[1:], np.mean(another_useful_list[1:]), the_type_rs)
            the_sigma_chi_rs = vf.STD_DEV(chi_vals_rs_list, np.mean(chi_vals_rs_list), the_type_rs)
            
            the_sigmas_list.append(the_sigma_fit_rs); the_sigmas_chi_list.append(the_sigma_chi_rs); another_list.append(np.array(another_useful_list))
          
        the_fit_data.create_dataset('Resampled', data = np.array(another_list))
        the_fit_data.create_dataset('Mean', data = np.array([the_ll + the_nt[0], [the_ul + the_nt[0]]*(len(the_ll)), the_energies_list, the_sigmas_list, the_chi_vals_list, the_sigmas_chi_list]))
        
        print('Minimization %s exp: DONE!'%the_type_fit)
        end_time = time.time()        
        print('Time taken: ' + str(round((end_time-begin_time)/60,2))+' min')
    
    
def FitMultiCorrelators(the_data, the_fit_data, the_type_rs, the_list_tmaxs, the_irreps, **kwargs):
    
    print("                     FITTING \n")
    
    ### Fits only one minimum time slice
    the_only_one_tmin = kwargs.get('one_tmin')
    
    ### Fits only one t0 chosen
    the_only_one_t0 = kwargs.get('one_t0')
    
    ### This is the type of fit
    the_type_fit = kwargs.get('type_fit')
    
    ### Correlated or uncorrelated fit
    the_type_correlated_fit = kwargs.get('type_correlation')
    
    ### The names of the irreps in this ensemble
    the_m_irreps = list(the_data.keys())
    
    ### If only one t0 was chosen, here it takes it
    if the_only_one_t0:
        the_t0_s = [int(kwargs.get('chosen_t0'))]
    else:
        the_t0_s = sorted([int(item[3:]) for item in list(the_data[the_m_irreps[0]+'/GEVP'])])
    
    ### Type of fit for the procedure. 
    ### 1-exponential fit
    if the_type_fit=='1':
        the_dof = np.zeros((1,2))
        da_minimization = vf.SINGLE_EXPONENTIAL
        the_fit_params = ('a0', 'e0')
    ### 2-exponential fit
    elif the_type_fit=='2':
        the_dof = np.zeros((1,4))
        da_minimization = vf.DOUBLE_EXPONENTIAL
        the_fit_params = ('a0', 'e0', 'a1','e1')
    elif the_type_fit=='g':
        the_dof = np.zeros((1,4))
        da_minimization = vf.GEOMETRIC_FORM
        the_fit_params = ('a0', 'e0', 'b','m')
    
    ### Loop over the irreps
    for the_irrep in the_m_irreps:           
         
        ### Getting the time slices of this dataset
        the_nt = np.array(the_data[the_irrep].get('Time_slices'))
        
        ### These are the operators used for the full matrix
        the_op_list = list(the_data[the_irrep+'/Operators'])
        
        print('----------------------------------------------------------------------------------------')
        print('IRREP (%s/'%str(the_irreps.index(the_irrep)+1) + '%s): '%str(len(the_irreps)) + str(the_irrep) + '\n   --->>   Operators list: ')
        for item in the_op_list: print('           ' + str(item.decode("utf-8")))
        print('----------------------------------------------------------------------------------------')
        
        ### This is creating a path for the irrep
        if the_irrep not in the_fit_data.keys(): dis_irrep = the_fit_data.create_group(the_irrep)
        else: dis_irrep = the_fit_data[the_irrep]
        
        ### It creates a path to store the operators
        if 'Operators' not in dis_irrep.keys():
            dis_irrep.create_dataset('Operators', data=the_op_list)
        
        ### This is in case of the ratio of correlators are chosen
        if kwargs.get('ratio_on')!=None:
            if 'Single_hadron_corrs' in dis_irrep.keys(): del the_fit_data[the_irrep+'/Single_hadron_corrs']
            dis_irrep.create_dataset('Single_hadron_corrs', data= list(the_data[the_irrep+'/Single_hadron_corrs']))
        
        ### Searching if the path for the data exists
        if '%sexp'%the_type_fit not in dis_irrep.keys(): 
            the_exp_fit = dis_irrep.create_group('%sexp'%the_type_fit)
        else:
            the_exp_fit = the_fit_data[the_irrep+'/%sexp'%the_type_fit]
            
            
        ### Checking the GEVP was done before
        if 'GEVP' in list(the_data[the_irrep].keys()) and kwargs.get('gevp')==True:            
            
            begin_time_tmin = time.time()
            for the_t0 in the_t0_s:
                print('T0 equal: ', the_t0)
                
                ### Checking if the path inside this file exists
                if 't0_%s'%the_t0 not in the_exp_fit.keys():   
                    t0_group = the_exp_fit.create_group('t0_%s'%the_t0)
                    tmin_data = t0_group.create_group('Tmin')
                else:
                    t0_group = the_exp_fit['t0_%s'%the_t0]
                    tmin_data = t0_group['Tmin']
                
                ### Retrieving data for the fits (eigenvalues and resamples)
                the_corr = the_data[the_irrep].get('GEVP/t0_%s'%the_t0)
                vf.DOING_THE_FITTING(the_corr, the_nt, the_type_rs, the_irreps, the_irrep, tmin_data, the_type_correlated_fit, the_type_fit, the_only_one_tmin, the_t0, the_list_tmaxs, da_minimization, the_fit_params)
            end_time_tmin=time.time()
            print('Time taken: ' + str(round((end_time_tmin-begin_time_tmin)/60,2))+' min')
            print('Minimization %s exp: E vs Tmin DONE!'%the_type_fit)
            
        ### Checking the Operators Analysis was done before
        if 'Operators_Analysis' in list(the_data[the_irrep].keys()) and kwargs.get('operators_analysis')==True:    
            
            if kwargs.get('the_operator_analysis_method')=='from_list': the_method = 'Ops_chosen_'
            
            elif kwargs.get('the_operator_analysis_method')=='adding': the_method = 'Add_Op_'
                
            elif kwargs.get('the_operator_analysis_method')=='removing': the_method = 'Remove_Op_'
            
            else: sys.exit("No method for the operators analysis chosen")
                
            the_list_of_chosen_ops = list(filter(lambda x: the_method in x, the_data[the_irrep+'/Operators_Analysis'].keys()))
            
                ### Loop over those elements
            for the_op_item in the_list_of_chosen_ops:
                
                if the_only_one_t0:
                    the_t0_s = [int(kwargs.get('chosen_t0'))]
                else:
                    the_t0_s = sorted([int(item[3:]) for item in list(the_data[the_irrep+'/Operators_Analysis/'+the_op_item])])
                
                ### Searching if the path for the data exists
                if the_op_item+'_%sexp'%the_type_fit not in dis_irrep.keys(): 
                    the_exp_fit = dis_irrep.create_group(the_op_item+'_%sexp'%the_type_fit)
                else:
                    the_exp_fit = the_fit_data[the_irrep+'/'+the_op_item+'_%sexp'%the_type_fit]
                
                begin_time_tmin = time.time()
                for the_t0 in the_t0_s:
                    print('Modified Eigenvalues -->> T0 equal: ', the_t0)
                    
                    ### Checking if the path inside this file exists
                    if 't0_%s'%the_t0 not in the_exp_fit.keys():   
                        t0_group = the_exp_fit.create_group('t0_%s'%the_t0)
                        tmin_data = t0_group.create_group('Tmin')
                    else:
                        t0_group = the_exp_fit['t0_%s'%the_t0]
                        tmin_data = t0_group['Tmin']
                    
                    ### Retrieving data for the fits (eigenvalues and resamples)
                    the_corr = the_data[the_irrep].get('Operators_Analysis/'+the_op_item+'/t0_%s'%the_t0)
                    vf.DOING_THE_FITTING(the_corr, the_nt, the_type_rs, the_irreps, the_irrep, tmin_data, the_type_correlated_fit, the_type_fit, the_only_one_tmin, the_t0, the_list_tmaxs, da_minimization, the_fit_params)
                end_time_tmin=time.time()
                print('Time taken: ' + str(round((end_time_tmin-begin_time_tmin)/60,2))+' min')
                print('Minimization %s exp: E vs Tmin DONE!'%the_type_fit)  
    
        

if __name__=="__main__":
    
    ### This is the ensemble you are analysing
    myEns = str(sys.argv[1]).upper()
    
    ### Single hadrons 's' or multihadrons 'm' correlators
    myWhichCorrelator = str(sys.argv[2]).lower()[0]
    
    ### type of resampling done before
    myTypeRs = str(sys.argv[3]).lower()
    
    ### Rebinning
    myRebinOn = str(sys.argv[4])
    
    myRb = 1
    myVersion = '_test'
    
    ### Type of fit, it could be 1-exp '1', 2-exp '2' or geometric 'g'
    myTypeFit = '1' # '2' # 'g'
    
    ### Correlated or uncorrelated fit
    myTypeCorrelation = 'Correlated' # 'Uncorrelated'
    
    ### One can choose only 1 tmin to do or all of them in a certain range.
    myOneTMin = True 
    
    ### ALso for the gevp results, one can do all t0s or just one
    myOneT0 =  True
    myT0 = 4
    
    ### Info for the fits form the ensembles files
    if myEns == 'N451': from files_n451 import listTMaxSingleHads, listTMaxMultiHads, name, name1
    elif myEns == 'N201': from files_n201 import listTMaxSingleHads, listTMaxMultiHads, name, name1 
    elif myEns == 'D200': from files_d200 import listTMaxSingleHads, listTMaxMultiHads, name, name1
    elif myEns == 'X451': from files_x451 import listTMaxSingleHads, listTMaxMultiHads, name, name1
    
    ### Root directory where the averaged correlators are stored
    myLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'$YOUR_OUTPUT_PATH(SAME_THAN_CORRS_SCRIPT_OUTPUT)$/%s/'%myEns)
    
    if myRebinOn=='rb':
        reBin ='_bin'+str(myRb)
    else:
        reBin=''  
    
    ### Just printing some info
    vf.INFO_PRINTING(myWhichCorrelator, myEns)
    
    
    if myWhichCorrelator=='s':
        myIrreps = name1
        
        myData = h5py.File(myLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion, 'r') 
        
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_SingleHadrons/')
        
        myFitData = h5py.File(myFitsLocation + 'Single_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        FitSingleCorrelators(myData, myFitData, myTypeRs, listTMaxSingleHads, myIrreps, one_tmin = myOneTMin, type_fit = myTypeFit, type_correlation = myTypeCorrelation)

        
    elif myWhichCorrelator=='m':
        myIrreps = name
        
        myOperatorAnalysisMethod = 'from_list' # 'adding' # 'removing' # 'from_list'
        
        myData = h5py.File(myLocation + '/Matrix_correlators_' +  myTypeRs + reBin + '_v%s.h5'%myVersion,'r')
        
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_Matrices/')
        
        myFitData = h5py.File(myFitsLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        FitMultiCorrelators(myData, myFitData, myTypeRs, listTMaxMultiHads, myIrreps,  type_fit = myTypeFit, type_correlation = myTypeCorrelation, one_tmin = myOneTMin, one_t0 = myOneT0, chosen_t0 = myT0, gevp=True, operators_analysis = False, the_operator_analysis_method = myOperatorAnalysisMethod)
        
    ### This ratio of correlators has not been fixed yet. Probably it doesnt work
    elif myWhichCorrelator=='mr':
        myData = h5py.File(myLocation + '/Matrix_correlators_ratios_' + ratioStr + myTypeRs + reBin + '_v%s.h5'%myVersion,'r')
        myFitsLocation = vf.DIRECTORY_EXISTS(myLocation + 'Fits_Ratios/')
        myFitData = h5py.File(myFitsLocation + 'Matrix_correlators_ratios_' + ratioStr + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        FitMultiCorrelators(myData, myFitData, myTypeRs, listTMaxMultiHads, type_fit = myTypeFit, type_correlation = myTypeCorrelation, one_tmin = myOneTMin, one_t0 = myOneT0, chosen_t0 = myT0, ratio_on='yes')

    myFitData.close()
    myData.close()
