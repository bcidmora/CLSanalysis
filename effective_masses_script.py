import numpy as np
import h5py
import time
import sys
import os
import set_of_functions as vf
import warnings
warnings.filterwarnings('ignore')


def SingleCorrelatorEffectiveMass(the_single_correlator_data, the_type_rs,**kwargs):   
    
    ### Defining distance between time-slice elements of the correlator
    if kwargs.get('dist_eff_mass')!=None and kwargs.get('dist_eff_mass')!=1:
        the_dist_eff_mass = int(kwargs.get('dist_eff_mass'))
    else:
        the_dist_eff_mass = 1 # Default is 1
    
    ### The irreps
    the_list_name_irreps = list(the_single_correlator_data.keys())
    
    begin_time = time.time()
    for j in range(len(the_list_name_irreps)):
        
        ### Extracting data from file
        this_data = the_single_correlator_data[the_list_name_irreps[j]]
                
        ### Getting the operators list and the time interval
        the_op_list, nt = list(this_data.get('Operators')), np.array(this_data.get('Time_slices'))
        
        ### The real part of the resampled data
        the_rs_real = np.array(this_data.get('Correlators/Real/Resampled'))
        the_rs_real = vf.NT_TO_NCFGS(the_rs_real)
        
        ### If the Effective Masses were already computed, this part gets deleted and created a new branch.
        if 'Effective_masses' in this_data.keys(): del the_single_correlator_data[the_list_name_irreps[j]+'/Effective_masses']
        
        group_em = this_data.create_group('Effective_masses')
        
        ### Mean values of the real part of the correlator to get the effective masses
        the_mean_corr_real = np.array(this_data.get('Correlators/Real/Mean'))
        
        ### Effective Mass computation
        the_em_rs_f =  vf.EFF_MASS(the_mean_corr_real, the_dist_eff_mass)
        
        ### Loop over the resampled data
        l, the_em_rs = 0, []
        for l in range(len(the_rs_real)):
            the_em_rs.append(vf.EFF_MASS(the_rs_real.real[l],the_dist_eff_mass))
        
        ### Reshaping data
        the_em_rs = vf.NCFGS_TO_NT(the_em_rs)
        
        ### Getting the mean value of the resamples
        the_mrs_f_real_rs  = []
        for tt in range(len(the_em_rs)):
            the_mrs_f_real_rs.append(np.mean(the_em_rs[tt]))
        the_mrs_f_real_rs = np.array(the_mrs_f_real_rs)
        
        ### Sigma values for the resamples
        the_sigma_eff_mass = vf.STD_DEV_MEAN(the_em_rs, the_mrs_f_real_rs, the_type_rs)
        
        group_em.create_dataset('Mean', data=the_em_rs_f)
        group_em.create_dataset('Sigmas', data=the_sigma_eff_mass)
        
        print('Irrep nr.: '+ str(j+1) + ' out of ' +str(len(the_list_name_irreps)))
    end_time = time.time()
    print('TOTAL TIME TAKEN: ' + str((end_time-begin_time)/60) +' mins')
            
            
def MultiCorrelatorEffectiveMass(the_matrix_correlator_data, the_type_rs, **kwargs):
    
    ### Defining distance between time-slice elements of the correlator
    if kwargs.get('dist_eff_mass')!=None and kwargs.get('dist_eff_mass')!=1:
        the_dist_eff_mass = int(kwargs.get('dist_eff_mass'))
    else:
        the_dist_eff_mass = 1
        
    ### The irreps
    the_list_name_irreps = list(the_matrix_correlator_data.keys())    
    
        
    begin_time = time.time()
    for j in range(len(the_list_name_irreps)): 
        
        ### Extracting data from file
        this_data = the_matrix_correlator_data[the_list_name_irreps[j]]
        
        ### Getting the operators list and the time interval
        the_op_list, nt = list(this_data.get('Operators')), np.array(this_data.get('Time_slices'))
        
        ### This is the size of the correlator matrices
        the_size_matrix = len(the_op_list)

        ### The real part of the resampled data
        the_rs_real = np.array(this_data.get('Correlators/Real/Resampled'))
        
        ### Mean values of the real part of the correlator to get the effective masses
        the_mean_corr_real = np.array(this_data.get('Correlators/Real/Mean'))
        
        ### Reshaping data
        the_reshaped_mean_corr = vf.RESHAPING_CORRELATORS(the_mean_corr_real)
        the_reshaped_rs_corr = vf.RESHAPING_CORRELATORS_RS(the_rs_real)
        
        ### Loop over the nr. of operators = size of the correlator matrix
        ii, the_efm_mass, the_sigma_efm = 0, [], []
        for ii in range(the_size_matrix):
            
            ### Effective Masses of the entral values of the correlators
            the_mean_eff = vf.EFF_MASS(the_reshaped_mean_corr[ii][ii],the_dist_eff_mass)
            
            the_efm_mass.append(the_mean_eff)
            
            ### Loop over the resamples
            the_rs_eff = []
            for xyz in range(len(the_rs_real[0])):
                the_rs_eff.append(vf.EFF_MASS(the_reshaped_rs_corr[ii][ii][xyz], the_dist_eff_mass))
            
            ### Reshaping data
            the_rs_eff = np.array(vf.NCFGS_TO_NT(the_rs_eff))
            
            ### MEan value of the resampled data to compute the sigma vals.
            the_rs_mean_eff = vf.MEAN(the_rs_eff)
            the_sigma_efm.append(vf.STD_DEV_MEAN(the_rs_eff, the_rs_mean_eff, the_type_rs))
                
            ### If the branch Effective Masses exists in the file, then it gets deleted and created a new one with the new values.
            if 'Effective_masses' in this_data['Correlators/Real'].keys(): del the_matrix_correlator_data[the_list_name_irreps[j]+'/Correlators/Real/Effective_masses']
            this_data.get('Correlators/Real').create_dataset('Effective_masses', data=np.array(the_efm_mass))
            
            if 'Effective_masses_sigmas' in this_data['Correlators/Real'].keys(): del the_matrix_correlator_data[the_list_name_irreps[j]+'/Correlators/Real/Effective_masses_sigmas']
            this_data.get('Correlators/Real').create_dataset('Effective_masses_sigmas', data=np.array(the_sigma_efm))
            
        ### If the GEVP was performed, then it will identify this section.
        if 'GEVP' in this_data.keys(): 
            
            gevp_group = this_data.get('GEVP')
                
            print("Effective Masses of GEVP eigenvalues in process...")
            
            vf.DOING_EFFECTIVE_MASSES_EIGENVALUES(gevp_group, the_dist_eff_mass, the_type_rs)
                
        ### If the Operator Analysis was performed, this section will also be identified. 
        if 'Operators_Analysis' in this_data.keys():
            
            ### Checks if the keys of ops_chosen is in the list of the operator analysis
            if any('Ops_chosen_' in the_keys for the_keys in this_data['Operators_Analysis'].keys()):
                
                the_list_of_chosen_ops = list(filter(lambda x: 'Ops_chosen' in x, this_data['Operators_Analysis'].keys()))
                
                for the_op_item in the_list_of_chosen_ops:
                    gevp_group = this_data['Operators_Analysis/'].get(the_op_item)
                
                    vf.DOING_EFFECTIVE_MASSES_EIGENVALUES(gevp_group, the_dist_eff_mass, the_type_rs)
                    
            if any('Add_Op' in the_keys for the_keys in this_data['Operators_Analysis'].keys()):
                the_list_of_chosen_ops = list(filter(lambda x: "Add_Op" in x, this_data['Operators_Analysis'].keys()))
                
                for the_op_item in the_list_of_chosen_ops:
                    gevp_group = this_data['Operators_Analysis/'].get(the_op_item)
                    
                    vf.DOING_EFFECTIVE_MASSES_EIGENVALUES(gevp_group, the_dist_eff_mass, the_type_rs)
#                     
            if any('Remove_Op' in the_keys for the_keys in this_data['Operators_Analysis'].keys()):
                the_list_of_chosen_ops = list(filter(lambda x: "Remove_Op" in x, this_data['Operators_Analysis'].keys()))
                
                for the_op_item in the_list_of_chosen_ops:
                    gevp_group = this_data['Operators_Analysis/'].get(the_op_item)
                    
                    vf.DOING_EFFECTIVE_MASSES_EIGENVALUES(gevp_group, the_dist_eff_mass, the_type_rs)
            
            print('Irrep nr.: '+ str(j+1) + ' out of ' +str(len(the_list_name_irreps)))
    end_time = time.time()
    print('TIME TAKEN: ' + str((end_time-begin_time)/60) +' mins')



### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------


if __name__=="__main__":
    myEns = str(sys.argv[1]).upper()
    myWhichCorrelator = str(sys.argv[2]).lower()
    myTypeRs = str(sys.argv[3]).lower()
    myRebinOn = str(sys.argv[4]).lower()
    myRb = 1
    myVersion = 'test'
    myEffMassDistance = 1 #None #2 #3
    
    if myRebinOn=='rb': 
        rb = int(myRb)
        reBin = '_bin'+str(rb)
    else:
        reBin = ''  
    
    if myEns == 'N451': from files_n451 import location
    elif myEns == 'N201': from files_n201 import location
    elif myEns == 'D200': from files_d200 import location
    elif myEns == 'X451': from files_x451 import location
    
    myLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'$YOUR_OUTPUT_PATH(TOTAL_SAME_THAN_CORRS_SCRIPT)$/%s/'%myEns)
    
    vf.INFO_PRINTING(myWhichCorrelator, myEns)
    
    if myWhichCorrelator=='s':
        myNameArchivo = myLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion
        mySingleCorrelatorData = h5py.File(myNameArchivo,'r+')            
        SingleCorrelatorEffectiveMass(mySingleCorrelatorData,  myTypeRs, dist_eff_mass = myEffMassDistance) 
        mySingleCorrelatorData.close()
    
    elif myWhichCorrelator=='m':        
        myNameArchivo = myLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion
        myMatrixCorrelatorData = h5py.File(myNameArchivo, 'r+')
        MultiCorrelatorEffectiveMass(myMatrixCorrelatorData, myTypeRs, dist_eff_mass = myEffMassDistance)
        myMatrixCorrelatorData.close()
    
    elif myWhichCorrelator=='mr':
        myNameArchivo = myLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin + '_v%s.h5'%myVersion
        myRatioMatrixCorrelatorData = h5py.File(myNameArchivo, 'r+')
        MultiCorrelatorEffectiveMass(myRatioMatrixCorrelatorData, myTypeRs)
        myRatioMatrixCorrelatorData.close()
    
    print('-'*(len(myNameArchivo)+1))
    print('Saved as: \n' + myNameArchivo)
    print('_'*(len(myNameArchivo)+1))
