import numpy as np
import h5py
import time
import sys
import os
import set_of_functions as vf
import warnings
warnings.filterwarnings('ignore')


def SingleCorrelatorEffectiveMass(the_single_correlator_data, the_type_rs,**kwargs):  
    
    if kwargs.get('dist_eff_mass')!=None and kwargs.get('dist_eff_mass')!=1:
        the_dist_eff_mass = int(kwargs.get('dist_eff_mass'))
    else:
        the_dist_eff_mass = 1 
        
    the_list_name_irreps = list(the_single_correlator_data.keys())
    
    begin_time = time.time()
    for j in range(len(the_list_name_irreps)):
        this_data = the_single_correlator_data[the_list_name_irreps[j]]
                
        the_op_list, nt = list(this_data.get('Operators')), np.array(this_data.get('Time_slices'))
        
        the_size_matrix = len(the_op_list)
        
        rs_real = np.array(this_data.get('Correlators/Real/Resampled'))
        rs_real = vf.NT_TO_NCFGS(rs_real)
        
        if 'Effective_masses' in this_data.keys(): del the_single_correlator_data[the_list_name_irreps[j]+'/Effective_masses']
        
        group_em = this_data.create_group('Effective_masses')
        mrs_f_real = np.array(this_data.get('Correlators/Real/Mean'))
        
        em_rs_f =  vf.EFF_MASS(mrs_f_real,the_dist_eff_mass)
        
        l=0; em_rs=[]
        for l in range(len(rs_real)):
            em_rs.append(vf.EFF_MASS(rs_real.real[l],the_dist_eff_mass))
        em_rs = vf.NCFGS_TO_NT(em_rs)
        
        mrs_f_real_rs  = []
        for tt in range(len(em_rs)):
            mrs_f_real_rs.append(np.mean(em_rs[tt]))
        mrs_f_real_rs = np.array(mrs_f_real_rs )
        
        sigma_eff_mass = vf.STD_DEV_MEAN(em_rs, mrs_f_real_rs, the_type_rs)
        
        group_em.create_dataset('Mean', data=em_rs_f)
        group_em.create_dataset('Sigmas', data=sigma_eff_mass)
        
        print('Irrep nr.: '+ str(j+1) + ' out of ' +str(len(the_list_name_irreps)))
    end_time = time.time()
    print('TOTAL TIME TAKEN: ' + str((end_time-begin_time)/60) +' mins')
            
            
def MultiCorrelatorEffectiveMass(the_matrix_correlator_data, the_type_rs, **kwargs):
    
    the_list_name_irreps = list(the_matrix_correlator_data.keys())    
    
    if kwargs.get('dist_eff_mass')!=None and kwargs.get('dist_eff_mass')!=1:
        the_dist_eff_mass = int(kwargs.get('dist_eff_mass'))
    else:
        the_dist_eff_mass = 1
        
    begin_time = time.time()
    for j in range(len(the_list_name_irreps)): 
        this_data = the_matrix_correlator_data[the_list_name_irreps[j]]
        
        the_op_list, nt = list(this_data.get('Operators')), np.array(this_data.get('Time_slices'))
        the_size_matrix = len(the_op_list)

        rs_real = np.array(this_data.get('Correlators/Real/Resampled'))
        mean_corr_real = np.array(this_data.get('Correlators/Real/Mean'))

        reshaped_mean_corr = vf.RESHAPING_CORRELATORS(mean_corr_real,the_size_matrix)
        reshaped_rs_corr = vf.RESHAPING_CORRELATORS_RS(rs_real,the_size_matrix)
        
        ii=0; efm_mass=[]; sigma_efm=[]
        for ii in range(the_size_matrix):
            mean_eff = vf.EFF_MASS(reshaped_mean_corr[ii][ii],the_dist_eff_mass)
            
            efm_mass.append(mean_eff)
            rs_eff=[]
            for xyz in range(len(rs_real[0])):
                rs_eff.append(vf.EFF_MASS(reshaped_rs_corr[ii][ii][xyz], the_dist_eff_mass))
            
            rs_eff = np.array(vf.NCFGS_TO_NT(rs_eff))
            
            rs_mean_eff = vf.MEAN(rs_eff)
            sigma_efm.append(vf.STD_DEV_MEAN(rs_eff, rs_mean_eff, the_type_rs))
                
            if 'Effective_masses' in this_data['Correlators/Real'].keys(): del the_matrix_correlator_data[the_list_name_irreps[j]+'/Correlators/Real/Effective_masses']
            this_data.get('Correlators/Real').create_dataset('Effective_masses', data=np.array(efm_mass))
            
            if 'Effective_masses_sigmas' in this_data['Correlators/Real'].keys(): del the_matrix_correlator_data[the_list_name_irreps[j]+'/Correlators/Real/Effective_masses_sigmas']
            this_data.get('Correlators/Real').create_dataset('Effective_masses_sigmas', data=np.array(sigma_efm))
            
        
        if 'GEVP' not in this_data.keys(): j+=1; print('Irrep nr.: '+ str(j) + ' out of ' +str(len(the_list_name_irreps))); continue
        else:
            
            gevp_group = this_data.get('GEVP')
                
            print("Effective Masses of GEVP eigenvalues in process...")
                
            for item in gevp_group.keys():
                if 'Effective_masses' in gevp_group.get(item).keys(): del gevp_group[item+'/Effective_masses']
                
                group_em_t0 = gevp_group.get(item).create_group('Effective_masses')
                evalues_rs_f = np.array(gevp_group[item+'/Eigenvalues/Resampled'])
                evalues_mean_f = np.array(gevp_group[item+'/Eigenvalues/Mean'])
                
                eff_mass_rs_f=[]; eff_mass_mean = []; cov_eff_mass=[]; 
                for ls in range(the_size_matrix):
                    average = np.array(vf.EFF_MASS(evalues_mean_f[ls], the_dist_eff_mass),dtype=np.float64)
                    eff_mass_mean.append(average)
                    
                    eff_mass_rs = []
                    for zz in range(evalues_rs_f.shape[1]):
                        eff_mass_rs.append(np.array(vf.EFF_MASS(evalues_rs_f[ls][zz], the_dist_eff_mass), dtype=np.float64))
                    eff_mass_rs = np.array(vf.NCFGS_TO_NT(eff_mass_rs))
                    eff_rs_mean = vf.MEAN(np.array(eff_mass_rs))
                    cov_eff_mass.append(vf.STD_DEV_MEAN(eff_mass_rs, eff_rs_mean, the_type_rs))
                
                group_em_t0.create_dataset('Mean', data=np.array(eff_mass_mean))
                group_em_t0.create_dataset('Sigmas',data=np.array(cov_eff_mass))
                
                
            if kwargs.get('op_analysis')!=None and kwargs.get('op_analysis')!=False:                
                
                print("Effective Masses of NEW GEVP eigenvalues in process...")
                
                for op in range(the_size_matrix):
                    
                    print("Eigenvalue = %s"%str(op+1))
                    
                    gevp_group = this_data.get('Operators_Analysis/Op_%s'%op)
                    for item in gevp_group.keys():
                        if 'Effective_masses' in gevp_group.get(item).keys(): del gevp_group[item+'/Effective_masses']
                        
                        group_em_t0 = gevp_group.get(item).create_group('Effective_masses')
                        evalues_rs_f = np.array(gevp_group[item+'/Eigenvalues/Resampled'])
                        evalues_mean_f = np.array(gevp_group[item+'/Eigenvalues/Mean'])
                        
                        eff_mass_rs_f=[]; eff_mass_mean = []; cov_eff_mass=[]; 
                        for ls in range(the_size_matrix-1):
                            average = np.array(vf.EFF_MASS(evalues_mean_f[ls], the_dist_eff_mass),dtype=np.float64)
                            eff_mass_mean.append(average)
                            
                            eff_mass_rs = []
                            for zz in range(evalues_rs_f.shape[1]):
                                eff_mass_rs.append(np.array(vf.EFF_MASS(evalues_rs_f[ls][zz], the_dist_eff_mass), dtype=np.float64))
                            eff_mass_rs = np.array(vf.NCFGS_TO_NT(eff_mass_rs))
                            eff_rs_mean = vf.MEAN(np.array(eff_mass_rs))
                            cov_eff_mass.append(vf.STD_DEV_MEAN(eff_mass_rs, eff_rs_mean, the_type_rs))
                        
                        group_em_t0.create_dataset('Mean', data=np.array(eff_mass_mean))
                        group_em_t0.create_dataset('Sigmas',data=np.array(cov_eff_mass))
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
    myRb = 2
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
