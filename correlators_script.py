import numpy as np
import h5py
import time
import sys
import set_of_functions as vf


### ------------------------------- START FUNCTIONS ----------------------------------------------------

def SingleCorrelatorAnalysis(the_archivo, the_location, the_version, the_type_rs, the_irreps, the_weight, **kwargs):
    if kwargs.get('rebin_on')=='rb': 
        if kwargs.get('rb')==None:
            rb, the_re_bin = 1, '' 
            print('Missing argument: bin size. Using bin size equals %s.'%rb)
        else:
            rb = int(kwargs.get('rb'))
            the_re_bin = '_bin'+str(rb)
    else:
        rb, the_re_bin = 1, ''  
    
    if kwargs.get('nr_irreps')!=None:
        the_nr_irreps = int(kwargs.get('nr_irreps'))
    else:
        the_nr_irreps = len(the_irreps)    
    
    op = list(the_archivo[the_irreps[0]].attrs.keys())[1]
    
    if the_type_rs=='jk':
        the_resampling_scheme = 'Jackknife'
    elif the_type_rs=='bt':
        the_resampling_scheme = 'Bootstrap'
        if kwargs.get('kbt')==None:
            k_bt = 500
            print('Missing argument: Bootstrap sample size. Using sample size equals to %s'%k_bt)
        else:
            k_bt = int(kwargs.get('kbt'))
        if kwargs.get('own_kbt_list')==None:
            if kwargs.get('number_cfgs')==None:
                the_number_cnfgs = np.array(the_archivo[the_irreps[0]+'/data']).shape[0]
            else:
                the_number_cnfgs = int(kwargs.get('number_cfgs'))
            bt_cfgs = vf.RANDOM_GENERATOR(k_bt, int(the_number_cnfgs/rb))
        else:
            bt_cfgs = np.array(kwargs.get('own_kbt_list'))

    binned_rw = vf.BINNING(the_weight, rb)
    norm_reweight = vf.RW_NORMALIZATION(binned_rw, len(binned_rw))
    
    the_single_correlator_data = h5py.File(the_location+ '/Single_correlators_' + the_type_rs + the_re_bin + '_v%s.h5'%the_version,'w')
    
    begin_time = time.time()
    for j in range(the_nr_irreps):
        the_op_list = list(the_archivo[the_irreps[j]].attrs[op])
        the_size_matrix = len(the_op_list)
        
        datos_raw = np.array(the_archivo[the_irreps[j]+'/data'])
        
        the_times = str(the_archivo[the_irreps[j]].attrs['Other_Info']).split(' \n ')
        the_min_nt = int(the_times[0][the_times[0].index('= ')+2:])
        the_max_nt = int(the_times[1][the_times[1].index('= ')+2:])
        the_nt = np.arange(the_min_nt,the_max_nt+1)
        
        print('-->   IRREP (%s/'%str(j+1) + str(len(the_irreps)) +'): ', the_irreps[j])
        
        rw_datos = vf.REWEIGHTED_CORR(datos_raw.real, the_weight)
        if kwargs.get('rebin_on')=='rb':
            the_datos=[]
            for tt in range(len(rw_datos)):
                the_binned_datos = np.array(vf.BINNING(rw_datos[tt], rb))
                the_datos.append(np.array(the_binned_datos))
            the_datos = np.array(the_datos)
        else: 
            the_datos = np.array(rw_datos)
        if the_type_rs=='jk':
            the_rs = []
            for tt in range(len(the_datos)):
                the_rs.append(vf.JACKKNIFE(the_datos[tt], norm_reweight))
        elif the_type_rs=='bt':
            the_rs=[]
            for tt in range(len(the_datos)):
                the_rs.append(np.array(vf.BOOTSTRAP(the_datos[tt], bt_cfgs, norm_reweight)))
        the_rs = np.array(the_rs)
        
        print('----------------------------------------------\n               DATA SHAPE \n----------------------------------------------\nNr. of gauge configurations: ' +  str(len(the_datos[0])) + '\n' + 'Time slices: '+str(the_nt[0])+' to '+str(the_nt[-1]) +'\nResampling data (%s): '%the_resampling_scheme + str(len(the_rs[0])) + '\n....................................')
        print('      OPERATORS LIST ')
        i=0
        for i in range(the_size_matrix):
            print('      -->>  '+str(the_op_list[i]))
        g_i = the_single_correlator_data.create_group(the_irreps[j]) 
        
        g_i.create_dataset('Time_slices', data=the_nt)
        group_corr= g_i.create_group('Correlators')
        group_corr_real = group_corr.create_group('Real')
        g_i.create_dataset('Operators', data=the_op_list) 
        
        mrs_f_rs = np.array(vf.MEAN(the_rs))
        
        mrs_f = np.array(vf.MEAN(the_datos))
        
        sigma_corr = np.array(vf.STD_DEV_MEAN(the_rs, mrs_f_rs, the_type_rs))
        
        cov_corr = np.array(vf.COV_MATRIX(the_rs, mrs_f_rs, the_type_rs))
        
        group_corr_real.create_dataset('Mean', data = mrs_f) 
        group_corr_real.create_dataset('Sigmas', data = sigma_corr)
        group_corr_real.create_dataset('Resampled', data = the_rs) 
        group_corr_real.create_dataset('Covariance_matrix', data = cov_corr)
        j+=1
        print('Irrep nr.: '+ str(j) + ' out of ' +str(len(the_irreps)))
    the_single_correlator_data.close()
    end_time = time.time()
    print('TIME TAKEN: ' + str((end_time-begin_time)/60) +' mins')
    return the_location + '/Single_correlators_' + the_type_rs + the_re_bin + '_v%s.h5'%the_version
    


def MultiCorrelatorAnalysis(the_archivo, the_location, the_version, the_type_rs, the_irreps, the_weight, **kwargs):
    if kwargs.get('rebin_on')=='rb': 
        if kwargs.get('rb')==None:
            rb, the_re_bin = 1, '' 
            print('Missing argument: bin size. Using bin size equals %s.'%rb)
        else:
            rb = int(kwargs.get('rb'))
            the_re_bin = '_bin'+str(rb)
    else:
        rb, the_re_bin = 1, ''  

    op = list(the_archivo[the_irreps[0]].attrs.keys())[1]
    
    if kwargs.get('nr_irreps')!=None:
        the_nr_irreps = int(kwargs.get('nr_irreps'))
    else:
        the_nr_irreps = len(the_irreps)    

    if the_type_rs=='jk':
        the_resampling_scheme = 'Jackknife'
    elif the_type_rs=='bt':
        the_resampling_scheme = 'Bootstrap'
        if kwargs.get('kbt')==None:
            k_bt = 500
            print('Missing argument: Bootstrap sample size. Using sample size equals to %s'%k_bt)
        else:
            k_bt = int(kwargs.get('kbt'))
        if kwargs.get('number_cfgs')==None:
            the_number_cnfgs = np.array(the_archivo[the_irreps[0]+'/data']).shape[0]
        else:
            the_number_cnfgs = int(kwargs.get('number_cfgs'))
        bt_cfgs = vf.RANDOM_GENERATOR(k_bt, int(the_number_cnfgs/rb))

    binned_rw = vf.BINNING(the_weight, rb)
    norm_reweight = vf.RW_NORMALIZATION(binned_rw, len(binned_rw))
    
    the_matrix_correlator_data = h5py.File(the_location + '/Matrix_correlators_' + the_type_rs + the_re_bin + '_v%s.h5'%the_version,'w')    
    
    begin_time = time.time()
    
    for j in range(the_nr_irreps):
        the_op_list = list(the_archivo[the_irreps[j]].attrs[op])
        the_size_matrix = len(the_op_list)
        datos_raw = np.array(the_archivo[the_irreps[j]+'/data'])    
            
        the_times = str(the_archivo[the_irreps[j]].attrs['Other_Info']).split(' \n ')
        the_min_nt = int(the_times[0][the_times[0].index('= ')+2:])
        the_max_nt = int(the_times[1][the_times[1].index('= ')+2:])
        the_nt = np.arange(the_min_nt, the_max_nt+1)
        
        print('\n----------------------------------------------')
        print(' -->   IRREP (%s/'%str(j+1) + str(len(the_irreps)) +'): ', the_irreps[j])
        
        re_datos = vf.RESHAPING(datos_raw.real ,the_size_matrix)
        the_datos_n1_n2=[]
        if kwargs.get('rebin_on')=='rb':
            for n1 in range(the_size_matrix):
                the_datos_n1=[]
                for n2 in range(the_size_matrix):
                    rw_datos =  vf.REWEIGHTED_CORR(re_datos[n1][n2], the_weight)
                    rw_t_datos = []
                    for tt in range(len(rw_datos)):
                        rw_t_datos.append(np.array(vf.BINNING(rw_datos[tt], rb)))
                    the_datos_n1.append(np.array(rw_t_datos))
                the_datos_n1_n2.append(np.array(the_datos_n1))
            the_datos = np.array(the_datos_n1_n2,dtype=np.float64)
        else:
            for n1 in range(the_size_matrix):
                the_datos_n1=[]
                for n2 in range(the_size_matrix):
                    rw_datos =  vf.REWEIGHTED_CORR(re_datos[n1][n2], the_weight)
                    the_datos_n1.append(np.array(rw_datos))
                the_datos_n1_n2.append(np.array(the_datos_n1))
            the_datos = np.array(the_datos_n1_n2,dtype=np.float64)
        
        if the_type_rs=='jk':
            the_rs=[]
            for n1 in range(the_size_matrix):
                rs_n1=[]
                for n2 in range(the_size_matrix):
                    dis_data = np.array(the_datos[n1][n2])
                    rs_n1_t = []
                    for tt in range(len(dis_data)):
                        rs_n1_t.append(vf.JACKKNIFE(dis_data[tt], norm_reweight))
                    rs_n1.append(np.array(rs_n1_t))
                the_rs.append(rs_n1)
            the_rs = np.array(the_rs)
        elif the_type_rs=='bt':
            the_rs=[]
            for n1 in range(the_size_matrix):
                rs_n1=[]
                for n2 in range(the_size_matrix):
                    dis_data = np.array(the_datos[n1][n2])
                    rs_n1_t = []
                    for tt in range(len(dis_data)):
                        rs_n1_t.append(vf.BOOTSTRAP(dis_data[tt], bt_cfgs, norm_reweight))
                    rs_n1.append(np.array(rs_n1_t))
                the_rs.append(rs_n1)
            the_rs = np.array(the_rs)     
        
        print('----------------------------------------------\n               DATA SHAPE \n----------------------------------------------\nNr. of gauge configurations: ' +  str(the_datos.shape[-1]) + '\nSize of the Correlation matrix: ' + str(the_size_matrix)+ 'x' + str(the_size_matrix) +  '\nTime slices: '+str(the_nt[0])+' to '+str(the_nt[-1]) + '\nResampling data (%s): '%the_resampling_scheme+ str(the_rs.shape[-1])+ '\n....................................')
        print('      OPERATORS LIST ')
        for i in range(the_size_matrix):
            print('      -->>  '+str(the_op_list[i]))
                
        g_i = the_matrix_correlator_data.create_group(the_irreps[j]) 
        g_i.create_dataset('Time_slices', data=the_nt)
        g_i.create_dataset('Operators', data=the_op_list) 
        group_corr = g_i.create_group('Correlators')
        group_corr_real = group_corr.create_group('Real')

        mrs_f_real=[]
        for n1 in range(the_size_matrix):
            mrs_f_real_n1=[]
            for n2 in range(the_size_matrix):
                mrs_f_real_n1.append(np.array(vf.MEAN(the_datos[n1][n2])))
            mrs_f_real.append(np.array(mrs_f_real_n1))
        mrs_f = np.array(mrs_f_real)
        
        re_rs = vf.RESHAPING_EIGENVALS_RS(the_rs, the_size_matrix)
        re_mean = vf.RESHAPING_EIGENVALS_MEAN(mrs_f, the_size_matrix)

        group_corr_real.create_dataset('Resampled',data=re_rs)
        group_corr_real.create_dataset('Mean', data= re_mean)
        j+=1
    the_matrix_correlator_data.close()
    end_time = time.time()
    print('TIME TAKEN: ' + str((end_time-begin_time)/60) +' mins')    
    return the_location + '/Matrix_correlators_' + the_type_rs + the_re_bin + '_v%s.h5'%the_version



### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------


if __name__== "__main__":
    myEns = str(sys.argv[1]).upper()
    myWhichCorrelator = str(sys.argv[2]).lower()
    myTypeRs = str(sys.argv[3]).lower()
    myRebinOn = str(sys.argv[4]).lower()
    myRb = 2
    myVersion = 'test'
    myKbt = 500
    myNrIrreps=1
    
    if myEns == 'N451': from files_n451 import *
    elif myEns == 'N201': from files_n201 import * 
    elif myEns == 'D200': from files_d200 import *
    
    myWeight = weight
    myLocation = vf.DIRECTORY_EXISTS(location + '$YOUR_OUTPUT_PATH$/%s/'%myEns)
    myCnfgs = ncfgs
    
    vf.INFO_PRINTING(myWhichCorrelator, myEns)
    
    if myWhichCorrelator =='s':        
        myArchivo, myIrreps = f1, name1 
        correlatorAnalysis = SingleCorrelatorAnalysis 
    
    elif myWhichCorrelator=='m':
        myArchivo, myIrreps = f, name
        correlatorAnalysis = MultiCorrelatorAnalysis
    else:
        print('NOT AN OPTION.\nQUITTING.')
        sys.exit()
    
    savedLocation = correlatorAnalysis(myArchivo, myLocation, myVersion, myTypeRs, myIrreps, myWeight, rebin_on = myRebinOn, rb = myRb, kbt = myKbt, number_cfgs = myCnfgs, nr_irreps=myNrIrreps)
    
    myArchivo.close()
    print('-'*(len(savedLocation)+1))
    print('Saved as: \n' + savedLocation)
    print('_'*(len(savedLocation)+1))
