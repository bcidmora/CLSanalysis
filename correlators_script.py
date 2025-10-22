import numpy as np
import h5py
import time
import sys
import set_of_functions as vf


### ------------------------------- START FUNCTIONS ----------------------------------------------------

def SingleCorrelatorAnalysis(the_archivo, the_location, the_version, the_type_rs, the_irreps, the_weight, **kwargs):
    
    print("                     CORRELATORS ANALYSIS \n")
    
    ### It chooses the rebin
    if kwargs.get('rebin_on')=='rb': 
        if kwargs.get('rb')==None:
            rb, the_re_bin = 1, '' 
            print('Missing argument: bin size. Using bin size equals %s.'%rb)
        else:
            rb = int(kwargs.get('rb'))
            the_re_bin = '_bin'+str(rb)
    else:
        rb, the_re_bin = 1, ''  
    
    ### How many irreps do you want to study    
    the_s_irreps = list(the_archivo.keys())
    ### How many irreps do you want to study    
    if kwargs.get('nr_irreps')!=None:
        the_first_irrep = 0
        the_last_irrep = int(kwargs.get('nr_irreps'))+1
    else:
        ### This one checks for an irrep in particular
        if kwargs.get('first_irrep')!=None and kwargs.get('last_irrep')!=None:
            the_first_irrep = int(kwargs.get('first_irrep'))-1
            the_last_irrep = int(kwargs.get('last_irrep'))
        elif kwargs.get('last_irrep')!=None and kwargs.get('first_irrep')==None:
            the_first_irrep = 0
            the_last_irrep = int(kwargs.get('last_irrep'))
        elif kwargs.get('first_irrep')!=None and kwargs.get('last_irrep')==None:
            the_first_irrep = int(kwargs.get('first_irrep'))-1
            the_last_irrep = len(the_s_irreps)
        elif kwargs.get('first_irrep')==None and kwargs.get('last_irrep')==None:
            the_first_irrep = 0
            the_last_irrep = len(the_s_irreps)
    
    s_irreps = the_s_irreps[the_first_irrep:the_last_irrep]
    
    ### If not all configs, this can be changed.
    if kwargs.get('number_cfgs')==None:
        the_number_cnfgs = np.array(the_archivo[the_irreps[0]+'/data']).shape[0]
    else:
        the_number_cnfgs = int(kwargs.get('number_cfgs'))
    
    ### The list of operators
    the_op = list(the_archivo[the_irreps[0]].attrs.keys())[1]
    
    ### Resampling scheme
    if the_type_rs=='jk':
        the_resampling_scheme = 'Jackknife'
    elif the_type_rs=='bt':
        the_resampling_scheme = 'Bootstrap'
        if kwargs.get('kbt')==None:
            k_bt = 500 # Default
            print('Missing argument: Bootstrap sample size. Using sample size equals to %s'%k_bt)
        else:
            k_bt = int(kwargs.get('kbt'))
        if kwargs.get('own_kbt_list')==None: ### You can also choose random numbers from a list
            bt_cfgs = vf.RANDOM_GENERATOR(k_bt, int(the_number_cnfgs/rb))
        else:
            bt_cfgs = np.array(kwargs.get('own_kbt_list'))

    ### The reweighting factors are renormalized according to the configs taken.
    the_weight = the_weight[:the_number_cnfgs]
    
    ### The binning now also includes the reweighting factors
    binned_rw = vf.BINNING(the_weight, rb)
    
    ### Normalized again with the binning
    norm_reweight = vf.RW_NORMALIZATION(binned_rw, len(binned_rw))
    
    ### This is the single correlators data
    the_single_correlator_data = h5py.File(the_location+ '/Single_correlators_' + the_type_rs + the_re_bin + '_v%s.h5'%the_version,'w')
    
    begin_time = time.time()
    ### Start of the analysis for the nr. of irreps.
    for the_irrep in s_irreps:
        
        ### The list of operators 
        the_op_list = list(the_archivo[the_irrep].attrs[the_op])
        
        ### The size of the matrix for the single hadrons is always 1
        the_size_matrix = len(the_op_list)
        
        ### Extracting the original/raw data for the analysis
        the_datos_raw = np.array(the_archivo[the_irrep+'/data'])[:the_number_cnfgs]
        
        ### Time slices
        the_times = str(the_archivo[the_irrep].attrs['Other_Info']).split(' \n ')
        
        ### Min time slice
        the_min_nt = int(the_times[0][the_times[0].index('= ')+2:])
        
        ### Max time slice
        the_max_nt = int(the_times[1][the_times[1].index('= ')+2:])
        
        ### Total range
        the_nt = np.arange(the_min_nt,the_max_nt+1)
        
        print('-->   IRREP (%s/'%str(the_s_irreps.index(the_irrep)+1) + str(len(the_irreps)) +'): ', the_irrep)
        
        ### Reweighted data set
        rw_datos = vf.REWEIGHTED_CORR(the_datos_raw.real, the_weight)
        
        ### If there is binning, then this applies to the  data and to the reweighting factors
        if kwargs.get('rebin_on')=='rb':
            the_datos=[]
            for tt in range(len(rw_datos)):
                the_binned_datos = np.array(vf.BINNING(rw_datos[tt], rb))
                the_datos.append(np.array(the_binned_datos))
            the_datos = np.array(the_datos)
        else: 
            the_datos = np.array(rw_datos)
        
        ### Resampling must be done with the normalized reweighting factors
        if the_type_rs=='jk':
            the_rs = []
            for tt in range(len(the_datos)):
                the_rs.append(vf.JACKKNIFE(the_datos[tt], norm_reweight))
        elif the_type_rs=='bt':
            the_rs=[]
            for tt in range(len(the_datos)):
                the_rs.append(np.array(vf.BOOTSTRAP(the_datos[tt], bt_cfgs, norm_reweight)))
                
        ### The resampled data
        the_rs = np.array(the_rs)
        
        ### Information about the ongoing analysis
        print('----------------------------------------------\n               DATA SHAPE \n----------------------------------------------\nNr. of gauge configurations: ' +  str(len(the_datos[0])) + '\n' + 'Time slices: '+str(the_nt[0])+' to '+str(the_nt[-1]) +'\nResampling data (%s): '%the_resampling_scheme + str(len(the_rs[0])) + '\n....................................')
        print('      OPERATORS LIST ')
        
        ### Printing operators
        i=0
        for i in range(the_size_matrix):
            print('      -->>  '+str(the_op_list[i]))
            
        
        g_i = the_single_correlator_data.create_group(the_irrep) 
        
        g_i.create_dataset('Time_slices', data=the_nt)
        group_corr= g_i.create_group('Correlators')
        group_corr_real = group_corr.create_group('Real')
        g_i.create_dataset('Operators', data=the_op_list) 
        
        ### This is the mean value of the resampled correlators
        the_mrs_f_rs = np.array(vf.MEAN(the_rs))
        
        ### This is the mean value of the original/raw correlators
        the_mrs_f = np.array(vf.MEAN(the_datos))
        
        ### The statistical error of the resampled data
        the_sigma_corr = np.array(vf.STD_DEV_MEAN(the_rs, the_mrs_f_rs, the_type_rs))
        
        ### The covariance matrix
        the_cov_corr = np.array(vf.COV_MATRIX(the_rs, the_mrs_f_rs, the_type_rs))
        
        group_corr_real.create_dataset('Mean', data = the_mrs_f) 
        group_corr_real.create_dataset('Sigmas', data = the_sigma_corr)
        group_corr_real.create_dataset('Resampled', data = the_rs) 
        group_corr_real.create_dataset('Covariance_matrix', data = the_cov_corr)
    the_single_correlator_data.close()
    end_time = time.time()
    print('TIME TAKEN: ' + str((end_time-begin_time)/60) +' mins')
    return the_location + '/Single_correlators_' + the_type_rs + the_re_bin + '_v%s.h5'%the_version
    

def MultiCorrelatorAnalysis(the_archivo, the_quantum_number, the_location, the_version, the_type_rs, the_irreps, the_weight, **kwargs):
    
    print("                     CORRELATORS ANALYSIS \n")
    
    ### It chooses the rebin
    if kwargs.get('rebin_on')=='rb': 
        if kwargs.get('rb')==None:
            rb, the_re_bin = 1, '' 
            print('Missing argument: bin size. Using bin size equals %s.'%rb)
        else:
            rb = int(kwargs.get('rb'))
            the_re_bin = '_bin'+str(rb)
    else:
        rb, the_re_bin = 1, ''  
    
    the_m_irreps = list(the_archivo.keys())
    ### How many irreps do you want to study    
    if kwargs.get('nr_irreps')!=None:
        the_first_irrep = 0
        the_last_irrep = int(kwargs.get('nr_irreps'))
    else:
        ### This one checks for an irrep in particular
        if kwargs.get('first_irrep')!=None and kwargs.get('last_irrep')!=None:
            the_first_irrep = int(kwargs.get('first_irrep'))-1
            the_last_irrep = int(kwargs.get('last_irrep'))
        elif kwargs.get('last_irrep')!=None and kwargs.get('first_irrep')==None:
            the_first_irrep = 0
            the_last_irrep = int(kwargs.get('last_irrep'))
        elif kwargs.get('first_irrep')!=None and kwargs.get('last_irrep')==None:
            the_first_irrep = int(kwargs.get('first_irrep'))-1
            the_last_irrep = len(the_m_irreps)
        elif kwargs.get('first_irrep')==None and kwargs.get('last_irrep')==None:
            the_first_irrep = 0
            the_last_irrep = len(the_m_irreps)
    
    m_irreps = the_m_irreps[the_first_irrep:the_last_irrep]
    
    ### If not all configs, this can be changed.
    if kwargs.get('number_cfgs')==None:
        the_number_cnfgs = np.array(the_archivo[the_irreps[0]+'/data']).shape[0]
    else:
        the_number_cnfgs = int(kwargs.get('number_cfgs'))

    ### The list of operators
    the_op = list(the_archivo[the_irreps[0]].attrs.keys())[1]
    
    ### Resampling scheme
    if the_type_rs=='jk':
        the_resampling_scheme = 'Jackknife'
    elif the_type_rs=='bt':
        the_resampling_scheme = 'Bootstrap'
        if kwargs.get('kbt')==None:
            k_bt = 500
            print('Missing argument: Bootstrap sample size. Using sample size equals to %s'%k_bt)
        else:
            k_bt = int(kwargs.get('kbt'))
        bt_cfgs = vf.RANDOM_GENERATOR(k_bt, int(the_number_cnfgs/rb))
    
    ### The reweighting factors are renormalized according to the configs taken.
    the_weight = the_weight[:the_number_cnfgs]

    ### The binning now also includes the reweighting factors
    binned_rw = vf.BINNING(the_weight, rb)
    
    ### Normalized again with the binning
    norm_reweight = vf.RW_NORMALIZATION(binned_rw, len(binned_rw))
    
    ### This is the single correlators data
    the_matrix_correlator_data = h5py.File(the_location + '/Matrix_correlators' + the_quantum_number + the_type_rs + the_re_bin + '_v%s.h5'%the_version,'w')      
    
    begin_time = time.time()
    ### Start of the analysis for the nr. of irreps.
    for the_irrep in m_irreps:
         ### The list of operators 
        the_op_list = list(the_archivo[the_irrep].attrs[the_op])
        
        ### The size of the matrix
        the_size_matrix = len(the_op_list)
        
        ### Extracting the original/raw data for the analysis
        the_datos_raw = np.array(the_archivo[the_irrep+'/data'])[:the_number_cnfgs]
        
        ### Time slices
        the_times = str(the_archivo[the_irrep].attrs['Other_Info']).split(' \n ')
        
        ### Min time slice
        the_min_nt = int(the_times[0][the_times[0].index('= ')+2:])
        
        ### Max time slice
        the_max_nt = int(the_times[1][the_times[1].index('= ')+2:])
        
         ### Total range
        the_nt = np.arange(the_min_nt, the_max_nt+1)
        
        print('\n----------------------------------------------')
        print(' -->   IRREP (%s/'%str(the_m_irreps.index(the_irrep)+1) + str(len(the_irreps)) +'): ', the_irrep)
        
        ### Reweighted data set
        re_datos = vf.RESHAPING(the_datos_raw.real)
        
        ### If there is binning, then this applies to the  data and to the reweighting factors
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
        
        ### Resampling must be done with the normalized reweighting factors
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
        
        ### The resampled data    
        the_rs = np.array(the_rs)     
        
        print('----------------------------------------------\n               DATA SHAPE \n----------------------------------------------\nNr. of gauge configurations: ' +  str(the_datos.shape[-1]) + '\nSize of the Correlation matrix: ' + str(the_size_matrix)+ 'x' + str(the_size_matrix) +  '\nTime slices: '+str(the_nt[0])+' to '+str(the_nt[-1]) + '\nResampling data (%s): '%the_resampling_scheme+ str(the_rs.shape[-1])+ '\n....................................')
        print('      OPERATORS LIST ')
        for i in range(the_size_matrix):
            print('      -->>  '+str(the_op_list[i]))
                
        g_i = the_matrix_correlator_data.create_group(the_irrep) 
        g_i.create_dataset('Time_slices', data=the_nt)
        g_i.create_dataset('Operators', data=the_op_list) 
        group_corr = g_i.create_group('Correlators')
        group_corr_real = group_corr.create_group('Real')
        
        ### Calculating the mean values of the datasets
        the_mrs_f_real, the_rs_mean_real = [], []
        for n1 in range(the_size_matrix):
            the_mrs_f_real_n1, the_rs_mean_real_n1 = [], []
            for n2 in range(the_size_matrix):
                ### This is the mean value of the original/raw correlators
                the_mrs_f_real_n1.append(np.array(vf.MEAN(the_datos[n1][n2])))
                ### This is the mean value of the resampled correlators
                the_rs_mean_real_n1.append(np.array(vf.MEAN(the_rs[n1][n2])))
            the_mrs_f_real.append(np.array(the_mrs_f_real_n1))
            the_rs_mean_real.append(np.array(the_rs_mean_real_n1))
        
        
        the_mrs_f = np.array(the_mrs_f_real)
        the_mrs_f_rs = np.array(the_rs_mean_real)
        
        ### The statistical error of the resampled data
        the_sigmas_corr = []
        for ss in range(the_size_matrix):
            the_sigmas_corr.append(np.array(vf.STD_DEV_MEAN(the_rs[ss][ss], the_mrs_f_rs[ss][ss], the_type_rs)))
        the_sigmas_corr=np.array(the_sigmas_corr)
        
        ### Reshaping the data for later diagonizalization and extraction of eigenvalues
        re_rs = vf.RESHAPING_EIGENVALS_RS(the_rs)
        re_mean = vf.RESHAPING_EIGENVALS_MEAN(the_mrs_f)
        re_mean_rs = vf.RESHAPING_EIGENVALS_MEAN(the_mrs_f_rs)
        
        group_corr_real.create_dataset('Resampled',data=re_rs)
        group_corr_real.create_dataset('Mean', data= re_mean)
        group_corr_real.create_dataset('Sigmas', data= the_sigmas_corr)
    the_matrix_correlator_data.close()
    end_time = time.time()
    print('TIME TAKEN: ' + str((end_time-begin_time)/60) +' mins')    
    return the_location + '/Matrix_correlators' + the_quantum_number + the_type_rs + the_re_bin + '_v%s.h5'%the_version

### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------



### ------------------------------- START EXECUTING --------------------------------------------------



if __name__== "__main__":
    
    ### The ensemble
    myEns = str(sys.argv[1]).upper()
    
    ### Single hadrons or multihadrons
    myWhichCorrelator = str(sys.argv[2]).lower()
    
    ### Type of resampling 'bt' or 'jk'
    myTypeRs = str(sys.argv[3]).lower()
    
    ### Rebinning or not
    myRebinOn = str(sys.argv[4]).lower()
    myRb = 1
    
    ### Name of the output file
    myVersion = '_test'
    
    ### Default bootstrap sampling
    myKbt = 500
    
    ### If you don't want to start with the very first irrep
    myNrIrreps = None # 2 # 1
    myFirstIrrep = None # 1 # 2
    myLastIrrep = None
    
    if myEns == 'N451': from files_n451 import *
    elif myEns == 'N201': from files_n201 import * 
    elif myEns == 'D200': from files_d200 import *
    elif myEns == 'X451': from files_x451 import *
    
    ### This information comes from the files_ens.py, so be careful how you define your directories.
    myWeight = weight
    myLocation = vf.DIRECTORY_EXISTS(location + '$YOUR_OUTPUT_PATH$/%s/'%myEns)
    myCnfgs = ncfgs
    myChosenIsospin = the_hadron_state
    
    vf.INFO_PRINTING(myWhichCorrelator, myEns)
    
    if myWhichCorrelator =='s':        
        myArchivo, myIrreps = f1, name1 
        
        savedLocation = SingleCorrelatorAnalysis(myArchivo, myLocation, myVersion, myTypeRs, myIrreps, myWeight, rebin_on = myRebinOn, rb = myRb, kbt = myKbt, number_cfgs = myCnfgs, nr_irreps=myNrIrreps, own_kbt_list = myKbtSamples, first_irrep = myFirstIrrep , last_irrep = myLastIrrep)
    
    elif myWhichCorrelator=='m':
        myArchivo, myIrreps = f, name
        
        savedLocation = MultiCorrelatorAnalysis(myArchivo, myChosenIsospin, myLocation, myVersion, myTypeRs, myIrreps, myWeight, rebin_on = myRebinOn, rb = myRb, kbt = myKbt, number_cfgs = myCnfgs, nr_irreps=myNrIrreps, own_kbt_list = myKbtSamples, first_irrep = myFirstIrrep , last_irrep = myLastIrrep)
    else:
        print('NOT AN OPTION.\nQUITTING.')
        sys.exit()
    
    myArchivo.close()
    print('-'*(len(savedLocation)+1))
    print('Saved as: \n' + savedLocation)
    print('_'*(len(savedLocation)+1))
