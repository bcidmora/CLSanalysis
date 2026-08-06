import numpy as np
import h5py
import time
import sys
import set_of_analysis_functions as vfa
import set_of_layout_functions as vfl


### ------------------------------- START FUNCTIONS ----------------------------------------------------

def SingleCorrelatorAnalysis(the_archivo, the_location, the_version, the_type_rs, the_irreps, the_weight, **kwargs):
    
    print("                     CORRELATORS ANALYSIS \n")
    
    ### It chooses the rebin
    if kwargs.get('rebin_on'): 
        if kwargs.get('rb')==None:
            rb, the_re_bin = 1, '' 
            print(f'Missing argument: bin size. Using bin size equals {rb}.')
        else:
            rb = int(kwargs.get('rb'))
            the_re_bin = f'_bin{rb}'
    else:
        rb, the_re_bin = 1, ''  

    ### If not all configs, this can be changed.
    if kwargs.get('number_cfgs')==None:
        the_number_cnfgs = np.asarray(the_archivo[f'{the_irreps[0]}/data']).shape[0]
    else:
        the_number_cnfgs = int(kwargs.get('number_cfgs'))
    
    ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(the_irreps)

    s_irreps = the_irreps[the_first_irrep:the_last_irrep]
    
    ### Resampling scheme
    if the_type_rs=='jk':
        the_resampling_scheme = 'Jackknife'
        if rb == 'C':
            jackknife_choice = vfa.JACKKNIFE_CUSTOM
        else: jackknife_choice = vfa.JACKKNIFE
    elif the_type_rs=='bt':
        the_resampling_scheme = 'Bootstrap'
        if kwargs.get('kbt')==None:
            k_bt = 500 # Default
            print(f'Missing argument: Bootstrap sample size. Using sample size equals to {k_bt}')
        else:
            k_bt = int(kwargs.get('kbt'))
        if kwargs.get('own_kbt_list') is None: ### You can also choose random numbers from a list
            bt_cfgs = vfa.RANDOM_GENERATOR(k_bt, int(the_number_cnfgs/rb))
        else:
            bt_cfgs = np.array(kwargs.get('own_kbt_list'))

    ### The reweighting factors are renormalized according to the configs taken.
    the_weight = the_weight[:the_number_cnfgs]
    
    ### Special binning routine
    if rb == 'C':
        binned_rw = vfa.BINNING_CUSTOM_CORR_WEIGHTS(the_weight)
        ### Normalized again with the binning
        norm_reweight = vfa.RW_NORMALIZATION(binned_rw, len(binned_rw))
        
    else:
        ### The binning now also includes the reweighting factors
        binned_rw = vfa.BINNING(the_weight, rb)
        
        ### Normalized again with the binning
        norm_reweight = vfa.RW_NORMALIZATION(binned_rw, len(binned_rw))
    
    ### This is the single correlators data
    the_single_correlator_data = h5py.File(f'{the_location}/Single_correlators_{the_type_rs}{the_re_bin}_{the_version}.h5','w')
    
    begin_time = time.time()
    ### Start of the analysis for the nr. of irreps.
    for the_irrep in s_irreps:
        
        ### The list of operators 
        the_op_list = list(the_archivo[the_irrep].attrs['op_list'])
        
        ### The size of the matrix for the single hadrons is always 1
        the_size_matrix = len(the_op_list)
        
        ### Extracting the original/raw data for the analysis
        the_datos_raw = np.asarray(the_archivo[f'{the_irrep}/data'])[:the_number_cnfgs]
        
        ### Time slices
        the_times = str(the_archivo[the_irrep].attrs['Other_Info']).split(' \n ')
        
        ### Min time slice
        the_min_nt = int(the_times[0][the_times[0].index('= ') + 2:])
        
        ### Max time slice
        the_max_nt = int(the_times[1][the_times[1].index('= ') + 2:])
        
        ### Total range
        the_nt = np.arange(the_min_nt, the_max_nt + 1)
        
        print('\n----------------------------------------------')
        print(f'-->   IRREP ({the_irreps.index(the_irrep)+1}/{len(the_irreps)}): {the_irrep}')
        
        ### Reweighted data set (keeping it complex)
        rw_datos = vfa.REWEIGHTED_CORR(the_datos_raw, the_weight)
        
        ### If there is binning, then this applies to the  data and to the reweighting factors
        if kwargs.get('rebin_on'):
            if rb == 'C':
                the_binned_datos, bin_weights = zip(*(vfa.BINNING_CUSTOM(rw_datos[tt]) for tt in range(len(rw_datos))))
                the_datos = np.asarray(the_binned_datos)
                bins = list(bin_weights)
            else:
                bins = None
                the_datos = np.asarray([vfa.BINNING(tt, rb) for tt in rw_datos])
        else: 
            bins = None
            the_datos = np.asarray(rw_datos)

        ### Resampling must be done with the normalized reweighting factors
        if the_type_rs=='jk':
            bin_rw = bins if bins is not None else None
            the_rs = np.asarray([jackknife_choice(tt, norm_reweight, bin_rw) for tt in the_datos])
        elif the_type_rs=='bt':
            the_rs = np.asarray([vfa.BOOTSTRAP(tt, bt_cfgs, norm_reweight) for tt in the_datos])
                
        ### Information about the ongoing analysis
        print(f'----------------------------------------------\n               DATA SHAPE \n----------------------------------------------\n    Nr. of gauge configurations: {len(the_datos[0])} \n     Time slices: {the_nt[0]} to {the_nt[-1]} \n     Resampling data ({the_resampling_scheme}): {len(the_rs[0])}\n....................................')
        print('      OPERATORS LIST ')
        
        ### Printing operators
        for i in range(the_size_matrix):
            print(f'      -->>  {the_op_list[i]}')
        
        g_i = the_single_correlator_data.create_group(the_irrep) 
        
        g_i.create_dataset('Time_slices', data = the_nt)
        group_corr= g_i.create_group('Correlators')
        group_corr_real = group_corr.create_group('Real')
        g_i.create_dataset('Operators', data = the_op_list) 
        
        ### This is the mean value of the resampled correlators
        the_mrs_f_rs = np.asarray(vfa.MEAN(the_rs))
        
        ### This is the mean value of the original/raw correlators
        if rb == 'C':
            the_mrs_f = np.asarray(vfa.MEAN_CUSTOM(the_datos, bins))
        else:
            the_mrs_f = np.asarray(vfa.MEAN(the_datos))
        
        ### The statistical error of the resampled data
        the_sigma_corr = np.asarray(vfa.STD_DEV_MEAN(the_rs.real, the_mrs_f_rs.real, the_type_rs))
        
        ### The covariance matrix
        the_cov_corr = np.asarray(vfa.COV_MATRIX(the_rs.real, the_mrs_f_rs.real, the_type_rs))
        
        group_corr_real.create_dataset('Mean', data = the_mrs_f) 
        group_corr_real.create_dataset('Sigmas', data = the_sigma_corr)
        group_corr_real.create_dataset('Resampled', data = the_rs) 
        group_corr_real.create_dataset('Covariance_matrix', data = the_cov_corr)
    the_single_correlator_data.close()
    end_time = time.time()
    print(f'TIME TAKEN: {(end_time-begin_time)/60} mins')
    return  f'{the_location}Single_correlators_{the_type_rs}{the_re_bin}_{the_version}.h5'
    


def MultiCorrelatorAnalysis(the_archivo, the_location, the_version, the_type_rs, the_irreps, the_weight, **kwargs):
    
    print("                     CORRELATORS ANALYSIS \n")
    
    ### It chooses the rebin
    if kwargs.get('rebin_on'): 
        if kwargs.get('rb')==None:
            rb, the_re_bin = 1, '' 
            print(f'Missing argument: bin size. Using bin size equals {rb}.')
        else:
            rb = int(kwargs.get('rb'))
            the_re_bin = f'_bin{rb}'
    else:
        rb, the_re_bin = 1, ''  
    
    ### If not all configs, this can be changed.
    if kwargs.get('number_cfgs')==None:
        the_number_cnfgs = np.array(the_archivo[f'{the_irreps[0]}/data']).shape[0]
    else:
        the_number_cnfgs = int(kwargs.get('number_cfgs'))
        
    ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(the_irreps)
    
    m_irreps = the_irreps[the_first_irrep:the_last_irrep]
    
    ### Resampling scheme
    if the_type_rs=='jk':
        the_resampling_scheme = 'Jackknife'
        if rb == 'C':
            jackknife_choice = vfa.JACKKNIFE_CUSTOM
        else: 
            jackknife_choice = vfa.JACKKNIFE
    elif the_type_rs=='bt':
        the_resampling_scheme = 'Bootstrap'
        if kwargs.get('kbt')==None:
            k_bt = 500 # Default
            print(f'Missing argument: Bootstrap sample size. Using sample size equals to {k_bt}')
        else:
            k_bt = int(kwargs.get('kbt'))
        bt_cfgs = vfa.RANDOM_GENERATOR(k_bt, int(the_number_cnfgs/rb))
    
    ### The reweighting factors are renormalized according to the configs taken.
    the_weight = the_weight[:the_number_cnfgs]
    
    ### The binning now also includes the reweighting factors
    if rb == 'C':
        ### The binning now also includes the reweighting factors
        binned_rw = vfa.BINNING_CUSTOM_CORR_WEIGHTS(the_weight)

        ### Normalized again with the binning
        norm_reweight = vfa.RW_NORMALIZATION(binned_rw, len(binned_rw))
    else:
        ### The binning now also includes the reweighting factors
        binned_rw = vfa.BINNING(the_weight, rb)

        ### Normalized again with the binning
        norm_reweight = vfa.RW_NORMALIZATION(binned_rw, len(binned_rw))
    
    ### This is the single correlators data
    the_matrix_correlator_data = h5py.File(f'{the_location}/Matrix_correlators_{the_type_rs}{the_re_bin}_{the_version}.h5','w')    
    
    begin_time = time.time()
    ### Start of the analysis for the nr. of irreps.
    for the_irrep in m_irreps:
         ### The list of operators 
        the_op_list = list(the_archivo[the_irrep].attrs['op_list'])
        
        ### The size of the matrix
        the_size_matrix = len(the_op_list)
        
        ### Extracting the original/raw data for the analysis
        the_datos_raw_pre = np.asarray(the_archivo[f'{the_irrep}/data'])[:the_number_cnfgs]
        
        ### Making all matrices Hermitian
        the_datos_raw = 0.5 * (the_datos_raw_pre + the_datos_raw_pre.conj().transpose(0,2,1,3))
        
        ### Time slices
        the_times = str(the_archivo[the_irrep].attrs['Other_Info']).split(' \n ')
        
        ### Min time slice
        the_min_nt = int(the_times[0][the_times[0].index('= ') + 2:])
        
        ### Max time slice
        the_max_nt = int(the_times[1][the_times[1].index('= ') + 2:])
        
         ### Total range
        the_nt = np.arange(the_min_nt, the_max_nt + 1)
        
        print('\n----------------------------------------------')
        print(f' -->   IRREP ({the_irreps.index(the_irrep)+1}/{len(the_irreps)}): {the_irrep}')
            
        the_datos = np.empty((the_size_matrix, the_size_matrix, the_datos_raw.shape[-1], len(binned_rw)), dtype=the_datos_raw.dtype)
        
        for n1, n2 in np.ndindex(the_size_matrix, the_size_matrix):
            
            ### Keeping the data complex
            the_n_data = the_datos_raw[:, n1, n2, :]
            rw_datos = vfa.REWEIGHTED_CORR(the_n_data, the_weight)
            if kwargs.get('rebin_on'):
                if rb == 'C':
                    rw_t_datos, bin_weights = np.asarray([vfa.BINNING_CUSTOM(rw_datos[tt], rb) for tt in range(rw_datos.shape[0])])
                    bins = list(rw_t_datos)
                else:
                    rw_t_datos = np.asarray([vfa.BINNING(rw_datos[tt], rb) for tt in range(rw_datos.shape[0])])
                    bins = None
                the_datos[n1, n2, :, :] = rw_t_datos
            else:
                bins = None
                the_datos[n1, n2, :, :] = rw_datos
        
        ### Resampling must be done with the normalized reweighting factors
        if the_type_rs=='jk':
            the_rs = np.empty((the_size_matrix, the_size_matrix, the_datos_raw.shape[-1], len(binned_rw)), dtype=the_datos.dtype)
            for n1, n2, tt in np.ndindex(the_datos.shape[:-1]):
                bin_rw = bins if bins is not None else None
                the_rs[n1, n2, tt, :] = jackknife_choice(the_datos[n1, n2, tt, :], norm_reweight, bin_rw)
        elif the_type_rs=='bt':
            the_rs = np.empty((the_size_matrix, the_size_matrix, the_datos_raw.shape[-1], len(bt_cfgs)), dtype=the_datos.dtype)
            for n1, n2, tt in np.ndindex(the_datos.shape[:-1]):
                the_rs[n1, n2, tt, :] = vfa.BOOTSTRAP(the_datos[n1, n2, tt, :], bt_cfgs, norm_reweight)
                

        ### Information about the ongoing analysis
        print('----------------------------------------------\n               DATA SHAPE \n----------------------------------------------')
        print(f'Nr. of gauge configurations: {the_datos.shape[-1]}\nSize of the Correlation matrix: {the_size_matrix}x{the_size_matrix}\nTime slices: {the_nt[0]} to {the_nt[-1]}\nResampling data ({the_resampling_scheme}): {the_rs.shape[-1]}\n....................................')
        
        print('      OPERATORS LIST ')
        for i in range(the_size_matrix):
            print(f'      -->>  {the_op_list[i]}')
                
        g_i = the_matrix_correlator_data.create_group(the_irrep) 
        g_i.create_dataset('Time_slices', data = the_nt)
        g_i.create_dataset('Operators', data = the_op_list) 
        group_corr = g_i.create_group('Correlators')
        group_corr_real = group_corr.create_group('Real')
        
        ### Calculating the mean values of the datasets
        the_mrs_f, the_mrs_f_rs = np.empty((the_datos.shape[:-1]), dtype=the_datos.dtype), np.empty((the_datos.shape[:-1]), dtype=the_datos.dtype)
        for n1, n2 in np.ndindex(the_size_matrix, the_size_matrix):
            if rb == 'C':
                the_mrs_f[n1, n2, :] = vfa.MEAN_CUSTOM(the_datos[n1, n2, :, :], bins)
            else:
                the_mrs_f[n1, n2, :] = vfa.MEAN(the_datos[n1, n2, :, :])
            the_mrs_f_rs[n1, n2, :] = vfa.MEAN(the_rs[n1, n2, :, :])
        
        ### The statistical error of the resampled data        
        the_sigmas_corr = np.empty((the_size_matrix, the_datos_raw.shape[-1]))
        for ss in range(the_size_matrix):
            the_sigmas_corr[ss,:] = vfa.STD_DEV_MEAN(the_rs[ss,ss,:,:], the_mrs_f_rs[ss,ss,:], the_type_rs)
        
        ### Reshaping the data for later diagonizalization and extraction of eigenvalues
        re_rs = vfa.RESHAPING_CORRELATORS_RS_NT(the_rs)
        re_mean = vfa.RESHAPING_EIGENVALS_MEAN(the_mrs_f)
        re_mean_rs = vfa.RESHAPING_EIGENVALS_MEAN(the_mrs_f_rs)
        
        group_corr_real.create_dataset('Resampled',data = re_rs)
        group_corr_real.create_dataset('Mean', data = re_mean)
        group_corr_real.create_dataset('Sigmas', data = the_sigmas_corr)
    the_matrix_correlator_data.close()
    end_time = time.time()
    print(f'TIME TAKEN: {(end_time-begin_time)/60} mins')    
    return f'{the_location}/Matrix_correlators_{the_type_rs}{the_re_bin}_{the_version}.h5'



def MultiCorrelatorAnalysisRatios(the_multi_hadrons_archivo, the_single_hadrons_archivo, the_location, the_list_non_interacting_hads, the_version, the_type_rs, the_m_irreps, the_s_irreps, **kwargs):
    
    print("                     CORRELATORS RATIO ANALYSIS \n")
    
     ### It chooses the rebin
    if kwargs.get('rebin_on'): 
        if kwargs.get('rb')==None:
            rb, the_re_bin = 1, '' 
            print(f'Missing argument: bin size. Using bin size equals {rb}.')
        else:
            rb = int(kwargs.get('rb'))
            the_re_bin = f'_bin{rb}'
    else:
        rb, the_re_bin = 1, ''  
        
    ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(the_m_irreps)
    
    m_irreps = the_m_irreps[the_first_irrep:the_last_irrep]
    
    ### This is the single correlators data
    the_ratio_matrix_correlator_data = h5py.File(f'{the_location}/Ratio_matrix_correlators_{the_type_rs}{the_re_bin}_{the_version}.h5','w')    
    
    begin_time = time.time()
    ### Start of the analysis for the nr. of irreps.
    for the_irrep in m_irreps:
        
        ### The matrix correlators
        this_data = the_multi_hadrons_archivo[the_irrep]
        
        ### The operators list in this irrep and the time slices
        the_op_list, the_nt = list(this_data['Operators']), np.array(this_data['Time_slices'])
        
        ### The number of operators in this irrep
        the_size_matrix = len(the_op_list)
        
        ### Printing some information in screen
        print('\n--------------------------------------------------------------------------------------')
        print(f'-->  IRREP({the_m_irreps.index(the_irrep)+1}/{len(the_m_irreps)}): {the_irrep}')
        print(f'     Size Matrix: {the_size_matrix}x{the_size_matrix}')
        print(f'     OPERATORS LIST \n----------------------------------------------')
        for ii in range(the_size_matrix):
            print(f'       {the_op_list[ii].decode("utf-8")}')

        the_list = list(the_list_non_interacting_hads)[the_m_irreps.index(the_irrep)]
        
        irreps_group = the_ratio_matrix_correlator_data.create_group(the_irrep)
        irreps_group.create_dataset('Time_slices', data = the_nt)
        irreps_group.create_dataset('Operators', data = the_op_list)
        irreps_group.create_dataset('Single_hadron_corrs', data = the_list)
        
        ### Getting all the t0's 
        try: 
            the_gevp_group = this_data.get('GEVP')
        except KeyError: 
            sys.exit("No GEVP done. No Ratios can be computed. ")
        
        the_nr_non_int = len(the_list)
        the_denom_mean, the_denom_rs = [], []
        
        ### Getting the non-int levels data (just once)
        for non_int in range(the_nr_non_int):
            ## Getting the Non-interacting levels
            the_non_int_states = vfa.NonInteractingLevels(the_list[non_int], the_s_irreps)
            the_first_non, the_second_non = the_non_int_states.First, the_non_int_states.Second
            
            ### Single-hadron correlators
            the_first_mean  = np.asarray(the_single_hadrons_archivo[f'{the_first_non}/Correlators/Real/Mean']).real
            the_first_rs = np.asarray(the_single_hadrons_archivo[f'{the_first_non}/Correlators/Real/Resampled']).real

            the_second_mean = np.asarray(the_single_hadrons_archivo[f'{the_second_non}/Correlators/Real/Mean']).real
            the_second_rs = np.asarray(the_single_hadrons_archivo[f'{the_second_non}/Correlators/Real/Resampled']).real
            
            ### Compute denominators
            the_denom_mean.append(the_first_mean * the_second_mean)
            the_denom_rs.append(the_first_rs * the_second_rs)  # shape (Nt, Ncfg)
            
        the_denom_mean = np.asarray(the_denom_mean)
        the_denom_rs = np.asarray(the_denom_rs)
            
        ## Loop over all the possible t0s
        for item in the_gevp_group:
            
            the_t0_group = irreps_group.create_group(f'GEVP/{item}')
            the_t0_group_mean = the_t0_group.create_group('Eigenvalues')
            
            the_eigenvals = this_data[f'GEVP/{item}/Eigenvalues']
            
            ### Getting the mean values and the resample values (Eigenvalues)
            the_eignvals_mean = np.asarray(the_eigenvals.get('Mean'))
            the_eignvals_rs = np.asarray(the_eigenvals.get('Resampled'))
            
            the_nr_eigs = len(the_eignvals_mean)         
            
            ### Creating the new arrays that will be saved in the ratio files
            the_ratio_eigenvalues_mean, the_ratio_eigenvalues_rs, the_ratio_cov_matrix = [], [], []
            
            ### Loop over all the eigenvalues
            for the_eign in range(the_nr_eigs):
                # lists of shape (n_non_int, nt)
                ratios_mean_i, ratios_rs_i, covs_i = [], [], []
                
                ### Transform eigenvalue resamples
                dis_rs_eign = vfa.NT_TO_NCFGS(the_eignvals_rs[the_eign])

                ### Loop over all the possible SH combinations
                for non_int in range(the_nr_non_int):

                    ### Ratio (mean)
                    the_ratio_mean = the_eignvals_mean[the_eign] / the_denom_mean[non_int]
                    
                    ### Ratio (resampled)
                    the_ratio_rs = dis_rs_eign / the_denom_rs[non_int]
                    
                    ### Covariance matrix
                    cov_i = vfa.COV_MATRIX(the_ratio_rs, the_ratio_mean, the_type_rs)
                    ratios_mean_i.append(the_ratio_mean)
                    ratios_rs_i.append(the_ratio_rs)
                    covs_i.append(cov_i)
                    
                the_ratio_eigenvalues_mean.append(np.asarray(ratios_mean_i))   # shape (n_non_int, nt)
                the_ratio_eigenvalues_rs.append(np.asarray(ratios_rs_i))
                the_ratio_cov_matrix.append(np.asarray(covs_i))
            the_t0_group_mean.create_dataset('Mean', data = np.asarray(the_ratio_eigenvalues_mean))
            the_t0_group_mean.create_dataset('Resampled', data = np.asarray(the_ratio_eigenvalues_rs))
            the_t0_group_mean.create_dataset('Covariance_matrix', data = np.asarray(the_ratio_cov_matrix))
                    
    the_ratio_matrix_correlator_data.close()
    end_time = time.time()
    print(f'TIME TAKEN: {(end_time-begin_time)/60} mins')    
    return f'{the_location}/Ratio_matrix_correlators_{the_type_rs}{the_re_bin}_{the_version}.h5'
    

### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------



### ------------------------------- START EXECUTING --------------------------------------------------


if __name__== "__main__":
    print("Nothing to do here.")
