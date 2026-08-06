import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.linalg import eig
from scipy.linalg import fractional_matrix_power
from iminuit import Minuit
from scipy.optimize import linear_sum_assignment
import re
import ast
import h5py
from collections import defaultdict

## ------------------- SOME USEFUL FUNCTIONS -----------------------------------
# Comments
# k_size: is the size of the bootstrap sample
# cfgs_nr: is the amount of gauge configs after rebinning
# it returns a list of list of numbers to plug into the boostrap sampling
def RANDOM_GENERATOR(k_size, cfgs_nr):
    np.random.seed(414557)
    return np.random.randint(0, cfgs_nr, size=(k_size, cfgs_nr))


# It receives a list [Ncfgs]
# It returns the normalization factor
def NORM_FACTOR(a_list):
    return np.double(np.mean(a_list, dtype=np.float128))


# It checks for hermiticity
# a_matrix: is the matrix [N,N]
def MAKES_HERMITIAN(a_matrix):
    return np.float128(1/2)*(np.matrix(a_matrix)+np.conj(np.matrix(a_matrix).T))   


def LINEAR_COMBINATION(the_list_of_operators, the_list_of_coeff, the_prefactor):
    the_new_correlator = np.float128(0.)
    ## OpLambda=(sqrt(3)/2)*(-Op27 +Op28 -Op33 -2*Op36 + Op37)
    for ii in range(len(the_list_of_operators)):
        the_new_correlator+=np.float128(the_prefactor*the_list_of_operators[ii]*the_list_of_coeff[ii])
    return the_new_correlator

### Comments:
# This receives a list d_list that is a file with the tmax for the fitting
# d_list: list of strings separated by '\n' and returns a list of integers
def T_MAX_LIST(d_list):
  d_list_0 = d_list[0]
  tmax_list = []
  for hh in range(len(d_list_0)):
      temp_list = d_list_0[hh].split(' ')
      other_temp = []
      for jj in range(len(temp_list)):
          if temp_list[jj]!='\n':
              other_temp.append(int(temp_list[jj]))
      tmax_list.append(other_temp)
  return tmax_list

### Comments:
# a_string: this is a string of the name of the correlator data without being analysis. It gets from the string if the correlator was modified. It returns the list (First Ncfg, Last Ncfg, step Ncfgs)
def GETTING_MIN_MAX_CONFIGS(a_string):
    the_short_string = a_string[:-5]
    the_split_string = list(the_short_string.split('_'))
    the_important_part = the_split_string[-1]
    the_new_string = list(the_important_part.split('-'))
    return the_new_string


## ------------------- RESHAPING CORRELATORS -----------------------------------

#This function is meant to reshape the multihadron correlators from [Ncfgs, N,N,nt] -> [N,N, Ncfgs,nt]. It receives:
# a: The original data with that shape
# s: the amount of operators to construct the NxN "matrix".
def RESHAPING(a):
    return np.asarray(np.transpose(a, (1,2,0,3)))

# This function reshapes the correlators as: [nt, N, N] --> [N,N,nt]
# a: the list to reshape
# s: size of the matrix
def RESHAPING_CORRELATORS(a):
    return np.asarray(np.transpose(a, (1,2,0)))

# This function reshapes the correlators as: [nt, Ncfgs, N, N] --> [N, N, Ncgfs, nt]
# a: the list to reshape
# s: size of the matrix
def RESHAPING_CORRELATORS_RS(a):
    return np.asarray(np.transpose(a, (2,3,1,0)))

# This function reshapes the correlators as: [nt, Ncfgs, N, N] --> [N, N, nt, Ncgfs] 
# a: the list to reshape
# s: size of the matrix
def RESHAPING_CORRELATORS_RS_NT(a):
    return np.asarray(np.transpose(a, (2,3,0,1)))

# This function reshapes the data in the way: [N, N, nt] -> [nt,N,N] This is used to later get the eigenvalues easier form the matrices.
# a: list with the data
def RESHAPING_EIGENVALS_MEAN(a):
    return np.asarray(np.transpose(a, (2,0,1)))


# This function reshapes the eigenvalues from [nt, Ncfgs, Neigens] -->> [Ncfgs, nt, Neigens]. This is used to obatin the fits easier later. 
# a: is the list of eigenvals, n configs, and time slices.
def RESHAPING_EIGEN_FOR_SORTING(a):
    return np.asarray(np.transpose(a, (1,0,2)))

# This function reshapes the eigenvalues from [nt, Ncfgs, N, N] -->> [Ncfgs, nt, N,N]. This is used to obatin the fits easier later. 
# a: is the list of eigenvals, n configs, and time slices.
def RESHAPING_EIGENVEC_FOR_SORTING(a):
    return np.asarray(np.transpose(a, (1,0,2,3)))

# This function reshapes the eigenvalues from [Neigens, Ncfgs, nt] -->> [nt, Ncfgs, Neigens] 
# a: is the list of eigenvals, n configs, and time slices.
def RESHAPING_EIGENVALS_NN(a):
    return np.asarray(np.transpose(a, (2,1,0)))

# This function reshapes the eigenvalues from [nt, Ncfgs] -->> [Ncfgs, nt]. This is used to obatin the fits easier later. 
# a: is the list of n configs, and time slices.
# s: number of data entries
def NT_TO_NCFGS(a):
    return np.asarray(np.transpose(a, (1,0)))

# This function reshapes the eigenvalues from  [Ncfgs, nt] -->> [nt, Ncfgs]. This is used to obatin the fits easier later. 
# a: is the list of n configs, and time slices.
# s: number of data entries
def NCFGS_TO_NT(a):
    eig_corrs = []
    for t_slices in range(len(a[0])):
        corr_n=[]
        for n1 in range(len(a)):
            corr_n.append(a[n1][t_slices])
        eig_corrs.append(np.array(corr_n))
    return np.array(eig_corrs)

#Comments
# This function below shrink the covariance matrix to the size of the nt desired. if originally is: Nt_totalxNt_total, then now it will be: (Nt_total-Nt_min-Nt_max)x(Nt_total-Nt_min-Nt_max).
def SHRINK_MATRIX(c,low,up):
    return c[low:up, low:up].astype(np.double)


## ------------------- OPERATOR REMOVAL FUNCTIONS ----------------------------

#Comments:
# This functions removes cols and rows from  a correlation matrix, leaving it squared as it should be. 
# c: this is the mean value of correlator matrix. It has a shape of [N, N, nt]
# r: This is the resampled correlator matrix. It has the shape [N, N, nt, Ncfgs]
# ss: this is the index of the row/column to be removed
def REMOVE_ROWS_COLS(c,r,ss):
    c = np.asarray(c)
    r = np.asarray(r)    
    the_new_c = np.delete(np.delete(c, ss, axis=0), ss, axis=1)
    the_new_r = np.delete(np.delete(r, ss, axis=0), ss, axis=1)
    return the_new_c, the_new_r

#Comments:
# This functions adds cols and rows from  a correlation matrix, leaving it squared as it should be. 
# c: this is the mean value of correlator matrix. It has a shape of [N, N, nt]
# r: This is the resampled correlator matrix. It has the shape [N, N, nt, Ncfgs]
# ss: this is the index of the row/column to be added
def ADD_ROWS_COLS(c,r,ss):
    c = np.asarray(c)
    r = np.asarray(r)
    return c[:ss+1, :ss+1], r[:ss+1, :ss+1]

#Comments:
# This functions adds cols and rows from  a correlation matrix, leaving it squared as it should be. 
# c: this is the mean value of correlator matrix. It has a shape of [N, N, nt]
# r: This is the resampled correlator matrix. It has the shape [N, N, nt, Ncfgs]
# ss: This is a list of operators to be included
def CHOOSE_OPS(c,r,ss):
    c = np.asarray(c)
    r = np.asarray(r)
    the_index = np.ix_(ss,ss)
    return c[the_index], r[the_index]


#################### ----- ZOE -----

### Comments:
# Splits the input from the ensemble dictionary 'operatorRemovalChoice' into the strings and the irrep name
def INPUT_SPLIT(a):
    return [block[1:] for block in a], [block[0] for block in a]

### Comments:
# Turns the string given from input split into a list of indices.
def STRING_TO_NUMBERS(a):
    b = [] # indices we want to remove
    c = [] # indices we want to keep
    for i in range(len(a)):
        if int(a[i]) == 0:
            b.append(i)
        else:
            c.append(i)
    return b,c

### Comments:
# copies the attributes to a given group
def COPY_ATTRIBUTES(group_name, attr1, attr2):
    group_name.attrs['Other_Info'] = attr1
    group_name.attrs['op_list'] = attr2

### Comments:
# matrix has the shape (ncfgs, nr ops, nr ops, nt)
# removes the operators according to the list indices_to_remove
def REMOVE_OPERATORS(indices_to_remove, matrix):
    matrix_new = matrix.copy()
    for i in reversed([1,2]):
        matrix_new = np.delete(matrix_new, indices_to_remove, axis=i)
    return matrix_new

### Comments:
# turns the input into a list
def IRREP_OP_LIST(the_Nr_Irreps, the_Nr_Ops):
    the_Nr_Irreps_return = ast.literal_eval(the_Nr_Irreps)
    the_Nr_Ops_return = ast.literal_eval(the_Nr_Ops)
    return the_Nr_Irreps_return, the_Nr_Ops_return


## ------------------- EFFECTIVE MASSES -----------------------------------

# This function receives a list "a" of time slices, and it calculates the effective mass, returning a list of effectives masses for each time slice (half integer numbers).
# a: shape [nt]
# This is the old version of the Effective Masses. It can still be used, but the other one is more general
def EFF_MASS_CLASSIC(a):
    meff = np.asarray([np.log(np.double(a[i])/np.double(a[i+1])) for i in range(len(a)-1)])
    return meff 

# This function receives a list "a" of time slices, and it calculates the effective mass, returning a list of effectives masses for each time slice (half integer numbers).
# a: shape [nt]
# d: distance of two points, by default this is 1
def EFF_MASS(a,d):
    meff = np.asarray([np.log(np.abs(np.double(a[i])/np.double(a[i+d]))) for i in range(len(a)-d)])
    return meff


# This function receives a list "a" of time slices, and it calculates the effective mass, returning a list of effectives masses for each time slice (half integer numbers).
# a: shape [nt]
# d: distance of two points, by default this is 1
def EFF_MASS_COSH(a,d):
    meff=[np.acosh(np.abs((np.double(a[i+d]) + np.double(a[i-d]))/np.double(2.*a[i]))) for i in range(len(a)-d)]
    return np.asarray(meff)


def DOING_EFFECTIVE_MASSES_EIGENVALUES(gevp_group, the_dist_eff_mass, the_type_rs):
    for item in gevp_group.keys():
        the_group_item = gevp_group[item]
        if 'Effective_masses' in the_group_item.keys(): 
            del gevp_group[f'{item}/Effective_masses']
        group_em_t0 = the_group_item.create_group('Effective_masses')
        
        the_evalues_rs_f = np.asarray(the_group_item['Eigenvalues/Resampled'])
        the_evalues_mean_f = np.asarray(the_group_item['Eigenvalues/Mean'])
        
        the_n_modes = the_evalues_mean_f.shape[0]
        ### Loop over the total number of eigenvalues
        the_eff_mass_mean = np.empty((the_n_modes,the_evalues_mean_f.shape[1]-1))
        the_cov_eff_mass = np.empty((the_n_modes,the_evalues_mean_f.shape[1]-1))
        
        for ls in range(the_n_modes):
            the_eff_mass_mean[ls] = EFF_MASS(the_evalues_mean_f[ls,:], the_dist_eff_mass)
            the_eff_mass_rs = np.array([EFF_MASS(the_evalues_rs_f[ls,zz,:], the_dist_eff_mass) for zz in range(the_evalues_rs_f.shape[1])])
            ### Reshaping the data
            the_eff_mass_rs = NT_TO_NCFGS(the_eff_mass_rs)
            
            ### Here the statistical errors of the resampled data are computed
            the_eff_rs_mean = MEAN(np.asarray(the_eff_mass_rs))
            the_cov_eff_mass[ls] = STD_DEV_MEAN(the_eff_mass_rs, the_eff_rs_mean, the_type_rs)
        
        group_em_t0.create_dataset('Mean', data = np.asarray(the_eff_mass_mean))
        group_em_t0.create_dataset('Sigmas',data = np.asarray(the_cov_eff_mass))


## ------------------- GEVP --------------------------------------

### Comments:
# This function solves the GEVP using the modified regular eigenvalue problem, since it is more stable:
# the_ct0_mean: This is the reference Correlation matrix [N, N]
# the_mean_corr: This is the correlation matrix at a certain time slice t [N, N]
# It uses the definition: C^{-1/2}(t_{0}) C(t) C^{-1/2}(t_{0})^{\dagger} \vec{v}(t)= \lambda(t) \vec{v}(t), where the eigenvalues are the same than the ones obtained from the GEVP directly, but the eigenvectors are the correct ones associated to physical states. 
def SOLVING_GEVP(the_ct0_mean, the_mean_corr):
    ### This matrix is to compute the modified GEVP
    the_ct0_root = fractional_matrix_power(the_ct0_mean, -1/2)
    the_mean_corr_mod = the_ct0_root @ the_mean_corr @ (the_ct0_root.conj().T)
    # return eigh(the_mean_corr_mod, eigvals_only=False) # HERMITIAN VERSION
    return eig(the_mean_corr_mod) # ANY MATRIX


def DOING_THE_GEVP_SINGLE_PIVOT(the_t0_min_max, the_nt, the_mean_corr, the_rs_real, the_type_rs, the_sorting, the_sorting_process, group_i, **kwargs):
    print("SINGLE PIVOT")
    if kwargs.get('t_diag')==None:
        the_td=int(len(the_nt)/3)
    else:
        the_td=int(kwargs.get('t_diag'))
    the_t0_init=0
    for the_t0_init in range(np.abs(the_t0_min_max[0] - the_nt[0]), (the_t0_min_max[1] - the_nt[0]) + 1): 
        
        ### This is the reference correlation matrix for the GEVP
        the_ct0_mean = np.array(the_mean_corr[the_t0_init])
    
        the_evals_mean, the_evecs_mean = [], []
        the_evalues_rs, the_evectors_rs =[], []
        
        ### Diagonalization only at one time slice the_td
        the_evs_mean_nongevp, the_evec_mean_nongevp = SOLVING_GEVP(the_ct0_mean, the_mean_corr[the_td])
        
        ### Diagonalization using the eigenvectors from solving the GEVP at only one time slice
        for ttt in range(len(the_mean_corr)):
            # if ttt!=the_td:
            the_evs_mean_nongevp_t = np.real(np.diag(np.linalg.multi_dot([the_evec_mean_nongevp.conj().T, the_mean_corr[ttt], the_evec_mean_nongevp] )))
            
            the_evals_mean.append(the_evs_mean_nongevp_t)

            the_evecs_mean.append(the_evec_mean_nongevp)
            
            ### Loop over the resamples
            the_evalues_rs_raw, the_evectors_rs_raw = [], []
            xyz = 0
            for xyz in range(the_rs_real.shape[1]):
                the_ew_rs_nongevp = np.real(np.diag(np.linalg.multi_dot([ the_evec_mean_nongevp.conj().T, the_rs_real[ttt][xyz], the_evec_mean_nongevp ])))
                
                the_evalues_rs_raw.append(np.array(the_ew_rs_nongevp))
                the_evectors_rs_raw.append(np.array(the_evec_mean_nongevp,dtype=np.float128))              
                
            the_evalues_rs.append(np.array(the_evalues_rs_raw))
            the_evectors_rs.append(np.array(the_evectors_rs_raw))
                

        the_evals_mean, the_evecs_mean = SORTING_EIGENVALUES(the_t0_init, the_evals_mean, the_evecs_mean)       
        
        ### Here the final eigenvalues and eigenvectors are saved
        if len(the_evalues_rs)>0:
            
            ### Reshaping eigenvectors and eigenvalues for sorting
            the_mod_evals_rs = RESHAPING_EIGEN_FOR_SORTING(np.array(the_evalues_rs))
            the_mod_evectors_rs = RESHAPING_EIGENVEC_FOR_SORTING(np.array(the_evectors_rs))
                
                ### Reshaping again to save them in a file
            the_evalues_rs = RESHAPING_EIGEN_FOR_SORTING(the_mod_evals_rs)
            the_evectors_rs = RESHAPING_EIGENVEC_FOR_SORTING(the_mod_evectors_rs)                
            
            group_t0 = group_i.create_group(f't0_{the_t0_init+the_nt[0]}')
            
            the_eigevals_final_mean = NT_TO_NCFGS(the_evals_mean)
            the_evals_fits_rs = np.array(RESHAPING_EIGENVALS_NN(np.array(the_evalues_rs)), dtype=np.float128)

            ### Getting the statistical error and the covariance matrix for each eigenvalue.
            the_l, the_sigma_2 = 0, []
            for the_l in range(len(the_evals_fits_rs)):
                dis_eign = NCFGS_TO_NT(the_evals_fits_rs[the_l])
                the_evals_fits_rs_mean = MEAN(dis_eign)
                the_sigma_2.append(COV_MATRIX(dis_eign, the_evals_fits_rs_mean, the_type_rs))
            
            ### Modified-GEVP eigenvectors (central values)
            the_evecs_mean = np.array(the_evecs_mean)
            
            group_eigvecs = group_t0.create_group('Eigenvectors')
            group_eigvecs.create_dataset('Mean', data=the_evecs_mean)
            group_eigvecs.create_dataset('Resampled', data=the_evectors_rs)
            
            group_eigns = group_t0.create_group('Eigenvalues')
            group_eigns.create_dataset('Mean', data = the_eigevals_final_mean)
            group_eigns.create_dataset('Resampled', data = the_evals_fits_rs)
            group_eigns.create_dataset('Covariance_matrix', data = np.array(the_sigma_2))
            print(f'T0 = {the_t0_init + the_nt[0]}...DONE')




### Comments:
# This function does the GEVP for each time slice for each chosen t0 and for each resample. Sorting everything by eigenvalue by default. If a different sorting is required, then it is done separately.
# the_t0_min_max: is a list with (t0 min, t0 max) to do the GEVP
# the_nt: range where the GEVP will be solved
# the_mean_corr: this is the data to be used for the GEVP, shape = [nt, N, N]
# the_rs_real: This is the resmapled data to be used for the GEVP: shape = [nt, Ncnfgs, N, N]
# the_type_rs: this is the resample type for the covariance matrix computation
# group_i: This is the "folder" in which the information will be saved in the hdf5 file.
# the_sorting: This is the method chosen for sorting
# the_sorting_process: This is the function that will be called with the sorting.
def DOING_THE_GEVP(the_t0_min_max, the_nt, the_mean_corr, the_rs_real, the_type_rs, the_sorting, the_sorting_process, the_rs_sorting_process, group_i):
    
    print("ROLLING PIVOT")
    
    ### Choosing the t0's for which to perform the GEVP
    the_t0_start = np.abs(the_t0_min_max[0] - the_nt[0])
    the_t0_end = (the_t0_min_max[1] - the_nt[0]) + 1

    the_time_slices = the_mean_corr.shape[0]
    the_ncnfgs = the_rs_real.shape[1]
    
    ### Loop over the t0's chosen
    for the_t0_init in range(the_t0_start, the_t0_end):
        ### This is the reference correlation matrix for the GEVP
        the_ct0_mean = the_mean_corr[the_t0_init]        
        print(f"Checking Hermiticity: {np.allclose(the_ct0_mean, the_ct0_mean.conj().T, atol=1e-12)}")

        ### ---------- MEAN VALUES ----------
        ### Loop over the time slices (diagonalization at each time slice)
        try:
            the_ev_mean = [SOLVING_GEVP(the_ct0_mean, the_mean_corr[ttt]) for ttt in range(the_time_slices)]
        except np.linalg.LinAlgError:
            print(f"WARNING: Matrix isn't positive definite anymore. Skipping T0 = {the_t0_init + the_nt[0]}")
            continue
        
        the_evals_mean = np.asarray([the_ev[0] for the_ev in the_ev_mean])
        the_evecs_mean = np.asarray([the_ev[1] for the_ev in the_ev_mean])
        
        ### First sort by eigenvalue
        the_evals_mean, the_evecs_mean = SORTING_EIGENVALUES(the_t0_init, the_evals_mean.real, the_evecs_mean)
        
        if the_sorting is not None and the_sorting!='eigenvals':
            the_evals_mean, the_evecs_mean = the_sorting_process(the_evals_mean, the_evecs_mean, the_t0_init)
        
        ### ---------- RESAMPLED VALUES ----------        
        the_evals_rs, the_evecs_rs = [], []
        the_ct0_rs = the_rs_real[the_t0_init]
        
        ### Loop over the time slices for the resamples
        for ttt in range(the_time_slices):
            try:
                the_ct_rs = the_rs_real[ttt]
                
                the_evals_rs_t , the_evecs_rs_t = [], []
                ### Loop over the resamples
                for xyz in range(the_ncnfgs):
                    the_ev_mean_rs, the_evec_mean_rs = SOLVING_GEVP(the_ct0_rs[xyz], the_ct_rs[xyz])
                    
                    the_evals_rs_t.append(np.asarray(the_ev_mean_rs))
                    the_evecs_rs_t.append(np.asarray(the_evec_mean_rs))
                
                the_evals_rs.append(np.asarray(the_evals_rs_t))
                the_evecs_rs.append(np.asarray(the_evecs_rs_t))
                
            except np.linalg.LinAlgError:
                break
        
        ### Here the final eigenvalues and eigenvectors are saved
        if len(the_evals_rs) == 0:
            continue
        
        ### Reshaping eigenvectors and eigenvalues for sorting
        the_mod_evals_rs = RESHAPING_EIGEN_FOR_SORTING(np.asarray(the_evals_rs))
        the_mod_evectors_rs = RESHAPING_EIGENVEC_FOR_SORTING(np.asarray(the_evecs_rs))
        
        ### Loop over the resamples (sorting)
        for xyz in range(len(the_mod_evals_rs)):
            the_mod_evals_rs[xyz], the_mod_evectors_rs[xyz] = SORTING_EIGENVALUES(the_t0_init, the_mod_evals_rs[xyz].real, the_mod_evectors_rs[xyz])

        the_eigevals_final_mean = NT_TO_NCFGS(the_evals_mean)
        the_evals_fits_rs = RESHAPING_EIGENVALS_MEAN(the_mod_evals_rs)

        ### Getting the statistical error and the covariance matrix for each eigenvalue.
        the_sigma_2 = [COV_MATRIX(dis_eign:= NT_TO_NCFGS(the_evals_fits_rs[the_l]), MEAN(dis_eign), the_type_rs) for the_l in range(len(the_evals_fits_rs))]

        ### Modified-GEVP eigenvectors (central values)
        the_evecs_mean = np.asarray(the_evecs_mean)
        
        ### Modified-GEVP eigenvectors (resamples)
        the_evecs_rs = np.asarray(the_mod_evectors_rs)        
        
        group_t0 = group_i.create_group(f't0_{the_t0_init + the_nt[0]}')

        group_eigvecs = group_t0.create_group('Eigenvectors')
        group_eigvecs.create_dataset('Mean', data = the_evecs_mean)
        group_eigvecs.create_dataset('Resampled', data = the_evecs_rs)

        group_eigns = group_t0.create_group('Eigenvalues')
        group_eigns.create_dataset('Mean', data = the_eigevals_final_mean)
        group_eigns.create_dataset('Resampled', data = the_evals_fits_rs)
        group_eigns.create_dataset('Covariance_matrix', data = np.asarray(the_sigma_2))
        print(f'T0 = {the_t0_init + the_nt[0]}...DONE')



### Comments:
# Choose T_0 so that it fulfills the constraints in Paper by Blossier et al.: t_0 >= t / 2 and t-t_0 = const
# Starts at T_0 = 3 in case of an extended source
def T0_RETURN(n_t):
    T_0 = int
    C = 1
    if n_t < 4:
        T_0 = 3
    else: T_0 = n_t - C
    return T_0

### Comments:
# Solves the GEVP but the time slice diagonalized to changes with n_t -> no loop over T0
# the_t0_sorting: time slice where the sorting should start, e.g. = 3
# everything else like in DOING_THE_GEVP
def DOING_THE_GEVP_RUNNING_T0(the_t0_sorting, the_nt, the_mean_corr, the_rs_real, the_type_rs, the_sorting, the_sorting_process, the_rs_sorting_process, group_i):
    print('WINDOW METHOD')

    the_evals_mean, the_evecs_mean = [], []
    the_evalues_rs, the_evectors_rs = [], []

    for ttt in range(len(the_mean_corr)):
        try:
            ### This matrix is to compute the modified GEVP
            the_ct0_mean = np.array(the_mean_corr[T0_RETURN(ttt)])
            the_evs_mean_nongevp, the_evec_mean_nongevp = SOLVING_GEVP(the_ct0_mean, the_mean_corr[ttt])

            the_evals_mean.append(the_evs_mean_nongevp)
            the_evecs_mean.append(the_evec_mean_nongevp)

        except np.linalg.LinAlgError:
            print(
                f"WARNING: Matrix isn't positive definite anymore. Skipping T0 = {T0_RETURN(ttt) + the_nt[0]}")
            break

    ### First sort by eigenvalue
    the_evals_mean, the_evecs_mean = SORTING_EIGENVALUES_NEW(the_evals_mean, the_evecs_mean, the_t0_sorting)
    ### Now if wanted, sorting by whatever else
    if the_sorting == 'vecs_var_rs_mean':
        the_evals_mean, the_evecs_mean = the_sorting_process(the_evals_mean, the_evecs_mean, the_t0_sorting)
    elif the_sorting != None or the_sorting != 'eigenvals':
        the_evals_mean, the_evecs_mean = the_sorting_process(the_evals_mean, the_evecs_mean, the_t0_sorting)

    ### Loop over the time slices for the resamples
    ttt = 0
    for ttt in range(len(the_mean_corr)):
        the_evalues_rs_raw, the_evectors_rs_raw = [], []
        xyz = 0
        try:
            ### Loop over the resamples
            for xyz in range(the_rs_real.shape[1]):
                the_ew_rs_nongevp, the_ev_rw_nongevp = SOLVING_GEVP(np.array(the_rs_real[T0_RETURN(ttt)][xyz]), the_rs_real[ttt][xyz])

                the_evalues_rs_raw.append(np.array(the_ew_rs_nongevp))
                the_evectors_rs_raw.append(np.array(the_ev_rw_nongevp, dtype=np.float128))

            the_evalues_rs.append(np.array(the_evalues_rs_raw))
            the_evectors_rs.append(np.array(the_evectors_rs_raw))

        except np.linalg.LinAlgError:
            break

    ### Here the final eigenvalues and eigenvectors are saved
    if len(the_evalues_rs) > 0:

        ### Reshaping eigenvectors and eigenvalues for sorting
        the_mod_evals_rs = RESHAPING_EIGEN_FOR_SORTING(np.array(the_evalues_rs))
        the_mod_evectors_rs = RESHAPING_EIGEN_FOR_SORTING(np.array(the_evectors_rs))

        ### Loop over the resamples (sorting)
        for xyz in range(len(the_mod_evals_rs)):
            the_mod_evals_rs[xyz], the_mod_evectors_rs[xyz] = SORTING_EIGENVALUES_NEW(the_mod_evals_rs[xyz], the_mod_evectors_rs[xyz])
            if the_rs_sorting_process == SORTING_EIGENVECTORS_RS_MEAN:
                the_mod_evals_rs[xyz], the_mod_evectors_rs[xyz] = SORTING_EIGENVECTORS_RS_MEAN(the_mod_evals_rs[xyz], the_mod_evectors_rs[xyz], the_t0_sorting, the_evecs_mean, rs=True)
            else:
                the_mod_evals_rs[xyz], the_mod_evectors_rs[xyz] = the_rs_sorting_process(the_mod_evals_rs[xyz], the_mod_evectors_rs[xyz], the_t0_sorting)

            ### Reshaping again to save them in a file
        the_evalues_rs = RESHAPING_EIGEN_FOR_SORTING(the_mod_evals_rs)
        the_evectors_rs = RESHAPING_EIGEN_FOR_SORTING(the_mod_evectors_rs)

        group_t0 = group_i.create_group(f't0_{T0_RETURN(1) + the_nt[0]}_run')

        the_eigevals_final_mean = NT_TO_NCFGS(the_evals_mean)
        the_evals_fits_rs = np.array(RESHAPING_EIGENVALS_NN(np.array(the_evalues_rs)), dtype=np.float128)

        ### Getting the statistical error and the covariance matrix for each eigenvalue.
        the_l, the_sigma_2 = 0, []
        for the_l in range(len(the_evals_fits_rs)):
            dis_eign = NCFGS_TO_NT(the_evals_fits_rs[the_l])
            the_evals_fits_rs_mean = MEAN(dis_eign)
            the_sigma_2.append(COV_MATRIX(dis_eign, the_evals_fits_rs_mean, the_type_rs))

        ### Modified-GEVP eigenvectors (central values)
        the_evecs_mean = np.array(the_evecs_mean)

        group_eigvecs = group_t0.create_group('Eigenvectors')
        group_eigvecs.create_dataset('Mean', data=the_evecs_mean)
        group_eigvecs.create_dataset('Resampled', data=the_evectors_rs)

        group_eigns = group_t0.create_group('Eigenvalues')
        group_eigns.create_dataset('Mean', data=the_eigevals_final_mean)
        group_eigns.create_dataset('Resampled', data=the_evals_fits_rs)
        group_eigns.create_dataset('Covariance_matrix', data=np.array(the_sigma_2))
        print(f'T0 WINDOW METHOD...DONE')
        
        
        
        

## ------------------- SORTING STATES --------------------------------------



def SORTING_PROCESS(the_sorting, ev = True):

    if the_sorting is None or the_sorting == 'eigenvals':
        if ev: print("Sorting states based on Eigenvalues.")
        ### This function returns the eigenvalues sorted from the largest to the smallest.
        the_sorting_process = SORTING_EIGENVALUES_NEW

    elif the_sorting == 'vecs_fix':
        if ev: print("Sorting states by Eigenvectors with a fixed reference time slice.")
        ### This function returns the eigenvalues sorted based on the orthogonality of the eigenvectors based on a reference time slice where the eigenstated are already sorted.
        the_sorting_process = SORTING_EIGENVECTORS

    elif the_sorting == 'vecs_fix_norm':
        if ev: print("Sorting states by normalized Eigenvectors with a fixed reference time slice.")
        ### This function returns the eigenvalues sorted based on the orthogonality of the normalized eigenvectors based on a reference time slice where the eigenstated are already sorted.
        the_sorting_process = SORTING_EIGENVECTORS_NORMALIZED

    elif the_sorting == 'vecs_var':
        if ev: print("Sorting states by Eigenvectors with a varying reference time slice.")
        ### This function returns the eigenvalues sorted based on the orthogonality of the eigenvectors based on the previous reference time slice.
        the_sorting_process = SORTING_EIGENVECTORS_CHANGING_TSLICE

    elif the_sorting == 'vecs_var_norm':
        if ev: print("Sorting states by normalized Eigenvectors with a varying reference time slice.")
        ### This function returns the eigenvalues sorted based on the orthogonality of the normalized eigenvectors based on the previous reference time slice.
        the_sorting_process = SORTING_EIGENVECTORS_NORMALIZED_CHANGING_TSLICE

    elif the_sorting == 'vecs_var_rs_mean':
        if ev: print("Sorting states by Eigenvectors with a varying reference time slice. The resamples are sorted in the same way as the mean.")
        ### This function returns the eigenvalues sorted based on the orthogonality of the eigenvectors based on the previous reference time slice. The order of the mean values are applied to the resamples.
        the_sorting_process = SORTING_EIGENVECTORS_RS_MEAN

    return the_sorting_process


# THIS SHOULD WORK
# def SORTING_EIGENVALUES(the_t0, the_eigenvals, the_eigenvecs):
#     the_final_eigens = list(the_eigenvals[:the_t0])
#     the_final_eigenvecs = list(the_eigenvecs[:the_t0])
#     for ii in range(the_t0, len(the_eigenvals)):
#         the_sorted_indices = sorted(range(len(the_eigenvals[ii])), key=lambda i: the_eigenvals[ii][i], reverse=True)
#         the_final_eigens.append(np.array([the_eigenvals[ii][i] for i in the_sorted_indices]))
#         the_final_eigenvecs.append(np.array([the_eigenvecs[ii][i] for i in the_sorted_indices]))
#     return [the_final_eigens, the_final_eigenvecs]


def SORTING_EIGENVALUES(the_t0, the_eigenvals, the_eigenvecs):
    the_final_eigens = np.array(the_eigenvals, copy=True)
    the_final_eigenvecs = np.array(the_eigenvecs, copy=True)
    ### Sorting according to the real part
    sorted_indices = np.argsort(the_eigenvals[:the_t0].real, axis=1)
    the_final_eigens[:the_t0] = np.take_along_axis(the_eigenvals[:the_t0], sorted_indices, axis=1)
    the_final_eigenvecs[:the_t0] = np.take_along_axis( the_eigenvecs[:the_t0], sorted_indices[:, :, np.newaxis], axis=1)
    sorted_indices = np.argsort(-the_eigenvals[the_t0:].real, axis=1)
    the_final_eigens[the_t0:] = np.take_along_axis(the_eigenvals[the_t0:], sorted_indices, axis=1)
    the_final_eigenvecs[the_t0:] = np.take_along_axis( the_eigenvecs[the_t0:], sorted_indices[:, :, np.newaxis], axis=1)
    return [the_final_eigens, the_final_eigenvecs]


def SORTING_EIGENVALUES_NEW(the_eigenvals, the_eigenvecs, the_t0 = 0):
    the_final_eigens = []
    the_final_eigenvecs = []
    for ii in range(len(the_eigenvals)):
        the_sorted_indices = sorted(range(len(the_eigenvals[ii])), key=lambda i: the_eigenvals[ii][i], reverse=True)
        the_final_eigens.append(np.asarray([the_eigenvals[ii][i] for i in the_sorted_indices]))
        the_final_eigenvecs.append(np.asarray([the_eigenvecs[ii][i] for i in the_sorted_indices]))
    return [the_final_eigens, the_final_eigenvecs]



### Comments: alternative sorting not relying on the sorting function but explicitely calculating the maximum
def SORTING_EIGENVALUES_ALTERNATIVE(the_t0, the_eigenvals, the_eigenvecs):
    the_final_eigens, the_final_eigenvecs, the_final_eigens_nt = list(the_eigenvals[:the_t0]), list(
        the_eigenvecs[:the_t0]), []
    for the_nt in range(the_t0, len(the_eigenvals)):
        the_eigenvals_nt, the_indices = list(the_eigenvals[the_nt]), []
        while len(the_eigenvals_nt) > 0:
            the_max_eigen = max(the_eigenvals_nt)
            the_final_eigens_nt.append(the_max_eigen)
            the_indices.append(list(the_eigenvals[the_nt]).index(the_max_eigen))
            the_eigenvals_nt.remove(the_max_eigen)
        the_final_eigens.append(the_final_eigens_nt)
        the_final_eigenvecs.append(np.array([the_eigenvecs[the_nt][nev] for nev in the_indices]))
        the_final_eigens_nt, the_indices = [], []
    return [the_final_eigens, the_final_eigenvecs]



### Comments: This function checks for orthogonality of eigenvectors and returns the new order with the eigenvectors ordered such that they are associated to the corresponding state time slice by time slice. No normalization of the vectors
# the_eigenvals: This is an array of the values of the eigenvals: shape [Nt, Neigens]
# the_eigenvecs: This is an array of the eigenvectors in the following shape: [Nt, N, N]
def SORTING_EIGENVECTORS(the_eigenvals, the_eigenvecs, the_t0):
    the_ref_tslice_eigenval = int((len(the_eigenvals) - the_t0) / 3)
    the_final_eigens, the_final_eigenvecs = SORTING_EIGENVALUES_NEW(the_eigenvals[:the_ref_tslice_eigenval + 1], the_eigenvecs[:the_ref_tslice_eigenval + 1])
    the_ref_eigenvec = the_eigenvecs[the_ref_tslice_eigenval]
    for ii in range(the_ref_tslice_eigenval + 1, len(the_eigenvals)):
        ckl_k = []
        for kk in range(len(the_ref_eigenvec)):
            ckl_l = []
            for ll in range(len(the_ref_eigenvec)):
                ckl_l.append(np.abs(np.dot(the_ref_eigenvec[kk], the_eigenvecs[ii][ll].T)))
            ckl_k.append(ckl_l)
        ckl_k = np.matrix(ckl_k)
        ckl = sorted([(ckl_k[i, j], i, j) for i in range(ckl_k.shape[0]) for j in range(ckl_k.shape[1])], reverse=True)
        the_used_rows, the_used_cols, the_selected = set(), set(), []
        the_top_n = None
        for the_value, the_i, the_j in ckl:
            if the_i not in the_used_rows and the_j not in the_used_cols:
                the_selected.append([the_i, the_j]);
                the_used_rows.add(the_i);
                the_used_cols.add(the_j)
                if the_top_n is not None and len(the_selected) >= the_top_n:
                    break
        the_sorted_indices = sorted(the_selected, key=lambda x: (x[0], x[1]))
        the_sorted_indices = [the_sorted_indices[i][1] for i in range(len(the_sorted_indices))][::-1] # reverse so that they also have the same order as the eigenvalue routine
        the_sorted_eigenvecs = np.array([the_eigenvecs[ii][i] for i in the_sorted_indices])
        the_final_eigens.append(np.array([the_eigenvals[ii][i] for i in the_sorted_indices]))
        the_final_eigenvecs.append(the_sorted_eigenvecs)
        
        
        

# Comments: This function checks for orthogonality of eigenvectors and returns the new order with the eigenvectors ordered such that they are associated to the corresponding state time slice by time slice. Each vector is normalized first.
# the_eigenvals: This is an array of the values of the eigenvals: shape [Nt, Neigens]
# the_eigenvecs: This is an array of the eigenvectors in the following shape: [Nt, N, N]
def SORTING_EIGENVECTORS_NORMALIZED(the_eigenvals, the_eigenvecs, the_t0):
    the_ref_tslice_eigenval = int((len(the_eigenvals) - the_t0) / 3)
    the_final_eigens, the_final_eigenvecs = list(the_eigenvals[:the_ref_tslice_eigenval + 1]), list(
        the_eigenvecs[:the_ref_tslice_eigenval + 1])
    the_ref_eigenvec = the_eigenvecs[the_ref_tslice_eigenval]
    for vec in range(len(the_ref_eigenvec)):
        the_ref_eigenvec[vec] = np.abs(the_ref_eigenvec[vec] / np.linalg.norm(the_ref_eigenvec[vec]))
    for ii in range(the_ref_tslice_eigenval + 1, len(the_eigenvals)):
        ckl_k = []
        for kk in range(len(the_ref_eigenvec)):
            ckl_l = []
            for ll in range(len(the_ref_eigenvec)):
                ckl_l.append(np.abs(
                    np.dot(the_ref_eigenvec[kk], the_eigenvecs[ii][ll] / np.linalg.norm(the_eigenvecs[ii][ll].T))))
            ckl_k.append(ckl_l)
        ckl_k = np.matrix(ckl_k)
        ckl = sorted([(ckl_k[i, j], i, j) for i in range(ckl_k.shape[0]) for j in range(ckl_k.shape[1])], reverse=True)
        the_used_rows, the_used_cols, the_selected = set(), set(), []
        the_top_n = None
        for the_value, the_i, the_j in ckl:
            if the_i not in the_used_rows and the_j not in the_used_cols:
                the_selected.append([the_i, the_j]);
                the_used_rows.add(the_i);
                the_used_cols.add(the_j)
                if the_top_n is not None and len(the_selected) >= the_top_n:
                    break
        the_sorted_indices = sorted(the_selected, key=lambda x: (x[0], x[1]))
        the_final_eigenvecs.append(the_sorted_eigenvecs)
        the_final_eigenvecs.append(the_sorted_eigenvecs)
        the_sorted_indices = [the_sorted_indices[i][1] for i in range(len(the_sorted_indices))][::-1] # reverse so that they also have the same order as the eigenvalue routine
        the_sorted_eigenvecs = np.array([the_eigenvecs[ii][i] for i in the_sorted_indices])
        the_final_eigens.append(np.array([the_eigenvals[ii][i] for i in the_sorted_indices]))
        the_final_eigenvecs.append(the_sorted_eigenvecs)
    return [the_final_eigens, the_final_eigenvecs]




# Comments: This function checks for orthogonality of eigenvectors and returns the new order with the eigenvectors ordered such that they are associated to the corresponding state time slice by time slice. No normalization of the vectors, and changing the reference time slice to compare with.
# the_eigenvals: This is an array of the values of the eigenvals: shape [Nt, Neigens]
# the_eigenvecs: This is an array of the eigenvectors in the following shape: [Nt, N, N]
def SORTING_EIGENVECTORS_CHANGING_TSLICE_OLD(the_eigenvals, the_eigenvecs, the_t0, rs = True):
    the_ref_tslice_eigenval = int((len(the_eigenvals) - the_t0) / 3)
    the_final_eigens, the_final_eigenvecs = SORTING_EIGENVALUES_NEW(the_eigenvals[:the_ref_tslice_eigenval + 1], the_eigenvecs[:the_ref_tslice_eigenval + 1])
    the_ref_eigenvec = the_eigenvecs[the_ref_tslice_eigenval][np.argmax(the_eigenvals[the_ref_tslice_eigenval])]
    all_sorted_indices = []
    for ii in range(the_ref_tslice_eigenval + 1, len(the_eigenvals)):
        ckl_k = []
        for kk in range(len(the_ref_eigenvec)):
            ckl_l = []
            for ll in range(len(the_ref_eigenvec)):
                ckl_l.append(np.abs(np.dot(the_ref_eigenvec[kk], the_eigenvecs[ii][ll].T)))
            ckl_k.append(ckl_l)
        ckl_k = np.matrix(ckl_k)
        ckl = sorted([(ckl_k[i, j], i, j) for i in range(ckl_k.shape[0]) for j in range(ckl_k.shape[1])], reverse=True)
        the_used_rows, the_used_cols, the_selected = set(), set(), []
        the_top_n = None
        for the_value, the_i, the_j in ckl:
            if the_i not in the_used_rows and the_j not in the_used_cols:
                the_selected.append([the_i, the_j])
                the_used_rows.add(the_i)
                the_used_cols.add(the_j)
                if the_top_n is not None and len(the_selected) >= the_top_n:
                    break
        the_sorted_indices = sorted(the_selected, key=lambda x: (x[0], x[1]))
        the_sorted_indices = [the_sorted_indices[i][1] for i in range(len(the_sorted_indices))] # reverse so that they also have the same order as the eigenvalue routine
        all_sorted_indices.append(the_sorted_indices)
        the_sorted_eigenvecs = np.array([the_eigenvecs[ii][i] for i in the_sorted_indices])
        the_final_eigens.append(np.array([the_eigenvals[ii][i] for i in the_sorted_indices]))
        the_final_eigenvecs.append(the_sorted_eigenvecs)
        the_ref_eigenvec = the_sorted_eigenvecs  # This is the part where the new sorted eigenvectors plays a role
    if not rs:
        return [the_final_eigens, the_final_eigenvecs], all_sorted_indices
    else:

        return [the_final_eigens , the_final_eigenvecs]





def SORTING_EIGENVECTORS_CHANGING_TSLICE(the_eigenvals, the_eigenvecs, the_t0):
    the_ref_tslice = int((len(the_eigenvals) - the_t0) / 3)

    sorted_indices_ref = np.argsort(the_eigenvals[the_ref_tslice])[::-1]  # größter Eigenwert zuerst
    the_ref_eigenvec = np.array([the_eigenvecs[the_ref_tslice][i] for i in sorted_indices_ref])

    the_final_eigens = [np.array([the_eigenvals[tt][i] for i in sorted_indices_ref])
                        for tt in range(the_ref_tslice + 1)]
    the_final_eigenvecs = [np.array([the_eigenvecs[tt][i] for i in sorted_indices_ref])
                           for tt in range(the_ref_tslice + 1)]

    all_sorted_indices = []

    for ii in range(the_ref_tslice + 1, len(the_eigenvals)):
        ckl_k = np.array([[np.abs(np.dot(the_ref_eigenvec[kk], the_eigenvecs[ii][ll].T))
                           for ll in range(len(the_ref_eigenvec))]
                          for kk in range(len(the_ref_eigenvec))], dtype=np.float128)

        ckl = sorted([(ckl_k[i, j], i, j)
                      for i in range(ckl_k.shape[0])
                      for j in range(ckl_k.shape[1])], reverse=True)

        the_used_rows, the_used_cols, the_selected = set(), set(), []

        for value, i_row, j_col in ckl:
            if i_row not in the_used_rows and j_col not in the_used_cols:
                the_selected.append([i_row, j_col])
                the_used_rows.add(i_row)
                the_used_cols.add(j_col)

        the_sorted_indices = [sel[1] for sel in sorted(the_selected, key=lambda x: x[0])]
        all_sorted_indices.append(the_sorted_indices)

        sorted_vecs = np.array([the_eigenvecs[ii][i] for i in the_sorted_indices])
        sorted_vecs = sorted_vecs / np.linalg.norm(sorted_vecs, axis=1, keepdims=True)

        for k in range(len(sorted_vecs)):
            if np.dot(the_ref_eigenvec[k], sorted_vecs[k]) < 0:
                sorted_vecs[k] *= -1

        the_final_eigenvecs.append(sorted_vecs)
        the_final_eigens.append(np.array([the_eigenvals[ii][i] for i in the_sorted_indices]))
        the_ref_eigenvec = sorted_vecs

    return [the_final_eigens, the_final_eigenvecs]



def SORTING_EIGENVECTORS_CHANGING_TSLICE_TEST(the_eigenvals, the_eigenvecs, the_t0):
    the_ref_tslice = int((len(the_eigenvals) - the_t0) / 3)

    sorted_indices_ref = np.argsort(the_eigenvals[the_ref_tslice])[::-1]
    the_ref_eigenvec = np.array([the_eigenvecs[the_ref_tslice][i] for i in sorted_indices_ref])
    
    the_final_eigens = [np.array([the_eigenvals[tt][i] for i in sorted_indices_ref])
                        for tt in range(the_ref_tslice + 1)]
    the_final_eigenvecs = [np.array([the_eigenvecs[tt][i] for i in sorted_indices_ref])
                           for tt in range(the_ref_tslice + 1)]

    for ii in range(the_ref_tslice + 1, len(the_eigenvals)):
        curr_vecs = np.array(the_eigenvecs[ii])

        the_ref_eigenvec = the_ref_eigenvec / np.linalg.norm(the_ref_eigenvec, axis=1, keepdims=True)
        curr_vecs = curr_vecs / np.linalg.norm(curr_vecs, axis=1, keepdims=True)

        for k in range(len(curr_vecs)):
            sign = np.sign(np.dot(the_ref_eigenvec[k], curr_vecs[k]))
            curr_vecs[k] *= sign

        overlap = np.abs(np.dot(the_ref_eigenvec, curr_vecs.T))

        row_ind, col_ind = linear_sum_assignment(-overlap)  # maximize overlap
        sorted_indices = col_ind

        sorted_vecs = curr_vecs[sorted_indices]
        sorted_vals = np.array([the_eigenvals[ii][i] for i in sorted_indices])

        the_final_eigens.append(sorted_vals)
        the_final_eigenvecs.append(sorted_vecs)

        the_ref_eigenvec = sorted_vecs

    return [the_final_eigens, the_final_eigenvecs]



### Comments:
# sorts the mean values by eigenvector like in SORTING_EIGENVECTOR_CHANGING_TSLICE if rs = False
# sorts the resamples by the same indices as the mean value for every time slice
def SORTING_EIGENVECTORS_RS_MEAN(the_eigenvals, the_eigenvecs, the_t0, mean_eigvec =None, rs=False):

    if not rs:
        results = SORTING_EIGENVECTORS_CHANGING_TSLICE(the_eigenvals, the_eigenvecs, the_t0)
        return results

    else:
        the_ref_tslice = int((len(the_eigenvals) - the_t0) / 3)
        the_final_eigens = list(the_eigenvals[:the_ref_tslice + 1])
        the_final_eigenvecs = list(the_eigenvecs[:the_ref_tslice + 1])
        for ii in range(the_ref_tslice +1, len(the_eigenvals)):
            ref = mean_eigvec[ii]
            curr = the_eigenvecs[ii]

            ref = ref / np.linalg.norm(ref, axis=1, keepdims=True)
            curr = curr / np.linalg.norm(curr, axis=1, keepdims=True)

            for i in range(len(curr)):
                sign = np.sign(np.dot(ref[i], curr[i]))
                curr[i] *= sign

            overlap = np.abs(np.dot(ref, curr.T))
            sorted_indices = list(range(len(curr)))

            for i in range(len(curr)):
                j_best = np.argmax(overlap[i])

                if overlap[i, j_best] > overlap[i, i]:
                    sorted_indices[i], sorted_indices[j_best] = sorted_indices[j_best], sorted_indices[i]

            sorted_eigenvecs = np.array([the_eigenvecs[ii][i] for i in sorted_indices])
            sorted_eigenvals = np.array([the_eigenvals[ii][i] for i in sorted_indices])

            the_final_eigens.append(sorted_eigenvals)
            the_final_eigenvecs.append(sorted_eigenvecs)
        return [the_final_eigens, the_final_eigenvecs]
    
    
    
    
    
# Comments: This function checks for orthogonality of eigenvectors and returns the new order with the eigenvectors ordered such that they are associated to the corresponding state time slice by time slice
# the_eigenvals: This is an array of the values of the eigenvals: shape [Nt, Neigens]
# the_eigenvecs: This is an array of the eigenvectors in the following shape: [Nt, N, N]
def SORTING_EIGENVECTORS_NORMALIZED_CHANGING_TSLICE(the_eigenvals, the_eigenvecs, the_t0):
    the_ref_tslice_eigenval = int((len(the_eigenvals) - the_t0) / 3)
    the_final_eigens, the_final_eigenvecs = list(the_eigenvals[:the_ref_tslice_eigenval + 1]), list(
        the_eigenvecs[:the_ref_tslice_eigenval + 1])
    the_ref_eigenvec = the_eigenvecs[the_ref_tslice_eigenval]
    for vec in range(len(the_ref_eigenvec)):
        the_ref_eigenvec[vec] = np.abs(the_ref_eigenvec[vec] / np.linalg.norm(the_ref_eigenvec[vec]))
    for ii in range(the_ref_tslice_eigenval + 1, len(the_eigenvals)):
        ckl_k = []
        for kk in range(len(the_ref_eigenvec)):
            ckl_l = []
            for ll in range(len(the_ref_eigenvec)):
                ckl_l.append(np.abs(
                    np.dot(the_ref_eigenvec[kk], the_eigenvecs[ii][ll] / np.linalg.norm(the_eigenvecs[ii][ll].T))))
            ckl_k.append(ckl_l)
        ckl_k = np.matrix(ckl_k)
        ckl = sorted([(ckl_k[i, j], i, j) for i in range(ckl_k.shape[0]) for j in range(ckl_k.shape[1])], reverse=True)
        the_used_rows, the_used_cols, the_selected = set(), set(), []
        the_top_n = None
        for the_value, the_i, the_j in ckl:
            if the_i not in the_used_rows and the_j not in the_used_cols:
                the_selected.append([the_i, the_j]);
                the_used_rows.add(the_i);
                the_used_cols.add(the_j)
                if the_top_n is not None and len(the_selected) >= the_top_n:
                    break
        the_sorted_indices = sorted(the_selected, key=lambda x: (x[0], x[1]))
        the_sorted_indices = [the_sorted_indices[i][1] for i in range(len(the_sorted_indices))][::-1] # reverse so that they also have the same order as the eigenvalue routine
        the_sorted_eigenvecs = np.array([the_eigenvecs[ii][i] for i in the_sorted_indices])
        the_final_eigens.append(np.array([the_eigenvals[ii][i] for i in the_sorted_indices]))
        the_final_eigenvecs.append(the_sorted_eigenvecs)
    return [the_final_eigens, the_final_eigenvecs]


## ------------------- BINNING -------------------------------------------

# Comments: This can be used for the reweighting factors and for the data
# a_list: data with shape [Ncfgs]
# bin_size: size of the rebinning
# it returns a smaller list len(a_list)/ bin_size
def BINNING(a_list, bin_size):
    the_a = np.asarray(a_list)
    the_n = (the_a.size // bin_size) * bin_size
    the_a_trimmed = the_a[:the_n]
    the_rebinned = the_a_trimmed.reshape(-1, bin_size).mean(axis=1)
    return the_rebinned


### Comments:
# Bin the middle segment (126:325) pairwise, handle leftover,
#     and keep first and last segments unbinned. Return binned data
#     and binning weights for the full array. Here, the binning weights account for how many configs are in a bin
#     to keep the mean unchanged later on.
# Parameters: data : 1D array; bin_size : number of configs per bin in the middle segment
# Returns: binning_weights : list of number of configs in each bin
def BINNING_CUSTOM(data, bin_size=2):

    # Segments
    set_1 = data[:126]
    set_2 = data[126:325]
    set_3 = data[325:]

    # Binned middle segment
    binned_set_2 = []
    weights_set_2 = []
    n2 = len(set_2)

    for i in range(0, n2 - n2 % bin_size, bin_size):
        chunk = set_2[i:i + bin_size]
        binned_set_2.append(np.mean(chunk, dtype=np.float128))
        weights_set_2.append(len(chunk))

    # Handle leftover
    if n2 % bin_size != 0:
        last_chunk = set_2[-1:]
        binned_set_2.append(last_chunk[0])
        weights_set_2.append(1)

    binned_set_2 = np.array(binned_set_2, dtype=np.float128)
    # Combine all segments
    binned_data = np.concatenate([set_1, binned_set_2, set_3], axis=0)

    # Binning weights: 1 for unbinned configs
    weights_set_1 = [1] * len(set_1)
    weights_set_3 = [1] * len(set_3)
    binning_weights = weights_set_1 + weights_set_2 + weights_set_3

    return np.array(binned_data, dtype=np.float128), binning_weights


### Comments:
def BINNING_CUSTOM_CORR_WEIGHTS(a):
    set_1 = a[0:126]
    set_2 = a[126:324]
    set_3 = a[324:]
    new_set = []
    for i in range(0, len(set_2) - 1, 2):  # step = 2
        new_set.append(set_2[i] + set_2[i + 1])
    new_set = np.array(new_set, dtype=np.float128)
    return np.concatenate((set_1, new_set, set_3), axis=0)



## ------------------- REWIGHTING FACTORS --------------------------------
# Comments
# rw: list of reweighting factors [Ncfgs]
# it returns a list with the normalized reweighting factors
def RW_NORMALIZATION(rw, nfs):
    len_rw = np.double(nfs)
    nrm_fktr_rw = 0. ; new_rw = []
    for ll in range(int(len_rw)):
        nrm_fktr_rw += np.double(rw[ll])
    nrm_fktr_rw = np.double( nrm_fktr_rw / len_rw)
    for ll in range(int(len_rw)):
        new_rw.append(np.double(rw[ll]) / nrm_fktr_rw)
    return np.array(new_rw, dtype=np.double)

# It reweights the correlator with the normalized reweights
# da_corr: is the correlator, it must have the shape [Ncfgs, nt]
# rw:  is the reweighting factors, and it must have the shape [Ncfgs]
# It gives back a list of the shape [nt, Ncfgs]
def REWEIGHTED_CORR(da_corr, rw):
    return np.array(da_corr.T[:,] * rw)


# rw: list of reweighting factors [Ncfgs]
# it returns a list with the normalized reweighting factors
def REWEIGHTS(the_rw_list, the_nfs):
    def process_file(the_file):
        data = np.loadtxt(the_file, unpack=True)
        if data.ndim == 1:
            return data
        ncols = data.shape[0]
        if ncols == 1:
            return data[0]
        elif ncols == 2:
            return data[1]
        elif ncols == 3:
            return data[1] * data[2]
        elif ncols == 4:
            return data[3]
        else:
            raise ValueError(f"Unsupported structure in {the_file}: {ncols} columns")
    processed = [process_file(f) for f in the_rw_list]
    if len(processed) == 1:
        the_weight = processed[0]
    else:
        the_weight = np.concatenate(processed, axis=0)
    return RW_NORMALIZATION(the_weight, the_nfs)


### ------------- RESAMPLING ------------------------------------------------------

# Comments:
# a_list: is the list of correlators that needs to be bootsrapped, shape [Ncfgs]
# c_conf: is the list generated to choose those specific configs, shape [K samples, Ncgfs]
# dis_rw: is the list of normalized rebinned reweights, shape [Ncfgs]
# It normalizes the reweights chosen, and calculates de average over the samples, mutplied by this factor
# It returns the list with the new samples, it has the shape [K bt samples]
def BOOTSTRAP(a_list, c_conf, dis_rw):
    the_a, the_rw, k_th = np.asarray(a_list), np.asarray(dis_rw), np.asarray(c_conf, dtype=int)
    the_bt_samples = the_a[k_th]
    the_rw_samples = the_rw[k_th]
    the_weighted = the_bt_samples * the_rw_samples
    k_mean_corr = np.mean(the_weighted, axis=1)
    k_norm_factor = np.mean(the_rw_samples, axis=1)
    return (k_mean_corr / k_norm_factor)


# Comments:
# it receives a list of gauge configs: [Ncfgs]
# it returns a list of configs averaged
def JACKKNIFE(a, dis_rw, bin_rw = None):
    the_a, the_rw = np.asarray(a), np.asarray(dis_rw)
    the_w = the_a * the_rw
    the_n = the_a.size
    return ((the_w.sum() - the_w) / (the_n - 1)) * ((the_rw.sum() - the_rw) / (the_n - 1))

### Comments:
# it receives a list of gauge configs: [Ncfgs]
# it returns a list of configs averaged considering the different
def JACKKNIFE_CUSTOM(a, dis_rw, bin_rw):
    a = np.array(a, dtype=np.float128)
    w = np.array(dis_rw, dtype=np.float128) * np.array(bin_rw, dtype=np.float128)

    total_sum = np.sum(a * w)
    total_weight = np.sum(w)

    jk_samples = (total_sum - a * w) / (total_weight - w)
    return np.array(jk_samples, dtype=np.double)



### ------------- STATISTICS -------------------------------------------------
#  MEANS 
# Comments:
# a:  list of [time slices (nt), nr configs (Nf)]
# it returns an array of nt slice with a mean value for each, shape [nt slices]
def MEAN(a):
    the_mean_val = np.asarray(np.mean(a, axis=1, dtype=a.dtype))
    return the_mean_val


### Comments::
# Mean for the custom bin option where not all bins are of same size
# a: list of [time slices (nt), nr configs (Nf)]
# bins: list of length nt that keeps track of how many values went into one bin
def MEAN_CUSTOM(a, bins):
    mean_val = []
    for tt in range(len(a)):
        numerator = np.sum(np.float128(a[tt]) * np.float128(bins[tt]))
        denominator = np.sum(np.float128(bins[tt]))
        mean_val.append(numerator / denominator)
    return np.array(mean_val, dtype=np.float128)




# STATISTICAL ERROR
# a: list generated from resampling, it has the shape [ nt, nconfigs]
# b: list of averaged correlator, it has the shape of [nt]
# c: receives the type of resampling and calculates de prefactor associated
# sigma: list of statistical error per each time slice, it has shape [nt]
def STD_DEV_MEAN(a,b,c):
    len_a = np.double(len(a[0]))
    if c=='jk':
        pre_factor = np.double((len_a - 1.)/len_a)
    elif c=='bt':
        pre_factor = np.double(1./len_a)
    the_diffs = a - b[:, None]               
    the_s = np.sum(the_diffs**2, axis=1)      
    the_sigma = np.sqrt(the_s * pre_factor)  
    return the_sigma


# SIMPLEST STATISTICAL ERROR
# a: is a list shape [Ncnfgs]
# b: is the mean value shape (1,)
# c: is multiplicative factor
def STD_DEV(a,b,c):
    len_a = np.double(len(a))
    if c=='jk':
        pre_factor = np.double((len_a-1.)/len_a)
    elif c=='bt':
        pre_factor = np.double(1./len_a)
    sigma = []
    for ii in range(0,len(a)):
        sigma.append(np.double((np.double(a[ii]) - np.double(b))**2))
    sigma = np.sum(sigma)
    return np.sqrt(sigma*pre_factor)



# COVARIANCE MATRIX
# a: list generated from resampling, it has the shape [nt, nconfigs]
# b: list of averaged correlator, it has the shape of [nt]
# c: prefactor associated to resampling, is a type, and it calculates the prefactor
# sigma: list of statistical error per each time slice
def COV_MATRIX(a,b,c):
    len_a = np.double(len(a[0]))
    if c=='jk':
        pre_factor = np.double((len_a - 1.)/len_a)
    elif c=='bt':
        pre_factor = np.double(1./(len_a))
    the_diffs = a - b[:, None]               
    the_s = the_diffs @ the_diffs.T      
    the_sigma = the_s * pre_factor
    return np.matrix(the_sigma)



### ------------- FITTING FUNCTIONS AND STUFF --------------------------------------------

# Comments:
# This is the function for a single exponential fit. Ae^{-E0*(nt-t0))}
# x: is the t_slices
# the_pars: is a list (amplitud, energy)
# # *a: this is a variable size arguments, in thie case corresponds to t0.
def SINGLE_EXPONENTIAL(x,the_pars,*a): 
    a0, e0 = the_pars
    return a0 * np.exp(-e0 * (x - a))


### Comments:
# This is the function for an alternative double exponential fit. Ae^{-E0*(nt-t0)}(1+ B e^{-nt (E1-E0)})
# x: is the t_slices
# e0: is a list (amplitud, energy E0, Amplitude shift of energy, DeltaE**2)
# *a: this is a variable size arguments, in thie case corresponds to t0.
def DOUBLE_EXPONENTIAL_ALTERNATIVE(x, e0, *a):
    return e0[0] * np.exp(-(x - a) * e0[1]) * (1. + e0[2] * np.exp(-x * (e0[3] - e0[1])))


# Comments:
# This is the function for a geometric fit form: Ae^{-E0(nt-t0)}/ (1 - B e^{- nt*M})
# x: is the t_slices
# the_pars: is a list (amplitud, energy E0, Amplitude shift of energy)
# *a: this is a variable size arguments, in thie case corresponds to t0.
def GEOMETRIC_FORM(x,the_pars,*a):
    a0, e0, b, dm = the_pars
    return a0 * np.exp(-(x-a) * e0) / (1. - b * np.exp(-x * dm))


### Comments:
# This is the function for a double exponential fit. Ae^{-E0*(nt-t0)}(1+ Be^{-nt*D^{2}})
# x: is the t_slices
# the_pars: is a list (amplitud, energy E0, Amplitude shift of energy, DeltaE**2)
# *a: this is a variable size arguments, in thie case corresponds to t0.
def DOUBLE_EXPONENTIAL(x, the_pars, *a):
    a0, e0, b, dm = the_pars
    return a0 * np.exp(-(x - a) * e0) * (1. + b * np.exp(-x * (dm ** 2)))


### Comments:
# This function tries to find a good guess for the fit to have a prior, so it would in principle take less time. It uses a simple polynomial fit of order 1. 
def BEST_GUESS(c,t_i,tipo_fit):
    da_fits = np.polyfit(t_i, np.log(c), deg=1)
    amp0 = np.double(np.exp(da_fits[1]))
    eng0 = np.double(-da_fits[0])
    guess = []
    guess.append(amp0)
    guess.append(eng0)
    if tipo_fit=='2' or tipo_fit=='g':
        guess.append(np.double(0.1))
        guess.append(np.double(0.5))
    return guess


### Comments: 
# This function gets the Difference between a Chi^{2} of one time slice compared to the next time slice value of Chi^{2}. This in order to check for stability.
def DELTA_CHI(a):
    delta_chi = [a[ii ] - a[ii + 1] for ii in range(len(a)-1)]
    return np.asarray(delta_chi)


### Comments:
# This function gets the total Chi^{2}, not only per degree of freedom. 
# a: this is the list with the chi^{2} per degree of freedom
# b: the list of time slices
# c: the total number of time slices
def TOTAL_CHI(a,b,c,nrp):
    total_a = [a[ii] * (np.double(c[ii] - np.double(b[ii])) - np.double(nrp)) for ii in range(len(a))]
    return np.asarray(total_a)



### Comments:
# This class gets the function to do the Chi^{2} fit. 
class My_Fits:
    def __init__(self, model, x, y, cov, dgof, a):
        self.model = model  # model predicts y for given x
        self.x = np.asarray(x, dtype = float)
        self.y = np.asarray(y, dtype = float)
        self.cov = np.matrix(cov, dtype = float)
        
        self.arg = a
        self.nrPars = np.double(len(dgof))
        self.dof = len(self.x) - self.nrPars

    def __call__(self, *par):  # we accept a variable number of model parameters
        ym = self.model(self.x, *par, self.arg)
        r = self.y - ym
        
        the_chi2 = r @ self.cov @ r
        the_val_chi2 = float(the_chi2 / self.dof)
        
        return the_val_chi2
    
    
    
### Comments:
# This class gets the function to do the Chi^{2} fit.
class My_Fits_Update:
    def __init__(self, model, x, y, cov, dgof, a):
        self.model = model
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.cov = np.asarray(cov, dtype=float)

        self.npar = np.double(len(dgof))
        self.arg = a
        self.dof = len(self.x) - self.npar

    def __call__(self, *par):
        ym = self.model(self.x, *par, self.arg)
        r = self.y - ym

        chi2 = r @ self.cov @ r
        val = float(chi2 / self.dof)

        return val
    
    
    
def DOING_THE_FITTING(the_corr, the_nt, the_type_rs, the_irreps, the_irrep, tmin_data, the_type_correlated_fit, the_type_fit, the_only_one_tmin, the_t0, the_list_tmaxs, da_minimization, the_fit_params):
    
    ### Loop over all the t0 chosen
    the_corr_fit = np.asarray(the_corr['Eigenvalues/Mean']).real 
    the_corr_fit_rs = np.asarray(the_corr['Eigenvalues/Resampled']).real 
    
    ### Retrieving the covariance matrix for the fits
    the_cov_matrix = np.asarray(the_corr['Eigenvalues/Covariance_matrix']).real 
    
    ### Effective Masses are used as a prior
    the_eff_energy_hint = np.asarray(the_corr['Effective_masses/Mean'])
    
    ### Creating a folder for correlated or uncorrelated fits
    if the_type_correlated_fit=='Correlated':
        if 'Correlated' in tmin_data.keys(): del tmin_data['Correlated']
        the_fit_data = tmin_data.create_group('Correlated')
    elif the_type_correlated_fit=='Uncorrelated':
        if 'Uncorrelated' in tmin_data.keys(): del tmin_data['Uncorrelated']
        the_fit_data = tmin_data.create_group('Uncorrelated')
    
    ### Creating a folder for the fits of the central values
    the_mean_data = the_fit_data.create_group('Mean')
    
    ### Creating a folder for the resampled data fits
    the_rs_data = the_fit_data.create_group('Resampled')
    
    ### Loop over the nr. of eigenvalues
    for ls in range(len(the_corr_fit)):
        
        ### This is the modified time range to go from item 0 until the end without caring about labeling
        nt_mod = np.arange(0,len(the_nt))
        if the_only_one_tmin: 
            ### Upper limit for the fit
            the_ul = [x-the_nt[0] for x in the_list_tmaxs[the_irreps.index(the_irrep)]]
            
            ### Lower limit for the fit
            the_ll = [nt_mod[5]]                    
        else: 
            ### Upper limir for the fit
            the_ul = [x-the_nt[0] for x in the_list_tmaxs[the_irreps.index(the_irrep)]]
            
            ### Lower limit for the fit depends on the upper limit
            the_ll = np.arange(nt_mod[0]+1, int(nt_mod[-1]*.75))
            
        ### Choosing the covariance matrix depending on the type of correlated fit
        if the_type_correlated_fit=='Correlated':
            the_cov_matrix_fit = the_cov_matrix[ls]
        elif the_type_correlated_fit=='Uncorrelated':
            the_cov_matrix_fit = np.diag(np.diag(the_cov_matrix[ls]))        
        
        another_list = []
        
        the_results = {'the_energies': [], 'the_sigmas': [], 'the_chi_vals': [], 'the_sigmas_chi': [], 'the_resampled': []}
        
        ### Loop over all the tmins
        for yy in the_ll:
            print(f'Tmin = {yy+the_nt[0]}  || TMax = {the_ul[ls]+the_nt[0]}')
            another_useful_list = []
            
            the_nt_slice = the_nt[yy:the_ul[ls]]
            the_corr_fit_slice = the_corr_fit[ls,yy:the_ul[ls]]
            
            ### This is finding a good guess to make the fit converge easier.
            da_hint = BEST_GUESS(the_corr_fit_slice, the_nt_slice, the_type_fit) 
            if False in np.isnan(da_hint):
                the_dof = da_hint
            else: 
                the_dof = np.zeros((1,len(da_hint)));
                the_dof = the_dof[0]
                the_dof[0] = np.float64(0.1)
                the_dof[1] = np.float64(the_eff_energy_hint[ls,yy])
            
            ### Reduces the covariance matrix to the size of the time range chosen and takes the inverse            
            the_inverse_cov_m = np.linalg.inv(SHRINK_MATRIX(the_cov_matrix_fit, yy, the_ul[ls]))
            
            ### This chooses the fit function to use 
            the_fit_choice = My_Fits(da_minimization, the_nt_slice, the_corr_fit_slice, the_inverse_cov_m, the_dof, np.float64(the_t0))
            
            ### Fitting started
            the_fit = Minuit(the_fit_choice, the_dof, name = the_fit_params)
            
            the_fit.errordef, the_fit.tol = 1e-8, 1e-10
            the_fit.scan()
            the_fit.migrad(iterate=10,ncall=5000)
            
            ### Energy values
            e0 = np.float128(the_fit.values['e0'])
            
            ### Amplitude values
            a0 = np.float128(the_fit.values['a0'])
            
            ### The fitted energy results from the central values are used as an initial guess for the resamples
            the_dof_rs = the_dof
            the_dof_rs[0], the_dof_rs[1] = a0, e0
            
            another_useful_list.append(e0)
            
            the_results['the_energies'].append(e0)
            the_results['the_chi_vals'].append(the_fit.fval)

            chi_vals_rs_list = []
            the_corr_fit_rs_eigen = the_corr_fit_rs[ls]
            ### Loop over the resamples
            for zz in range(the_corr_fit_rs.shape[1]):
                my_fit_choice_rs = My_Fits(da_minimization, the_nt_slice, the_corr_fit_rs_eigen[zz,yy:the_ul[ls]], the_inverse_cov_m, the_dof_rs, np.float64(the_t0))
                the_fit_rs = Minuit(my_fit_choice_rs, the_dof_rs, name = the_fit_params)
                
                the_fit_rs.errordef, the_fit_rs.tol = 1e-8, 1e-7
                the_fit_rs.scan()
                the_fit_rs.migrad(iterate=10, ncall=5000)
                
                e0_rs = np.float64(the_fit_rs.values['e0']); chi_vals_rs_list.append(the_fit_rs.fval); another_useful_list.append(e0_rs)
            
            ### This is the sigma for the fittings
            sigma_fit_rs = STD_DEV(another_useful_list[1:], np.mean(another_useful_list[1:]), the_type_rs)
            sigma_chi_rs = STD_DEV(chi_vals_rs_list, np.mean(chi_vals_rs_list), the_type_rs)

            another_list.append(np.array(another_useful_list))
            
            the_results['the_sigmas'].append(sigma_fit_rs)
            the_results['the_sigmas_chi'].append(sigma_chi_rs)
        print(f'E = {ls} READY')    
        
        the_rs_data.create_dataset(f'lambda_{ls}', data = np.asarray(another_list))
        the_mean_data.create_dataset(f'lambda_{ls}', data = np.asarray([the_ll + the_nt[0], [the_ul[ls]+the_nt[0]]*len(the_ll), the_results['the_energies'], the_results['the_sigmas'], the_results['the_chi_vals'], the_results['the_sigmas_chi']]))    
        
        


### --------------------------------- DISPERSION RELATION --------------------------

### Comments:
def CONTINUUM_DISP_REL(p, E0, norm):
    E = np.sqrt(E0 ** 2 + int(p) * norm) # p is already squared so no need to square it again
    E_norm = E / E0  # aE
    return E_norm

### Comments:
# calculates the fraction of the measured energy E (lat units) and the values of the continuum dispersion relation for the measured E0
def RELATIVE_DISP_REL(p, E0, E, norm):
    E_disp = np.sqrt(E0 ** 2 + int(p) * norm) # p is already squared so no need to square it again
    return E / E_disp
    
    

###---------------------------------CORRELATED DIFFERNCES -------------------------

### Comments:
def CENTER_OF_MASS(p, E_lab, norm):
    E = np.sqrt(E_lab ** 2 - int(p) * norm)  # p is already squared so no need to square it again
    return np.double(E)

### Comments:
def CONTINUUM_DISP_REL_NO_NORM(p, E0, norm):
    E = np.sqrt(E0 ** 2 + int(p) * norm)
    return E

### Comments:
# WHAT IS THIS?
def EXTRACT_HADS_MOM(non_int_list):
    return None

### Comments:
def NON_INTERACTING_LEVELS(non_int_list, t_mins_range, singles_files, t_mins_shift, norm):

    single_hadrons_dict_keys = list(singles_files.keys()) # e.g. 'K', 'P'
    single_hadrons_all_irreps = { key: list(singles_files[key].keys()) for key in single_hadrons_dict_keys }

    mesons = single_hadrons_dict_keys
    singles_dict = {}
    nr_rs = int()
    for meson in mesons:
        singles_dict[meson] = {}
        for i, irrep in enumerate(single_hadrons_all_irreps[meson]):
            dataset = singles_files[meson][irrep]
            try:
                singles_dict[meson][irrep] = {
                    'Mean': dataset.get('1exp/Tmin/Correlated/Mean')[()][2, t_mins_range[meson][i] - t_mins_shift],
                    'Resampled': dataset.get('1exp/Tmin/Correlated/Resampled')[()][t_mins_range[meson][i] - t_mins_shift]}
                nr_rs = dataset.get('1exp/Tmin/Correlated/Resampled')[()][t_mins_range[meson][i] - t_mins_shift].size
            except KeyError:
                raise KeyError(f'Irrep {irrep} does not exist in {meson} file.')
            except TypeError:
                raise TypeError(f'Fit data in irrep {irrep} is missing in {meson} file.')
            except IndexError:
                raise IndexError(f'Chosen fit range does not exist in {meson} file.')
        singles_files[meson].close()


    # calculate the dispersion energies
    non_int_energies_disp = []
    non_int_energies_disp_rs = []
    non_int_energies_mean = []
    non_int_energies_rs = []
    for i in range(len(non_int_list)):
        non_int_energies_irrep_disp = []
        non_int_energies_irrep_disp_rs = []
        non_int_energies_mean_i = []
        non_int_energies_rs_i = []
        for string in non_int_list[i]:
            split =  re.findall(r'[A-Za-z]|\d+', string)
            letters = [x for x in split if x.isalpha()]
            numbers = [int(x) for x in split if x.isdigit()]
            single_level_disp = 0
            single_level_mean = 0
            single_level_disp_rs = np.zeros([nr_rs])
            single_level_rs = np.zeros([nr_rs])
            for ii, letter in enumerate(letters):
                irrep = list(singles_dict[letter].keys())
                val = singles_dict[letter][irrep[0]]['Mean']
                single_level_disp += CONTINUUM_DISP_REL_NO_NORM(numbers[ii], singles_dict[letter][irrep[0]]['Mean'], norm)
                single_level_disp_rs += np.array([CONTINUUM_DISP_REL_NO_NORM(numbers[ii], x, norm) for x in singles_dict[letter][irrep[0]]['Resampled']])
                match = next(
                    (irrep for irrep in single_hadrons_all_irreps[letter]
                     if re.search(rf'PSQ{numbers[ii]}_', irrep)),
                    None
                )

                if match is None:
                    print(f"No match for {letter} with PSQ{numbers[ii]}")
                    continue
                single_level_mean += singles_dict[letter][match]['Mean']
                resample = singles_dict[letter][match]['Resampled']
                single_level_rs += resample
            non_int_energies_mean_i.append(single_level_mean)
            non_int_energies_rs_i.append(list(single_level_rs))
            non_int_energies_irrep_disp.append(single_level_disp)
            non_int_energies_irrep_disp_rs.append(list(single_level_disp_rs))
        non_int_energies_disp.append(non_int_energies_irrep_disp)
        non_int_energies_disp_rs.append(non_int_energies_irrep_disp_rs)
        non_int_energies_mean.append(non_int_energies_mean_i)
        non_int_energies_rs.append(non_int_energies_rs_i)

    return non_int_energies_disp, non_int_energies_disp_rs, non_int_energies_mean, non_int_energies_rs

def UNIT_CONVERSION(val, a, MeV):
    return val / a * MeV




### ---------------------------------- RATIO OF CORRELATORS -------------------------------
### Comments:
# This class gets the info of the non-interacting levels. It only works when there is a threshold of 2 states nearby. It needs modification for the 3particle threshold.
class NonInteractingLevels:
    def __init__(self,nombre,all_hads):
        self.FirstState = str(nombre[:4])
        self.SecondState =  str(nombre[4:])
        FirstRaw = 'PSQ'+ str(self.FirstState[2])+'_'+ str(self.FirstState[0])
        SecondRaw = 'PSQ'+ str(self.SecondState[2])+'_'+ str(self.SecondState[0])
        for item in all_hads:
            if FirstRaw[:4] in item and FirstRaw[-1]==item[-1].capitalize(): first_name = item;break
            else: continue
        for item in all_hads:
            if SecondRaw[:4] in item and SecondRaw[-1]==item[-1].capitalize(): second_name = item;break
            else: continue
        self.First = first_name
        self.Second = second_name




### ---------------------------------- FILE READ OUT -------------------------------
### Comments:
# Returns the keys and data sets stored in the file and their type
def READ_STRUCTURE(file):
    with h5py.File(file, 'r') as f:
        print("Analyzing the structure of saved file...")

        def describe_type(value):
            """Return a human-readable type description for HDF5 attributes."""
            if isinstance(value, (str, np.str_)):
                return "string (unicode)"
            elif isinstance(value, (bytes, np.bytes_)):
                return "string (bytes/ASCII)"
            elif isinstance(value, np.ndarray) and value.dtype.kind in {"S", "U"}:
                return f"array of strings (dtype={value.dtype})"
            else:
                return type(value).__name__
        # Function to print the structure of the HDF5 file
        def print_structure(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"Dataset: {name} - Shape: {obj.shape} - Datatype: {obj.dtype}")
            elif isinstance(obj, h5py.Group):
                print(f"Group: {name}")
            # Print attributes if any
            if obj.attrs:
                print(f"  Attributes for {name}:")
                for key, value in obj.attrs.items():
                    if isinstance(value, np.ndarray):
                        print(f"    {key}: shape={value.shape}, dtype={value.dtype}, "
                              f"type={describe_type(value)}, value={value}")
                    else:
                        print(f"    {key}: type={describe_type(value)}, value={value}")

        # Visit each object in the file and print its name and type
        f.visititems(print_structure)
        print("Reading Complete.")
        


if __name__=="__main__":
   print('Nothing to run here, unless you want to change something.')
