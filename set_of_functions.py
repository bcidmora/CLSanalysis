import numpy as np
import os
from pathlib import Path


## ------------------- SOME USEFUL FUNCTIONS -----------------------------------
# Comments
# k_size: is the size of the bootstrap sample
# cfgs_nr: is the amount of gauge configs after rebinning
# it returns a list of list of numbers to plug into the boostrap sampling
def RANDOM_GENERATOR(k_size, cfgs_nr):
    np.random.seed(414557)
    random_list = []
    for ii in range(k_size):
        this_random = []
        for jj in range(cfgs_nr):
            i_cfg = int(np.random.random() * cfgs_nr)
            this_random.append(i_cfg)
        random_list.append(np.array(this_random))
    return np.array(random_list)

# It receives a list [Ncfgs]
# It returns the normalization factor
def NORM_FACTOR(a_list):
    # return np.double(np.mean(a_list)) #OLD
    return np.double(np.mean(a_list, dtype=np.float128))


# It checks for hermiticity
# a_matrix: is the matrix [N,N]
def MAKES_HERMITIAN(a_matrix):
    return np.float128(0.5)*(np.matrix(a_matrix)+np.conj(np.matrix(a_matrix).T))     
    


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
# This function prints the information of the ensemble and the type of correlators one have.
def INFO_PRINTING(the_corr_type, the_ensemble):
    if the_corr_type=='s':
        print('.............................................................................')
        print('                         SINGLE HADRON CORRELATORS')
        print('                               ENSEMBLE '+ the_ensemble)
        print('.............................................................................')
    elif the_corr_type=='m':
        print('.............................................................................')
        print('                         MULTIHADRON CORRELATORS')
        print('                               ENSEMBLE '+ the_ensemble)
        print('.............................................................................')
    elif the_corr_type=='mr':
        print('.............................................................................')
        print('                       MULTIHADRON RATIO CORRELATORS')
        print('                               ENSEMBLE '+ the_ensemble)
        print('.............................................................................')

## Comments:
# This function checks if a directory exists, if not then it creats it.
def DIRECTORY_EXISTS(a_dir):
    if not os.path.isdir(a_dir):
        Path(a_dir).mkdir(parents=True, exist_ok=True)
        new_dir = a_dir
    else:
        new_dir=a_dir
    return new_dir
        

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
    ncfgs = len(a)
    s = len(a[0])
    reshaped_corr=[]
    for n1 in range(s):
        reshaped_corr_n1=[]
        for n2 in range(s):
            reshaped_corr_n2=[]
            for nf in range(ncfgs):
                reshaped_corr_n2.append(np.array(a[nf][n1][n2]))
            reshaped_corr_n1.append(np.array(reshaped_corr_n2))
        reshaped_corr.append(np.array(reshaped_corr_n1))
    return np.array(reshaped_corr)


# This function reshapes the correlators as: [nt, N, N] --> [N,N,nt]
# a: the list to reshape
# s: size of the matrix
def RESHAPING_CORRELATORS(a):
    new_corr=[]
    s=len(a[0])
    for n1 in range(s):
        new_corr_n =[]
        for n2 in range(s):
            new_corr_t = []
            for dis_t in range(len(a)):
                new_corr_t.append(a[dis_t][n1][n2])
            new_corr_n.append(np.array(new_corr_t))
        new_corr.append(np.array(new_corr_n))
    return np.array(new_corr)

# This function reshapes the correlators as: [nt, Ncfgs, N, N] --> [N, N, Ncgfs, nt]
# a: the list to reshape
# s: size of the matrix
def RESHAPING_CORRELATORS_RS(a):
    new_corr=[]
    s = a.shape[-1]
    for n1 in range(s):
        new_corr_n =[]
        for n2 in range(s):
            new_corr_nfs = []
            for nfs in range(len(a[0])):
                new_corr_t=[]
                for dis_t in range(len(a)):
                    new_corr_t.append(a[dis_t][nfs][n1][n2])
                new_corr_nfs.append(np.array(new_corr_t))
            new_corr_n.append(np.array(new_corr_nfs))
        new_corr.append(np.array(new_corr_n))
    return np.array(new_corr)

# This function reshapes the correlators as: [nt, Ncfgs, N, N] --> [N, N, nt, Ncgfs]
# a: the list to reshape
# s: size of the matrix
def RESHAPING_CORRELATORS_RS_NT(a):
    s=a.shape[-1]
    new_corr=[]
    for n1 in range(s):
        new_corr_n =[]
        for n2 in range(s):
            new_corr_nfs = []
            for dis_t in range(len(a)):
                new_corr_t=[]
                for nfs in range(len(a[0])):
                    new_corr_t.append(a[dis_t][nfs][n1][n2])
                new_corr_nfs.append(np.array(new_corr_t))
            new_corr_n.append(np.array(new_corr_nfs))
        new_corr.append(np.array(new_corr_n))
    return np.array(new_corr)



# This function reshapes the data in the way: [N, N, nt] -> [nt,N,N] This is used to later get the eigenvalues easier form the matrices.
# a: list with the data
def RESHAPING_EIGENVALS_MEAN(a,s):
    t_slices = a.shape[-1]
    eig_corrs = []
    for da_times in range(t_slices):
        time_corrs = []
        for n1 in range(s):
            corr_n1=[]
            for n2 in range(s):
                corr_n1.append(a[n1][n2][da_times])
            time_corrs.append(np.array(corr_n1))
        eig_corrs.append(np.array(time_corrs))
    return np.array(eig_corrs)


# This function reshapes the data in the way: [N, N, nt, Ncfgs] -> [nt,Ncfgs N,N]
# a: list with the data
def RESHAPING_EIGENVALS_RS(a,s):
    t_slices = a.shape[-2]
    ncfgs = a.shape[-1]
    eig_corrs = []
    for da_times in range(t_slices):
        time_corrs = []
        for nf in range(ncfgs):
            corr_n1=[]
            for n1 in range(s):
                corr_n2=[]
                for n2 in range(s):
                    corr_n2.append(a[n1][n2][da_times][nf])
                corr_n1.append(np.array(corr_n2))
            time_corrs.append(np.array(corr_n1))
        eig_corrs.append(np.array(time_corrs))
    return np.array(eig_corrs)

# This function reshapes the eigenvalues from [nt, Ncfgs, Neigens] -->> [Neigens, Ncfgs, nt]. This is used to obatin the fits easier later. 
# a: is the list of eigenvals, n configs, and time slices.
# s: size of the matrix or amount of eigenvalues
def RESHAPING_EIGENVALS_FOR_FITS(a,s):
    t_slices = a.shape[0]
    ncfgs = a.shape[1]
    eig_corrs = []
    for n1 in range(s):
        corr_n=[]
        for nf in range(ncfgs):
            corr_t=[]
            for t_s in range(t_slices):
                corr_t.append(a[t_s][nf][n1])
            corr_n.append(np.array(corr_t))
        eig_corrs.append(np.array(corr_n))
    return np.array(eig_corrs)


# This function reshapes the eigenvalues from [nt, Ncfgs, Neigens] -->> [Ncfgs, nt, Neigens]. This is used to obatin the fits easier later. 
# a: is the list of eigenvals, n configs, and time slices.
def RESHAPING_EIGEN_FOR_SORTING(a):
    t_slices = a.shape[0]
    ncfgs = a.shape[1]
    eig_corrs = []
    for nf in range(ncfgs):
        corr_n=[]
        for t_s in range(t_slices):
            corr_n.append(np.array(a[t_s][nf]))
        eig_corrs.append(np.array(corr_n))
    return np.array(eig_corrs)


# This function reshapes the eigenvalues from [Ncfgs, nt, Neigens] -->> [nt, Ncfgs, Neigens] . This is used to obatin the fits easier later. 
# a: is the list of eigenvals, n configs, and time slices.
def RESHAPING_EIGEN_FOR_SORTING_REVERSE(a):
    t_slices = a.shape[1]
    ncfgs = a.shape[0]
    eig_corrs = []
    for t_s in range(t_slices):
        corr_n=[]
        for nf in range(ncfgs):
            corr_n.append(np.array(a[nf][t_s]))
        eig_corrs.append(np.array(corr_n))
    return np.array(eig_corrs) 

# This function reshapes the eigenvalues from [Neigens, Ncfgs, nt] -->> [nt, Ncfgs, Neigens] 
# a: is the list of eigenvals, n configs, and time slices.
def RESHAPING_EIGENVALS_NN(a):
    t_slices = a.shape[2]
    ncfgs = a.shape[1]
    n_eigvals = a.shape[0]
    eig_corrs = []
    for t_s in range(t_slices):
        corr_n=[]
        for nf in range(ncfgs):
            corr_t=[]
            for n1 in range(n_eigvals):
                corr_t.append(a[n1][nf][t_s])
            corr_n.append(np.array(corr_t))
        eig_corrs.append(np.array(corr_n))
    return np.array(eig_corrs)

# This function reshapes the eigenvalues from [nt, Ncfgs] -->> [Ncfgs, nt]. This is used to obatin the fits easier later. 
# a: is the list of n configs, and time slices.
# s: number of data entries
def NT_TO_NCFGS(a):
    eig_corrs = []
    for n1 in range(len(a[0])):
        corr_n=[]
        for t_slices in range(len(a)):
            corr_n.append(a[t_slices][n1])
        eig_corrs.append(np.array(corr_n))
    return np.array(eig_corrs)

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
    modified_cov_mat = []
    for x in range(low,up):
        cov_ij=[]
        for z in range(low,up):
            cov_ij.append(c[x][z])
        modified_cov_mat.append(np.array(cov_ij))
    out_small_cov = np.array(modified_cov_mat, dtype=np.double)
    return out_small_cov

#Comments:
# This functions removes cols and rows from  a correlation matrix, leaving it squared as it should be. 
# c: this is the mean value of correlator matrix. It has a shape of [N, N, nt]
# r: This is the resampled correlator matrix. It has the shape [N, N, nt, Ncfgs]
# ss: this is the index of the row/column to be removed
def REMOVE_ROWS_COLS(c,r,ss):
    modified_mean_corr = []
    modified_rs_corr  = []
    for ij in range(len(c)):
        the_mean_corr_ij = []
        the_rs_corr_ij = []
        if ij!=ss:
            for ji in range(len(c)):
                if ji!=ss:
                    the_mean_corr_ij.append(np.array(c[ij][ji]))
                    the_rs_corr_ij.append(np.array(r[ij][ji]))
            modified_mean_corr.append(np.array(the_mean_corr_ij))
            modified_rs_corr.append(np.array(the_rs_corr_ij))
    return [np.array(modified_mean_corr), np.array(modified_rs_corr)]


#Comments:
# This functions adds cols and rows from  a correlation matrix, leaving it squared as it should be. 
# c: this is the mean value of correlator matrix. It has a shape of [N, N, nt]
# r: This is the resampled correlator matrix. It has the shape [N, N, nt, Ncfgs]
# ss: this is the index of the row/column to be added
def ADD_ROWS_COLS(c,r,ss):
    modified_mean_corr = []
    modified_rs_corr  = []
    for ij in range(ss):
        the_mean_corr_ij = []
        the_rs_corr_ij = []
        for ji in range(ss):
            the_mean_corr_ij.append(np.array(c[ij][ji]))
            the_rs_corr_ij.append(np.array(r[ij][ji]))
        modified_mean_corr.append(np.array(the_mean_corr_ij))
        modified_rs_corr.append(np.array(the_rs_corr_ij))
    return [np.array(modified_mean_corr), np.array(modified_rs_corr)]
            

## ------------------- EFFECTIVE MASSES -----------------------------------

# This function receives a list "a" of time slices, and it calculates the effective mass, returning a list of effectives masses for each time slice (half integer numbers).
# a: shape [nt]
# This is the old version of the Effective Masses. It can still be used, but the other one is more general
def EFF_MASS_CLASSIC(a):
    meff=[]
    for i in range(len(a)-1):
        meff.append(np.log(np.double(a[i])/np.double(a[i+1])))
    return np.array(meff)


# This function receives a list "a" of time slices, and it calculates the effective mass, returning a list of effectives masses for each time slice (half integer numbers).
# a: shape [nt]
# d: distance of two points, by default this is 1
def EFF_MASS(a,d):
    meff=[]
    for i in range(len(a)-d):
        meff.append(np.log(np.abs(np.double(a[i])/np.double(a[i+d]))))
    return np.array(meff)

# This function receives a list "a" of time slices, and it calculates the effective mass, returning a list of effectives masses for each time slice (half integer numbers).
# a: shape [nt]
# d: distance of two points, by default this is 1
def EFF_MASS_COSH(a,d):
    meff=[]
    for i in range(len(a)-d):
        meff.append(np.acosh(np.abs((np.double(a[i+d]) + np.double(a[i-d]))/np.double(2.*a[i]))))
    return np.array(meff)


## ------------------- SORTING STATES --------------------------------------


# Comments: This function returns the new order based on the eigenvalues. They are ordered such that the largest eigenvalue and its corresponding eigenvector is at the right spot for all time slices based on a reference time slice. 
# the_t0: this is the t0 used for the GEVP, the eigenvalues must be sorted after this reference time slice.
# the_eigenvals: This is an array of the values of the eigenvals: shape [Nt, Neigens]
# the_eigenvecs: This is an array of the eigenvectors in the following shape: [Nt, N, N]
def SORTING_EIGENVALUES(the_t0, the_eigenvals, the_eigenvecs):
    the_final_eigens = list(the_eigenvals[:the_t0+1])
    the_final_eigenvecs = list(the_eigenvecs[:the_t0+1])
    for ii in range(the_t0+1, len(the_eigenvals)):
        the_sorted_indices = sorted(range(len(the_eigenvals[ii])), key=lambda i: the_eigenvals[ii][i], reverse=True)  
        the_final_eigens.append(np.array([the_eigenvals[ii][i] for i in the_sorted_indices]))
        the_final_eigenvecs.append(np.array([the_eigenvecs[ii][i] for i in the_sorted_indices]))
    return [the_final_eigens, the_final_eigenvecs]


# Comments: This function checks for orthogonality of eigenvectors and returns the new order with the eigenvectors ordered such that they are associated to the corresponding state time slice by time slice. No normalization of the vectors
# the_eigenvals: This is an array of the values of the eigenvals: shape [Nt, Neigens]
# the_eigenvecs: This is an array of the eigenvectors in the following shape: [Nt, N, N]
def SORTING_EIGENVECTORS(the_t0, the_eigenvals, the_eigenvecs):
    the_ref_tslice_eigenval = int((len(the_eigenvals)-the_t0)/3)
    the_final_eigens, the_final_eigenvecs = list(the_eigenvals[:the_ref_tslice_eigenval+1]), list(the_eigenvecs[:the_ref_tslice_eigenval+1])
    the_ref_eigenvec = the_eigenvecs[the_ref_tslice_eigenval]
    for ii in range(the_ref_tslice_eigenval+1, len(the_eigenvals)):
        ckl_k = []
        for kk in range(len(the_ref_eigenvec)):
            ckl_l = []
            for ll in range(len(the_ref_eigenvec)):
                ckl_l.append(np.abs(np.dot(the_ref_eigenvec[kk], the_eigenvecs[ii][ll].T)))
            ckl_k.append(ckl_l)
        ckl_k = np.matrix(ckl_k)
        ckl = sorted([(ckl_k[i, j], i, j) for i in range(ckl_k.shape[0]) for j in range(ckl_k.shape[1])], reverse=True)
        the_used_rows, the_used_cols, the_selected = set(), set(), []
        the_top_n=None
        for the_value, the_i, the_j in ckl:
            if the_i not in the_used_rows and the_j not in the_used_cols:
                the_selected.append([the_i, the_j]); the_used_rows.add(the_i) ; the_used_cols.add(the_j)
                if the_top_n is not None and len(the_selected) >= the_top_n:
                    break
        the_sorted_indices = sorted(the_selected, key=lambda x: (x[0], x[1]))
        the_sorted_indices = [the_sorted_indices[i][1] for i in range(len(the_sorted_indices))]
        the_sorted_eigenvecs = np.array([the_eigenvecs[ii][i] for i in the_sorted_indices])
        the_final_eigens.append(np.array([the_eigenvals[ii][i] for i in the_sorted_indices]))
        the_final_eigenvecs.append(the_sorted_eigenvecs)
    return [the_final_eigens, the_final_eigenvecs]



# Comments: This function checks for orthogonality of eigenvectors and returns the new order with the eigenvectors ordered such that they are associated to the corresponding state time slice by time slice. Each vector is normalized first.
# the_eigenvals: This is an array of the values of the eigenvals: shape [Nt, Neigens]
# the_eigenvecs: This is an array of the eigenvectors in the following shape: [Nt, N, N]
def SORTING_EIGENVECTORS_NORMALIZED(the_t0, the_eigenvals, the_eigenvecs):
    the_ref_tslice_eigenval = int((len(the_eigenvals)-the_t0)/3)
    the_final_eigens, the_final_eigenvecs = list(the_eigenvals[:the_ref_tslice_eigenval+1]), list(the_eigenvecs[:the_ref_tslice_eigenval+1])
    the_ref_eigenvec = the_eigenvecs[the_ref_tslice_eigenval]
    for vec in range(len(the_ref_eigenvec)):
        the_ref_eigenvec[vec] = np.abs(the_ref_eigenvec[vec]/np.linalg.norm(the_ref_eigenvec[vec]))
    for ii in range(the_ref_tslice_eigenval+1, len(the_eigenvals)):
        ckl_k = []
        for kk in range(len(the_ref_eigenvec)):
            ckl_l = []
            for ll in range(len(the_ref_eigenvec)):
                ckl_l.append(np.abs(np.dot(the_ref_eigenvec[kk], the_eigenvecs[ii][ll]/np.linalg.norm(the_eigenvecs[ii][ll].T))))
            ckl_k.append(ckl_l)
        ckl_k = np.matrix(ckl_k)
        ckl = sorted([(ckl_k[i, j], i, j) for i in range(ckl_k.shape[0]) for j in range(ckl_k.shape[1])], reverse=True)
        the_used_rows, the_used_cols, the_selected = set(), set(), []
        the_top_n=None
        for the_value, the_i, the_j in ckl:
            if the_i not in the_used_rows and the_j not in the_used_cols:
                the_selected.append([the_i, the_j]); the_used_rows.add(the_i) ; the_used_cols.add(the_j)
                if the_top_n is not None and len(the_selected) >= the_top_n:
                    break
        the_sorted_indices = sorted(the_selected, key=lambda x: (x[0], x[1]))
        the_sorted_indices = [the_sorted_indices[i][1] for i in range(len(the_sorted_indices))]
        the_sorted_eigenvecs = np.array([the_eigenvecs[ii][i] for i in the_sorted_indices])
        the_final_eigens.append(np.array([the_eigenvals[ii][i] for i in the_sorted_indices]))
        the_final_eigenvecs.append(the_sorted_eigenvecs)
    return [the_final_eigens, the_final_eigenvecs]


# Comments: This function checks for orthogonality of eigenvectors and returns the new order with the eigenvectors ordered such that they are associated to the corresponding state time slice by time slice. No normalization of the vectors, and changing the reference time slice to compare with.
# the_eigenvals: This is an array of the values of the eigenvals: shape [Nt, Neigens]
# the_eigenvecs: This is an array of the eigenvectors in the following shape: [Nt, N, N]
def SORTING_EIGENVECTORS_CHANGING_TSLICE(the_t0, the_eigenvals, the_eigenvecs):
    the_ref_tslice_eigenval = int((len(the_eigenvals)-the_t0)/3)
    the_final_eigens, the_final_eigenvecs = list(the_eigenvals[:the_ref_tslice_eigenval+1]), list(the_eigenvecs[:the_ref_tslice_eigenval+1])
    the_ref_eigenvec = the_eigenvecs[the_ref_tslice_eigenval]
    for ii in range(the_ref_tslice_eigenval+1, len(the_eigenvals)):
        ckl_k = []
        for kk in range(len(the_ref_eigenvec)):
            ckl_l = []
            for ll in range(len(the_ref_eigenvec)):
                ckl_l.append(np.abs(np.dot(the_ref_eigenvec[kk], the_eigenvecs[ii][ll].T)))
            ckl_k.append(ckl_l)
        ckl_k = np.matrix(ckl_k)
        ckl = sorted([(ckl_k[i, j], i, j) for i in range(ckl_k.shape[0]) for j in range(ckl_k.shape[1])], reverse=True)
        the_used_rows, the_used_cols, the_selected = set(), set(), []
        the_top_n=None
        for the_value, the_i, the_j in ckl:
            if the_i not in the_used_rows and the_j not in the_used_cols:
                the_selected.append([the_i, the_j]); the_used_rows.add(the_i) ; the_used_cols.add(the_j)
                if the_top_n is not None and len(the_selected) >= the_top_n:
                    break
        the_sorted_indices = sorted(the_selected, key=lambda x: (x[0], x[1]))
        the_sorted_indices = [the_sorted_indices[i][1] for i in range(len(the_sorted_indices))]
        the_sorted_eigenvecs = np.array([the_eigenvecs[ii][i] for i in the_sorted_indices])
        the_final_eigens.append(np.array([the_eigenvals[ii][i] for i in the_sorted_indices]))
        the_final_eigenvecs.append(the_sorted_eigenvecs)
        the_ref_eigenvec = the_sorted_eigenvecs # This is the part where the new sorted eigenvectors plays a role
    return [the_final_eigens, the_final_eigenvecs]



# Comments: This function checks for orthogonality of eigenvectors and returns the new order with the eigenvectors ordered such that they are associated to the corresponding state time slice by time slice
# the_eigenvals: This is an array of the values of the eigenvals: shape [Nt, Neigens]
# the_eigenvecs: This is an array of the eigenvectors in the following shape: [Nt, N, N]
def SORTING_EIGENVECTORS_NORMALIZED_CHANGING_TSLICE(the_t0, the_eigenvals, the_eigenvecs):
    the_ref_tslice_eigenval = int((len(the_eigenvals)-the_t0)/3)
    the_final_eigens, the_final_eigenvecs = list(the_eigenvals[:the_ref_tslice_eigenval+1]), list(the_eigenvecs[:the_ref_tslice_eigenval+1])
    the_ref_eigenvec = the_eigenvecs[the_ref_tslice_eigenval]
    for vec in range(len(the_ref_eigenvec)):
        the_ref_eigenvec[vec] = np.abs(the_ref_eigenvec[vec]/np.linalg.norm(the_ref_eigenvec[vec]))
    for ii in range(the_ref_tslice_eigenval+1, len(the_eigenvals)):
        ckl_k = []
        for kk in range(len(the_ref_eigenvec)):
            ckl_l = []
            for ll in range(len(the_ref_eigenvec)):
                ckl_l.append(np.abs(np.dot(the_ref_eigenvec[kk], the_eigenvecs[ii][ll]/np.linalg.norm(the_eigenvecs[ii][ll].T))))
            ckl_k.append(ckl_l)
        ckl_k = np.matrix(ckl_k)
        ckl = sorted([(ckl_k[i, j], i, j) for i in range(ckl_k.shape[0]) for j in range(ckl_k.shape[1])], reverse=True)
        the_used_rows, the_used_cols, the_selected = set(), set(), []
        the_top_n=None
        for the_value, the_i, the_j in ckl:
            if the_i not in the_used_rows and the_j not in the_used_cols:
                the_selected.append([the_i, the_j]); the_used_rows.add(the_i) ; the_used_cols.add(the_j)
                if the_top_n is not None and len(the_selected) >= the_top_n:
                    break
        the_sorted_indices = sorted(the_selected, key=lambda x: (x[0], x[1]))
        the_sorted_indices = [the_sorted_indices[i][1] for i in range(len(the_sorted_indices))]
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
    rebinned_list = []
    len_a = len(a_list)
    for hh in range(0,len_a, bin_size):
        if hh<len_a and len(rebinned_list)<int(len_a / bin_size):
            rebinned_list.append(np.mean(a_list[hh:hh+bin_size], dtype=np.float128))
    return np.array(rebinned_list, dtype=np.double)



## ------------------- REWIGHTING FACTORS --------------------------------
# Comments
# rw: list of reweighting factors [Ncfgs]
# it returns a list with the normalized reweighting factors
def RW_NORMALIZATION(rw, nfs):
    len_rw = np.double(nfs)#len(rw))
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
    da_rw_corr = []
    nt_corr = da_corr.shape[1]
    nfs = da_corr.shape[0]
    for tt in range(nt_corr):
        da_rw_corr_nt = []
        for nf in range(nfs):
            da_rw_corr_nt.append(np.double(da_corr[nf][tt]) * np.double(rw[nf]))
        da_rw_corr.append(np.array(da_rw_corr_nt,dtype=np.double))
    return np.array(da_rw_corr)

### ------------- RESAMPLING ------------------------------------------------------

# Comments:
# a_list: is the list of correlators that needs to be bootsrapped, shape [Ncfgs]
# c_conf: is the list generated to choose those specific configs, shape [K samples, Ncgfs]
# dis_rw: is the list of normalized rebinned reweights, shape [Ncfgs]
# It normalizes the reweights chosen, and calculates de average over the samples, mutplied by this factor
# It returns the list with the new samples, it has the shape [K bt samples]
def BOOTSTRAP(a_list, c_conf, dis_rw):
    bt_length = len(c_conf)
    bt_corr = []
    for k_conf in range(bt_length):
        bt_corr_k = []; rw_bt_k = []
        for cc in range(len(a_list)):
            k_th = int(c_conf[k_conf][cc])
            bt_corr_k.append( np.double(a_list[k_th]) * np.double(dis_rw[k_th]))
            rw_bt_k.append(np.double(dis_rw[k_th]))
        k_norm_factor = NORM_FACTOR(rw_bt_k)
        k_mean_corr = np.double(np.mean(bt_corr_k,dtype=np.float128))
        bt_corr.append(k_mean_corr / k_norm_factor)
    return np.array(bt_corr, dtype=np.double)


# Comments:
# it receives a list of gauge configs: [Ncfgs]
# it returns a list of configs averaged
def JACKKNIFE(a, dis_rw):
    corr_jack=[]
    for j in range(len(a)):
        corr=[]; new_rw = []
        for i in range(len(a)):
            if i!=j:
                corr.append(np.double(a[i]) * np.double(dis_rw[i]))
                new_rw.append(dis_rw[i])
        j_norm_factor = np.double(NORM_FACTOR(new_rw))
        corr_jack.append(np.double(np.mean(corr,dtype=np.float128)) * j_norm_factor)
    return np.array(corr_jack,dtype=np.double)


### ------------- STATISTICS -------------------------------------------------
#  MEANS 
# Comments:
# a:  list of [time slices (nt), nr configs (Nf)]
# it returns an array of nt slice with a mean value for each, shape [nt slices]
def MEAN(a):
    tt=0; mean_val = []
    for tt in range(len(a)):
        mean_val.append(np.double(np.mean(a[tt], dtype=np.float128)))
    return np.array(mean_val,dtype=np.float128)

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
    sigma=[]
    for ii in range(len(b)):
        s=[]
        for jj in range(len(a[0])):
            s.append(np.double((np.double(a[ii][jj]) - np.double(b[ii]))**2))
        s = np.sum(s)
        sigma.append(np.sqrt(s*pre_factor))
    return sigma

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


# STATISTICAL ERROR
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
    sigma=[]
    for ii in range(len(b)):
        sigma_jk=[]
        for jj in range(len(b)):
            s=[]
            for kk in range(len(a[0])):
                s.append(( np.double(a[ii][kk]) - np.double(b[ii])) * (np.double(a[jj][kk]) - np.double(b[jj])))
            s = np.sum(s)
            sigma_jk.append(s * pre_factor )
        sigma.append(np.array(sigma_jk, dtype=np.double))
    return np.matrix(sigma)


### ------------- FITTING FUNCTIONS AND STUFF --------------------------------------------

# Comments:
# This is the function for a single exponential fit. Ae^{-E0*(nt-t0))}
# x: is the t_slices
# e0: is a list (amplitud, energy)
# *a: this is a variable size arguments, in thie case corresponds to t0.
def SINGLE_EXPONENTIAL(x,e0,*a): 
    return e0[0] * np.exp((-e0[1]) * (x - a))

# Comments:
# This is the function for a double exponential fit. Ae^{-E0*(nt-t0)}(1+ Be^{-nt*D^{2}})
# x: is the t_slices
# e0: is a list (amplitud, energy E0, Amplitude shift of energy, DeltaE**2)
# *a: this is a variable size arguments, in thie case corresponds to t0.
def DOUBLE_EXPONENTIAL(x,e0,*a):
    return e0[0] * np.exp(-(x-a) * e0[1]) * (1. + e0[2] * np.exp(-x * (e0[3]**2)))


# Comments:
# This is the function for an alternative double exponential fit. Ae^{-E0*(nt-t0)}(1+ B e^{-nt (E1-E0)})
# x: is the t_slices
# e0: is a list (amplitud, energy E0, Amplitude shift of energy, DeltaE**2)
# *a: this is a variable size arguments, in thie case corresponds to t0.
def DOUBLE_EXPONENTIAL_ALTERNATIVE(x,e0,*a):
    return e0[0] * np.exp(-(x-a) * e0[1]) * (1. + e0[2] * np.exp(-x * (e0[3]-e0[1])))


### Comments:
# This function tries to find a good guess for the fit to have a prior, so it would in principle take less time. It uses a simple polynomial fit of order 1. 
def BEST_GUESS(c,t_i,tipo_fit):
    da_fits = np.polyfit(t_i, np.log(c), deg=1)
    amp0 = np.double(np.exp(da_fits[1]))
    eng0 = np.double(-da_fits[0])
    guess = []
    guess.append(amp0)
    guess.append(eng0)
    if tipo_fit=='2':
        guess.append(np.double(0.1))
        guess.append(np.double(0.5))
    return guess


### Comments: 
# This function gets the Difference between a Chi^{2} of one time slice compared to the next time slice value of Chi^{2}. This in order to check for stability.
def DELTA_CHI(a):
    delta_chi = []
    for ii in range(len(a)-1):
        delta_chi.append(a[ii ] - a[ii + 1])
    return np.array(delta_chi)

### Comments:
# This function gets the total Chi^{2}, not only per degree of freedom. 
def TOTAL_CHI(a,b,c,nrp):
    total_a = []
    for ii in range(len(a)):
        total_a.append(a[ii] * (np.double(c[ii] - b[ii]) - np.double(nrp)))
    return np.array(total_a)


### Comments:
# This class gets the function to do the Chi^{2} fit. 
class My_Fits:
    def __init__(self, model, x, y, cov, dgof, a):
        self.model = model  # model predicts y for given x
        self.x = np.array(x, dtype=np.double)
        self.y = np.array(y, dtype=np.double)
        self.cov = np.matrix(cov, dtype=np.double)
        self.arg = a
        self.nrPars = np.double(len(dgof))

    def __call__(self, *par):  # we accept a variable number of model parameters
        ym = self.model(self.x, *par, self.arg)
        return np.dot(np.dot(self.y - ym, self.cov), self.y - ym)/np.double((np.double(len(self.x))-self.nrPars))
    
    

### ------------------  PLOTTING STUFF --------------------------------------------------


def SQUARED_MOM(the_mom_str):
    the_mod_str = list(the_mom_str.split(','))
    the_sqrd_mom = (int(the_mod_str[0][the_mod_str[0].index('(')+1:])**2) + (int(the_mod_str[1])**2) + (int(the_mod_str[2][:the_mod_str[2].index(')')])**2)
    return the_sqrd_mom
    

### Comments:
# This function receives an operator name and returns a string in a nice way to put in the plots. 
def OPERATORS_SH(operator_name):
    new_op = list(operator_name.split(' '))
    OperatorPlot = ''
    if 'GI{' not in operator_name:
        if new_op[0].lower()=='pion':
            OperatorPlot = 'P[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='kaon':
            OperatorPlot = 'k[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='nucleon':
            OperatorPlot = 'N[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='lambda':
            OperatorPlot = 'L[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='sigma':
            OperatorPlot = 'S[%s]'%new_op[-1].replace('_','')
        elif new_op[0].lower()=='xi':
            OperatorPlot = 'X[%s]'%new_op[-1].replace('_','')
    elif 'GI{' in operator_name:
        OperatorPlot = str(new_op[-2])
    return str(OperatorPlot)


### Comments:
# This function receives an operator name and returns a string in a nice way to put in the plots. 
def OPERATORS_MH(the_operator_name):
    new_op = list(the_operator_name.split(' '))
    if "CG" in the_operator_name: the_shift = 1
    else: the_shift=0
    OperatorPlot = ''
    if 'GI{' not in the_operator_name:
        if '_' in new_op[0]:
            the_hads = list(new_op[0].split('_'))
            for ii in range(1,len(the_hads)):
                the_mom = str(SQUARED_MOM(new_op[2+the_shift+((ii-1)*3)]))
                the_irrep = str(new_op[3+the_shift+((ii-1)*3)])
                the_site = str(new_op[4+the_shift+((ii-1)*3)][:-1]).replace('_','')
                if the_hads[ii].lower()=='pion':
                    OperatorPlot += 'P[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
                elif the_hads[ii].lower()=='kaon':
                    OperatorPlot += 'k['+ the_mom + '_' + the_irrep + '_' + the_site + ']'
                elif the_hads[ii].lower()=='nucleon':
                    OperatorPlot += 'N[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
                elif the_hads[ii].lower()=='lambda':
                    OperatorPlot += 'L[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
                elif the_hads[ii].lower()=='sigma':
                    OperatorPlot += 'S[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
                elif the_hads[ii].lower()=='xi':
                    OperatorPlot += 'X[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
            
        elif '_' not in new_op[0]:
            the_mom = str(SQUARED_MOM(new_op[1]))
            the_irrep = str(new_op[2][:new_op[2].index('_')]) 
            the_site = str(new_op[-1]).replace('_','')
            if new_op[0].lower()=='pion':
                OperatorPlot = 'P[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
            elif new_op[0].lower()=='kaon':
                OperatorPlot = 'k['+ the_mom + '_' + the_irrep + '_' + the_site + ']'
            elif new_op[0].lower()=='nucleon':
                OperatorPlot = 'N[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
            elif new_op[0].lower()=='lambda':
                OperatorPlot = 'L[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
            elif new_op[0].lower()=='sigma':
                OperatorPlot = 'S[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
            elif new_op[0].lower()=='xi':
                OperatorPlot = 'X[' + the_mom + '_' + the_irrep + '_' + the_site + ']'
    elif 'GI{' in the_operator_name:
        if '_' in new_op[4]:
            OperatorPlot = new_op[4]
        else:
            the_mom = str(SQUARED_MOM(new_op[2]))
            the_irrep = new_op[3]
            if '}' in new_op[4]: new_op[4] = new_op[4][:-1]
            OperatorPlot = new_op[4][:new_op[4].index('[')+1] + the_mom + '_' + the_irrep + '_' + new_op[4][new_op[4].index('[')+1:-1] + ']'
    return str(OperatorPlot)



### Comments:
# This class rewrites the names of the irreps in TeX type of text such that they can be put in the plots in a nice way.
class IrrepInfo:
     def __init__(self,nombre):
        self.name = nombre.split('_')            
        self.Name = self.name[1]
        if self.Name[-1]=='m':
            nombre_plot = self.name[1][0]+r'$^{-}_{%s}$'%self.name[1][1:-1]
        elif self.Name[-1]=='p':
            nombre_plot = self.name[1][0]+r'$^{+}_{%s}$'%self.name[1][1:-1]
        else: 
            nombre_plot = self.name[1][0]+r'$_{%s}$'%self.name[1][1:]
        self.NamePlot = nombre_plot
        self.Momentum = self.name[0]
        self.TotalMomPlot = self.name[0][0]+ r'$^{2}=%s$'%self.name[0][-1]



### Comments:
# This class gets the info of the non-interacting levels. It only works when there is a threshold of 2 states nearby. It needs modification for the 3particle threshold.
class NonInteractingLevels:
    def __init__(self,nombre,all_hads):
        self.FirstState = str(nombre[:4])
        self.SecondState =  str(nombre[4:])
        FirstRaw = 'PSQ'+ str(self.FirstState[2])+'__'+ str(self.FirstState[0])
        SecondRaw = 'PSQ'+ str(self.SecondState[2])+'__'+ str(self.SecondState[0])
        for item in all_hads:
            if FirstRaw[:4] in item and FirstRaw[-1]==item[-1]: first_name = item;break
            else: continue
        for item in all_hads:
            if SecondRaw[:4] in item and SecondRaw[-1]==item[-1]: second_name = item;break
            else: continue
        self.First = first_name
        self.Second = second_name


### Comments:
# This function writes the errors in the plot in the way of parenthesis according to the precision given.
def WRITTING_ERRORS_PLOTS(an_error, the_precision):
    an_error = str(f'{an_error:.20f}')
    the_error_string = "("
    out_precision = True
    if len(an_error)>(the_precision+2): 
        an_error=an_error[:the_precision+2]
    for ii in range(len(an_error)):
        if len(the_error_string)>=2:
            the_error_string+=an_error[ii]
        else:
            if an_error[ii]!="0" and an_error[ii]!=".": 
                the_error_string+=an_error[ii]
            else: continue
    if len(the_error_string)>3:
        the_error_string=the_error_string[:-1]
        out_precision=False
    return [the_error_string+")", out_precision]



if __name__=="__main__":
   print('Nothing to run here, unless you want to change something.')



