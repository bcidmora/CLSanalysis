# CLS Analysis Scripts
This repository contains all the scripts I used to analyse CLS correlators in a certain format. Please look at the README files to understand each file. 


------------------------***** MAIN SCRIPTS *****-----------------------

Notice that the names written as: $NAME$ must be fully replaced $NAME = n451, for example.


----------------- SOME REQUIREMENTS ----------------------------

1. You must have all the following scripts in the same folder:
   FOR THE "main_analysis_script.py" and for "main_plot_script.py"
    - correlators_script.py
    - effective_masses_script.py
    - eigenvalues_script.py
    - fitting_script.py
    - files_$name.py
    - add_rem_operators.py
    - plot_correlators_script.py
    - plot_effective_masses_script.py
    - plot_fits_script.py
    - set_of_functions.py

3. You must have installed the following packages of python3 to make it run smoothly:
    - Numpy
    - Scipy
    - Minuit
    - h5py
    - time
    - sys
    - os
    - Matplotlib
    - PdfMerger

4. You gotta go to each script mentioned in (1.) and modify the PATH where your data will be stored, the input data, meaning the correlators that are going to be averaged and so on. (DO NOT FORGET THIS STEP! OTHERWISE NOTHING WILL WORK.)

5. How to run it: (It is pretty simple and intuitive)
    4.1. * Very First variables: ("main_analysis_script.py")
        - runCorrs = True/False:
            + This runs the analysis of the correlators, they must be already averaged over the irrep rows and over equivalent momentum. The name of the keys in the original file must have something like "PSQ0_A1um".
        - runEffMass = True/False:
            + This calculates the effective masses after averaging over the correlators. For the correlator matrices, it obtains the effective masses of the diagonals if there is no eigenvalues yet. The correlators part must have been ran once at least before. 
        - runEigenvals = True/False:
            + This computes the eigenvalues of the correlation matrices using GEVP. It asks for a tmin and tmax for which the GEVP will be done, this is t0min/t0max.
        - runRowsCols: True/False
            + This will do the operators analysis. One can choose the operators from a list defined in the files file_$name.py or adding/removing operators from the original basis.
        - runFits = True/False:
            + This does the fits of the correlators for whatever form you pick 1/2 exponentials.
        * Very First variables: ("main_plot_script.py")
        - plotCorrs: True/False
           + This variable defines the plotting of the correlators.
        - plotEffMass: True/False
           + It plots and saves the effective masses
        - plotFits: True/False
           + It plots the fits for a chosen time interval that must be defined in "file_$name.py"
   
    4.2. Then you have the variables that are input in the Terminal:
        - myEns: The ensemble to be analysed (n451, n201, etc.)
        - myWhichCorrelator: Which correlator (single=s ; multihadron=m ; multihadron ratios=mr )
        - myTypeRs: Type of resampling (Jackknife=jk or Bootstrap=bt)
        - myRebinOn: if you want to do rebin of the data=rb, otherwise=nr
        - myLocation: string. The location where the output will be saved. Try to chose a common folder with the raw data.
        
    4.3 You have variables that are just other numbers chosen by default: (YOU CAN CHANGE THEM)
        - myRb: Size of the rebin, this is taken equals to 1 if you choose nr, otherwise you can change it here. (This must be an integer number)
        - myVersion: By default is '_test', but it serves to compare different versions when you make changes and stuff. (It is a string)
        - myKbt: Amount of Bootstrap samples by default in case you choose bt, otherwise it has no use. (It must be an integer number)
        - myNrIrreps: None/1/2... This is the number of irreps one wants to analysis, it starts from the zero-th until this number.
        - myFirstIrrep: None/1/2... This is the first irrep to analyse in case you don't want to start from the very first one in the file.
        - myLastIrrep: None/1/2... This is the last irrep to analyse in case you don't want to finish at the very end of the file list.
        - myTypeFit: '1' for one-exponential fit; '2' for a two-exponential fit and 'g' for geometric fit. (It must be a number but coded as string)
        - myTypeCorrelation: 'Correlated' if you want a correlated fit; 'Uncorrelated' if you want a uncorrelated fit. (string)
        - myOneTMin: True/False. If you want to do only one Tmin for the fit, then True.
        - myOneT0: True/False. If you want to run only one T0 in the fit, then True.
        - myT0: integer. If the above equals True, then you gotta enter which T0 you want to fit.
        - mySorting: 'eigenvals'/'eigenvals'/'vecs_fix'/'vecs_fix_norm'/'vecs_var'/'vecs_var_norm' this is how the eigenstates are sorted. Once can use the eigenvectors or the eigenvalues.
        - myKbtSamples: If you have a list of random numbers to do the Bootstrap, the file goes here. Else, keep it as None.
        - myEffMassDistance: This is the distance between data to get the effective masses.
        - myOperatorAnalysisMethod: 'from_list'/'adding'/'removing'. These are the options if you choose from a list (included in "file_$name.py") or adding and removing operators from the original basis.
        - myLocation: directory of the analysed data.
        - myWeight: reweighting factors
        - myCnfgs: This can be modified if you don't want to do all the configs
       
    4.4. Finally it prints the location of the output file so you do not lose it.



---------------------------***** CORRELATORS SCRIPT *****--------------------------

THIS IS A BRIEF EXPLANATION OF HOW THE CORRELATORS SCRIPT WORKS

FUNCTIONALITY:

This script is the main code to do the following:
    1. Average the correlators over the nr. of gauge configurations.
    OBS: The averages ALWAYS correspond to the mean value of the original data, which will later be called "central values".
    2. Implements the reweighting factors.
    3. Does binning if you want to do it.
    4. Does resampling (Jakcknife jk or Bootstrap bt).
    5. Calculates the variances and covariance matrices according to the resampling scheme. 
    OBS: The Covariances and variances are calculated with respect to the resampling data, which is not unbiased. (CHECK CHAPTER 4.5 Gattringer & Lang, "Quantum Chromodynamics on the Lattice".)
    OBS: The Covariance and variance are not calculated at correlator level of the multi hadron correlators, because what it's important here is the eigenvalues. 
    6. Saves everything in a new hdf5 file.

NOTE: This file contains two functions: One for the single hadrons (SingleCorrelatorAnalysis) and one for the Multiple hadrons (MultiCorrelatorAnalysis). This was made in order to reduce errors, so the correlators can be treated separately. Although, the functions are the same, as they do the same and use the same routines, but one implemented in a difference size array. They both received the same input parameters.

NOTE: You can run this script by itself, or it can be run along with other parts of the data analysis in the script called "main_analysis_script.py".

To Run this script alone:
    1. Check that all the information in the part "__main__" of the code is correct. Meaning, that your file with your ensemble is imported if needed. (myEns)
    2. Change the directory where the file output will be saved, this is changed in the scripts "file_ens.py". (myLocation)    
    3. Modify the version of the analysis, so you can make diferent version with different choices. By default is '0' (myVersion)    
    4. Modify which number of rebin you want. By default is 1==no-binning. (myRb)    
    5. Modify the variable "myKbt", which is the k-number of Bootstrap samples you want. By default is myKbt=500.    
    6. Go to the Terminal, and to the folder where this script is and type:
    blah@blah~ python3 (or python, whatever name you have for it) correlators_script.py $Ensemble of interest$ $s/m$(single or multihadron) $jk/bt$(Jackknife or Bootstrap) $rb/nr$(rebin or no-rebin) 
    7. It will be shown in the terminal where the file was saved.    
    8. If you do not run this file by itself, all these parameters must be included in the "main_analysis_script.py".

DESCRIPTION OF THE FUNCTIONS: 
    1.1. the_archivo: hdf5 File, contains the correlation functions data. (MANDATORY)
    1.2. the_location: String. Is the location where the output file will be saved.  Remember this location will be used later for effective masses, fits, plots. (MANDATORY)
    1.3. the_version: string. Anything you want, it is only to differentiate the versions when you want to do changes. (MANDATORY)
    1.4. the_type_rs: string ('bt' or 'jk') is the type of resampling you want. (MANDATORY)
    1.5. the_irreps: list. List of irreps to be analysed, this comes from the script: "files_ens.py". (MANDATORY)
    1.6. the_weight: array. It contains the weights of the gauge configs. (Based on a txt file.) (MANDATORY)
    1.7. kwargs: 
        1.7.1. rebin_on: string ('rb'=rebin or 'nr'=no rebin) This is to rebin data, if nothing is included here, then by default there is no rebining. (OPTIONAL)
        1.7.2. rb: integer (e.g. rb=10 ) This is mandatory oif you choose rebin_on='rb'.
        1.7.3. nr_irreps: integer. Number of irreps you want to analyse in case you do not want a full analysis. (OPTIONAL), if nothing entered, then it goes over all the data.
        1.7.4. kbt: integer. This is the amount of Bootstrap samples if you choose the_type_rs=='bt', if nothing entered here, it does 500 bootstrap samples. You can change this number here. ("OPTIONAL")
        1.7.5. own_kbt_list: array with your own random numbers for the Bootstrap. (OPTIONAL)



------------------------***** EFFECTIVE MASSES SCRIPT *****------------------------

THIS IS A BRIEF EXPLANATION OF HOW THE EFFECTIVE MASS SCRIPT WORKS

FUNCTIONALITY:

This script contains two functions: one for the single hadrons and one for the multiple hadrons. 
    1. SINGLE HADRONS: 
        1.1. It calculates the effective masses of the central value and of the resamples. 
        1.2. It also calculates the variance to plot the errors later (the variance depending on the type of resampling). This part is pretty easy and straightforward. 
        1.3. It saves the info in the same file than that averaged correlators are.
    2. MULTIHADRONS: 
        2.1. This calculates the effective masses of the diagonal of the correlators and their variance.
        2.2. It also calculates the effective masses of the eigenvalues of the correlators and their corresponding variances. If there is no eigenvalues, then it skips this part. 
        2.3. It saves the info in the same file than the averaged correlators and eigenvalues are.

NOTE: This script can be run alone or within the "main_analysis_script.py". If you run it alone, then you must consider changing the variables.
1. The input in the terminal are:
    - myEns: name of the ensemble ( e.g.: n451)
    - myWhichCorrelator: s=single hads or m=multihads
    - myTypeRs: Jackknife=jk or Bootstrap=bt
    - myRebinOn: rebin=rb if you want to rebin the data, or nr=no bining.
    
2. The info you gotta modify in the script itself
    - myRb: number of bins.
    - myVersion: String of the version you want to use. (e.g.: '0')
    - myLocation: You have to modify '$YOUR_PATH_TO_THE_AVERAGED_CORRELATORS$', this is where your output will go.

3. It prints in the screen the location of the output file.

DESCRIPTION OF THE FUNCTIONS:
    1.1. the_single_correlator_data: hdf5 File. Contains the info of the averaged correlators. 
    1.2. the_type_rs: string ('bt' or 'jk'). This is te type of resampling done to calculate the statistical errors.




---------------------------***** EIGENVALUES SCRIPT *****--------------------------

THIS IS A BRIEF EXPLANATION OF HOW THE EIGENVALUES SCRIPT WORKS

FUNCTIONALITY:

This script is the main code to do the following:
    1. It obtains the eigenvalues for a certain reference value of T0, for each tslice. 
    2. It saves the eigenvalues and the eigenvectors.
    3. It does this for the central value and for the bootstrap samples.
    4. It calculates the covariance matrices to later use them in the fittings.

NOTE: You can run this script by itself, or it can be run along with other parts of the data analysis in the script called "main_analysis_script.py".

To Run this script alone:
    
1. The input in the terminal are:
    - myEns: name of the ensemble ( e.g.: n451)
    - myTypeRs: Jackknife=jk or Bootstrap=bt
    - myRebinOn: rebin=rb if you want to rebin the data, or nr=no bining.
    
2. The info you gotta modify in the script itself
    - myRb: number of bins.
    - myVersion: String of the version you want to use. (e.g.: '0')
    - myLocation: You have to modify '$YOUR_PATH_TO_THE_AVERAGED_CORRELATORS$', this is where your output will go.

3. Additional Inputs:
    - myT0Min: Min t0 to use as a reference time slice for the GEVP.
    - myT0Max: Max t0 to use as a reference time slice for the GEVP.

    
DESCRIPTION OF THE FUNCTIONS:
    1.1. the_matrix_correlator_data: hdf5 File. Averaged correlators. (MANDATORY)
    1.2. the_type_rs: string ('bt' or 'jk'). Type of resampling you did before. 
    1.3. kwargs:
        1.3.1. t0_min: min T0 to do the GEVP. (MANDATORY) 
        1.3.2. t0_max: max T0 to do the GEVP (MANDATORY)    
    
NOTE: DO NOT FORGET TO BE CONSISTENT CHOOSING THE PATH OF YOUR OUTPUTS, BECAUSE THEY ARE REFERENCED LATER WHEN YOU WANT TO DO OTHER STEPS OR WHEN YOU DO THE PREVIOUS STEPS. THEY MUST BE CHANGED ONLY ONCE, WHEN YOU START USING THIS CODE, THEN NEVER AGAIN UNLESS YOU WANT TO STORE THEM IN A DIFFERENT PLACE.



---------------------------***** FILES SCRIPT *****--------------------------

THIS IS A BRIEF EXPLANATION OF HOW THE FILES SCRIPT WORKS

This files contains all the information to do to all the steps later. You must create one of these for each ensemble you are using. 
NOTE: There must be a separate file for the single hadrons and for the correlation matrices. They must be averaged over irrep row and over momentum. 

DESCRIPTION:
    1. location: This must be changed to the location where you are storing the files. Keep it somewhat general, because it will also be used later to store the averaged correlators.
    2. f: hdf5 File. Contains the multihadron correlation functions. 
    3. f1: hdf5 File. Contains the single hadron correlation functions data.
    4. weight_raw: 'dat' file or 'txt' file that contains the reweighting factors. Notice they must be organized such that each column has a reweighting factor and each row corresponds to a gauge config. 
    5. name: list. List of irreps contained in f.
    6. name1: list. List of irreps contained in f1.
    7. listTMaxSingleHads: array. This is the array that has the tmax for which the fitting will be done. You can look at the effective masses to change this. The order goes as the irreps go. 
    8. listTMaxMultiHads: array of arrays. This array is similar to the above one, but the tmax for each irrep includes the tmax of each of the eigenvalues.
    9. singleTMinsFitPlots: this is a list of tmins to do plot later. This should be modified once you choose which tmins are better for each hadron/irrep.
    10. multiTMinsFitPlots: array of arrays. Same than above, but a different tmin for each eigenvalue of each irrep. (Look at "files_d200.py" for an explicit example)
    11. aLat: float. This is the lattice spacing for this irrep. 
    12. betaLat: beta value for this lattice/ensemble.
    13. fmToMev: conversion from fm to MeV. It is used later for the reconstruction or to plot stuff in Energy units.
    14. LatSize: This is the Lattice size (L^{3}xT, this is L)




------------------------***** FITTING SCRIPT *****------------------------

THIS IS A BRIEF EXPLANATION OF HOW THE FITTING SCRIPT WORKS

FUNCTIONALITY:

This script contains two functions: one for the single hadrons and one for the multiple hadrons. 
    1. SINGLE HADRONS: 
        1.1. It takes the correlators and searches for a good fit, using 1 or 2 exponentials. And based on Correlated fits or Uncorrelated fits. 
        1.2. Tmin can be changed, you can calculate the fit of only one-tmin or you can go through a range, this must be modified in the variable "one_tmin = True/False".
        1.3. It takes an ansatz of solution, it can be the results from BEST_GUESS or, if nothing comes out from this, it can be the effective mass at that time slice. 
        1.4. It takes the covariance matrix for that specific time interval and minimizes the ChiSqr. 
        1.5. This is done for the central value and for the resampling. 
        1.6. Later it calculated the error for that type of resampling and saves the following info:
            - [ tmins, tmaxs, energies, sigmas of these energies, chi^2 values, sigmas of these chi^2 vals]
    2. MULTIHADRONS: 
        2.1. It does pretty much the same than for the single hadrons, but it goes through each of the eigenvalues, for each t0 found in the file.

NOTE: This script can be run alone or within the "main_analysis_script.py". If you run it alone, then you must consider changing the variables.
1. The input in the terminal are:
    - myEns: name of the ensemble ( e.g.: n451)
    - myWhichCorrelator: s=single hads or m=multihads, mr=multihadrons ratio
    - myTypeRs: jk=Jackknife or bt=Bootstrap
    - myRebinOn: rb=rebin if you want to rebin the data, or nr=no bining.
    
2. The info you gotta modify it in the script itself
    - myRb: number of bins.
    - myVersion: String of the version you want to use. (e.g.: '0')
    
    - myTypeFit: This must be 1 or 2 (amount of exponentials).
    - myTypeCorrelation: Options are 'Correlated' or 'Uncorrelated'
    - myOneTMin: True=If you only want to fit one tmin, False=all Tmins to see stability plot.
    - myOneT0: True=If you want to fit only one T0, and you must change the value of myT0. False=All T0s for stability plot.
    - myT0: Value of T0 you want to do the fit.
    - myLocation: You have to modify '$YOUR_PATH_TO_THE_AVERAGED_CORRELATORS$', this is where your input data is retrieved.
    
3. What about the functions:
    3.1. FitSingleCorrelators(the_data, the_fit_data, the_type_rs, the_list_tmaxs, **kwargs)
        - the_data: hdf5 File that contains the Correlators averaged, effective masses and eigenvalues.
        - the_fit_data: hdf5 File that will contain all the fit information. This is created if it doesnt exist, and modified if it does. 
        - the_type_rs: the type of resampling scheme you choose ('jk' or 'bt', it must be string)
        - the_list_tmaxs: array that contains all the tmaxs for which the fits will be done. This should be included in the files_ens.py file.
        - kwargs: 
            + one_tmin: True/False. (MANDATORY)
            + type_fit: '1' or '2' (MANDATORY)
            + type_correlation: 'Correlated' or 'Uncorrelated' (MANDATORY)
    
    3.2. FitMultiCorrelators(the_data, the_fit_data, the_type_rs, the_list_tmaxs, **kwargs)
        - the_data: hdf5 File, contains all the info of the correlation matrices.
        - the_fit_data: hdf5 FIle, will save the info of the fits. It will create a file if it does not exists, or modify it if it does.
        - the_type_rs: string 'jk' or 'bt'.
        - the_list_tmaxs: array with all the tmaxs for the fits.
        - kwargs:
            + one_tmin: True/False. (MANDATORY)
            + one_t0: True/False. (MANDATORY)
            + chosen_t0: integer, if one_t0=True, then this must be entered.
            + type_fit: '1' or '2' exponentials. (MANDATORY)
            + type_correlation: 'Correlated' or 'Uncorrelated'. (MANDATORY)
            + ratio_on: 'yes' or None, 'yes'=for ratios of correlators, None=normal correlator matrices. (OPTIONAL)


-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


------------------------***** PLOT MAIN SCRIPT *****------------------------

THIS IS HOW THE MAIN PLOT SCRIPT WORKS. THERE ARE OTHER "README" FILES THAT EXPLAIN EACH OF THEM.

Notice that the names written as: $NAME$ must be fully replaced $NAME$ = n451, for example.


----------------- SOME REQUIREMENTS ----------------------------

1. You must have all the following scripts in the same folder: 
    - plot_correlators_script.py
    - plot_effective_masses_script.py
    - plot_fits_script.py
    - set_of_functions.py
    - files_$name$.py

2. You must have installed the following packages of python3 to make it run smoothly:
    - PdfMerger
    - h5py
    - sys
    - os

3. You gotta go to each script mentioned in (1.) and modify the PATH where your data will be stored, the input data, meaning the correlators that are going to be averaged and so on. (DO NOT FORGET THIS STEP! OTHERWISE NOTHING WILL WORK.)

4. How to run it: (It is pretty simple and intuitive)
    4.1. Very First variables: 
        - plotCorrs = True/False: This plots only the correlation functions.
        - plotEffMass = True/False: This plots the effective masses. 
        - plotFits = True/False: This plots the fits over tmins. Notice you gotta run it only if you have more than 1 tmin in the fits. 
        - joinPlots = True/False: This puts all the plots together in one PDF file by irrep name, and by irrep and eigenvalue. 
        
    4.2. Then you have the variables that are input in the Terminal:
        - myEns: The ensemble to be analysed (n451, n201, etc.)
        - myWhichCorrelator: Which correlator (single=s ; multihadron=m ; multihadron ratios=mr )
        - myTypeRs: Type of resampling (Jackknife=jk or Bootstrap=bt)
        - myRebinOn: if you want to do rebin of the data=rb, otherwise=nr
        
    4.3 You have variables that are just other numbers chosen by default: (YOU CAN CHANGE THEM)
        - myRb: Size of the rebin, this is taken equals to 1 if you choose nr, otherwise you can change it here. (This must be an integer number)
        - myVersion: By default is '0', but it serves to compare different versions when you make changes and stuff. (it is a string)
        - myNrExponentials: string ('1' or '2'). Number of exponential in the fits.
        - myTypeCorrelation: string ('Correlated' or 'Uncorrelated'). Type of correlated fit.
        - myOneTmin: If only onetmin was fitted, then DeltaChiQS will not exist. 
        - myT0: Which T0 will be mainly plotted.         
        - myDataLocation: string. The location where the data is saved.
        - myPlotLocation: The location of the output plots.
        
    4.4. Finally it prints the location of the output files so you do not lose them.
    


------------------------***** PLOT CORRELATORS SCRIPT *****------------------------

THIS IS A BRIEF EXPLANATION OF HOW THE PLOT CORRELATORS SCRIPT WORKS

This script has the functions to plot over all irreps contained in the averaged correlators file. 

1. Variables form the Terminal:
    1.1. myEns: string. The ensemble you want to plot (e.g. n451)
    1.2. myWhichCorrelator: string (single hadrons=s ; multihadrons=m ; ratio multihadrons=mr )
    1.3. myTypeRs: string ('bt' or 'jk'). Type of resampling.
    1.4. myRebinOn: string (rebinned='rb' or no rebinned='nr'). Plot rebinned data or non-rebinned.
2. Variables in the code:
    2.1. myRb: integer. If above equals myRebinOn='rb'.
    2.2. myVersion: string. The name of the version you want to plot.
    2.3. myT0: integer. which t0 you want to plot in case of multihadrons or ratio of correlators.
    2.4. myDataLocation: string. Checks the directory where the data to bt plotted is stored. 
    
DESCRIPTION OF THE FUNCTIONS:
    - the_single_correlator_data: hdf5 File. It contains the correlation functions. (MANDATORY)
    - the_type_rs: string ('bt' or 'jk'). Type of resampling to include in the figure.  (MANDATORY)
    - the_version: string. Version of the plotted data, to include in the name of the output figures. (MANDATORY)
    - the_location: string. Where it is going to be stored. (MANDATORY)
    - the_rebin: string. To include in the name of the plots. (MANDATORY)
    It produces several plots:
        + Correlators/Eiganevalues versus lattice time.
        + Log plot Correlators/Eigenvalues versus time.
        + Histogram of distribution of the resamples around the mean value. Only for one tslice.
    



---------------------***** PLOT EFFECTIVE MASSES SCRIPT *****---------------------

THIS IS A BRIEF EXPLANATION OF HOW THE PLOT EFFECTIVE MASSES SCRIPT WORKS

This script has the functions to plot over all irreps contained in the averaged correlators file. 

1. Variables form the Terminal:
    1.1. myEns: string. The ensemble you want to plot (e.g. n451)
    1.2. myWhichCorrelator: string (single hadrons=s ; multihadrons=m ; ratio multihadrons=mr )
    1.3. myTypeRs: string ('bt' or 'jk'). Type of resampling.
    1.4. myRebinOn: string (rebinned='rb' or no rebinned='nr'). Plot rebinned data or non-rebinned.
2. Variables in the code:
    2.1. myRb: integer. If above equals myRebinOn='rb'.
    2.2. myVersion: string. The name of the version you want to plot.
    2.3. myT0: integer. which t0 you want to plot in case of multihadrons or ratio of correlators.
    2.4. myDataLocation: string. Checks the directory where the data to bt plotted is stored. 
    
DESCRIPTION OF THE FUNCTIONS:
    - the_single_correlator_data: hdf5 File. It contains the correlation functions. (MANDATORY)
    - the_rs_scheme: string ('Bootstrap' or 'Jackknife'). Type of resampling scheme to print in the figure. (MANDATORY)
    - the_version: string. Version of the plotted data, to include in the name of the output figures. (MANDATORY)
    - the_location: string. Where it is going to be stored. (MANDATORY)
    - the_rebin: string. To include in the name of the plots. (MANDATORY)
    It produces one plot per each irrep:
        + Effective Masses of Correlators/Eigenvalues versus lattice time.
            



-----------------------------***** PLOT FITS SCRIPT *****------------------------

THIS IS A BRIEF EXPLANATION OF HOW THE PLOT FITS SCRIPT WORKS

This script has the functions to plot over all irreps contained in the fits of correlators file. 

1. Variables from the Terminal:
    1.1. myEns: string. The ensemble you want to plot (e.g. n451)
    1.2. myWhichCorrelator: string (single hadrons=s ; multihadrons=m ; ratio multihadrons=mr )
    1.3. myTypeRs: string ('bt' or 'jk'). Type of resampling.
    1.4. myRebinOn: string (rebinned='rb' or no rebinned='nr'). Plot rebinned data or non-rebinned.
2. Variables in the code:
    2.1. myTypeFit: string ('Correlated' or 'Uncorrelated'). Type of correlated fit you want to plot. 
    2.2. myNrExponentials: string ('1' or '2'). Number of exponential fit you want to plot.
    2.3. myRb: integer. If myRebinOn='rb'.
    2.2. myVersion: string. The name of the version you want to plot.
    2.3. myT0: integer. which t0 you want to plot in case of multihadrons or ratio of correlators.
    2.4. myDataLocation: string. Checks the directory where the data to bt plotted is stored. 
    
DESCRIPTION OF THE FUNCTIONS:
    - the_single_fit_data: hdf5 File. It contains the fits to the correlation functions. (MANDATORY)
    - the_type_fit: string (myTypeFit). Correlated or Uncorrelated. (MANDATORY)
    - the_nr_exps: string (myNrExponentials = '1' or '2'). Number of exponentials. (MANDATORY)
    - the_tmins: list. List of tmins to do the plotting. (MANDATORY)
    - the_version: string. Version of the plotted data, to include in the name of the output figures. (MANDATORY)
    - the_location: string. Where it is going to be stored. (MANDATORY)
    - the_rebin: string. To include in the name of the plots. (MANDATORY)
    It produces several plots per each irrep:
        + Energy fits versus tmins.
        + Zoom to Energy fits versus different tmins.
        + Chi^2 versus tmins.
        + Total Chi^2 versus tmins.
        + Delta Chi^2 versus tmins. This is to check how much the Chi^2 changes from tmin to tmin.
        + For Multihadrons it also plots: 
            ++ Energy versus T0s to check stability. Only if more than 2 T0s found in data file.
