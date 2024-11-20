# CLS Analysis Scripts
This repository contains all the scripts I used to analyse CLS correlators in a certain format. Please look at the README files to understand each file. 

 *** THIS IS HOW THE MAIN PLOT SCRIPT WORKS. THERE ARE OTHER "README" FILES THAT EXPLAIN EACH OF THEM. ***

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
    

