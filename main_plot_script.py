import plot_correlators_script as pcorr
import plot_effective_masses_script as peff
import plot_fits_script as pfit
import set_of_functions as vf
import sys
import os
import h5py
from PyPDF2 import PdfMerger

### ------- WHAT IT IS DONE --------

plotCorrs = True
plotEffMass = True
plotFits = True
joinPlots = True


### ------ MAIN VARIABLES ---------

myEns = str(sys.argv[1]).upper()
myWhichCorrelator = str(sys.argv[2]).lower()
myTypeRs = str(sys.argv[3]).lower()
myRebinOn = str(sys.argv[4]).lower()

myRb = 2
myVersion = 'test'
myNrExponentials = '1'
myTypeCorrelation = 'Correlated'
myOneTmin = True
myT0 = 4

myDataLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH(SAME_THAN_CORRS_SCRIPT_OUTPUT)$/%s/'%myEns)


if myEns == 'N451': from files_n451 import singleTMinsFitPlots, multiTMinsFitPlots
elif myEns == 'N201': from files_n201 import singleTMinsFitPlots, multiTMinsFitPlots 
elif myEns == 'D200': from files_d200 import singleTMinsFitPlots, multiTMinsFitPlots
    

if myRebinOn=='rb': 
    reBin = '_bin'+str(myRb)
else: 
    reBin = ''


if myTypeRs=='jk':
    myResamplingScheme='Jackknife'
elif myTypeRs=='bt':
    myResamplingScheme='Bootstrap' 



### -------- PRINTING INFO OF ENSEMBLE ---------

vf.INFO_PRINTING(myWhichCorrelator, myEns)



### ------------ START ----------------

if myWhichCorrelator =='s':
    mySingleCorrelatorData = h5py.File(myDataLocation + 'Single_correlators_' + myTypeRs + reBin + '_v%s.h5'%myVersion,'r')
    myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/SingleHadrons/'%myEns +  '%s/'%myResamplingScheme)
    
    if plotCorrs:        
        pcorr.PlotSingleHadronCorrelators(mySingleCorrelatorData, myTypeRs, myVersion, myPlotLocation, reBin)
    
    if plotEffMass: 
        peff.PlotSingleHadronsEffectiveMasses(mySingleCorrelatorData, myResamplingScheme, myVersion, myPlotLocation, reBin)
    
    if plotFits: 
        myFitsLocation = vf.DIRECTORY_EXISTS(myDataLocation + 'Fits_SingleHadrons/')
        myFitCorrelator =  h5py.File(myFitsLocation + 'Single_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        pfit.PlotSingleHadronsFits(myFitCorrelator, myTypeCorrelation, myNrExponentials,  singleTMinsFitPlots, myVersion, myPlotLocation, reBin)
        
        myFitCorrelator.close()
        
    if joinPlots:
        irreps = list(mySingleCorrelatorData.keys())
        for aa in irreps:
            ops = list(mySingleCorrelatorData[aa+'/Operators'])
            x=[]
            x.append(myPlotLocation + 'Correlator_'+ aa[:4] + '_' + aa[-1] + reBin +'_v%s.pdf'%myVersion)
            x.append(myPlotLocation + 'Correlator_'+ aa[:4] + '_' + aa[-1] + '_log' + reBin + '_v%s.pdf'%myVersion)
            x.append(myPlotLocation + 'Histogram_correlators_'+ aa[:4] + '_' + aa[-1] + reBin + '_v%s.pdf'%myVersion)
            x.append(myPlotLocation + 'EffectiveMass_'+ aa[:4] + '_' + aa[-1] + reBin + '_v%s.pdf'%myVersion)
            x.append(myPlotLocation + 'Tmin_Fits_'+ aa[:4] + '_' + aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
            x.append(myPlotLocation + 'Tmin_Fits_Zoom_'+ aa[:4] + '_' + aa[-1] +'_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
            x.append(myPlotLocation + 'Tmin_Chisqr_'+ aa[:4]+'_' +aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
            x.append(myPlotLocation + 'Tmin_TotalChisqr_'+ aa[:4]+'_' +aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
            if not myOneTmin:
                x.append(myPlotLocation +  'Tmin_DeltaChisqr_'+ aa[:4]+'_' +aa[-1] + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
            merger = PdfMerger()
            for pdf in x:
                merger.append(open(pdf, 'rb'))

            with open(myPlotLocation + "%s"%aa + "_%sexp"%myNrExponentials + reBin + "_v%s.pdf"%myVersion, "wb") as fout:
                merger.write(fout)
        print('Now all the plots are in one file')
        
    mySingleCorrelatorData.close()
    
elif myWhichCorrelator=='m':        
    myMatrixCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
    myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/Matrices/'%myEns +  '%s/'%myResamplingScheme)
    
    if plotCorrs: 
        pcorr.PlotMultiHadronCorrelators(myMatrixCorrelatorData, myTypeRs, myVersion, myT0, myPlotLocation, reBin)
    
    if plotEffMass: 
        peff.PlotMultiHadronsEffectiveMasses(myMatrixCorrelatorData, myResamplingScheme, myVersion, myT0, myPlotLocation, reBin)
        
    if plotFits:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myDataLocation + 'Fits_Matrices/')
        myFitCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        pfit.PlotMultiHadronsFits(myFitCorrelator, myTypeCorrelation, myNrExponentials, multiTMinsFitPlots, myT0, myVersion, myPlotLocation, reBin)
        
        myFitCorrelator.close()

    if joinPlots:
        irreps = list(myMatrixCorrelatorData.keys())
        for aa in irreps:
            ops = list(myMatrixCorrelatorData[aa+'/Operators'])
            x=[]
            for bb in range(len(ops)):
                x.append(myPlotLocation + 'Eigenvalues_' + aa + '_%s'%bb + reBin + '_v%s.pdf'%myVersion)
                x.append(myPlotLocation + 'Eigenvalues_' + aa +  '_%s_log'%bb + reBin + '_v%s.pdf'%myVersion)
                x.append(myPlotLocation +'Histogram_Eigenvalues_' + aa + '_%s'%bb + reBin +  '_v%s.pdf'%myVersion)
                x.append(myPlotLocation + 'EffectiveMass_Eigenvalues_'+ aa + '_%s'%bb + reBin + '_v%s.pdf'%myVersion)
                x.append(myPlotLocation + 'T0s_'+ aa + '_%s'%bb + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
                x.append(myPlotLocation + 'Tmin_Fits_'+ aa + '_%s'%bb + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
                x.append(myPlotLocation + 'Tmin_Fits_Zoom_'+ aa + '_%s'%bb + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
                x.append(myPlotLocation + 'Tmin_Chisqr_'+ aa + '_%s'%bb + '_%sexp'%myNrExponentials + reBin + '_v%s.pdf'%myVersion)
                x.append(myPlotLocation + 'Tmin_TotalChisqr_'+ aa + '_%s'%bb + '_%sexp'%myNrExponentials + reBin  + '_v%s.pdf'%myVersion)
                if not myOneTmin:
                    x.append(myPlotLocation + 'Tmin_DeltaChisqr_'+ aa + '_%s'%bb + '_%sexp'%myNrExponentials + reBin  + '_v%s.pdf'%myVersion)
            merger = PdfMerger()
            for pdf in x:
                merger.append(open(pdf, 'rb'))

            with open(myPlotLocation +"%s"%aa + "_%sexp"%myNrExponentials + reBin + "_v%s.pdf"%myVersion, "wb") as fout:
                merger.write(fout)
        print('Now all the plots are in one file')
        
    myMatrixCorrelatorData.close()
    
    
elif myWhichCorrelator=='mr':
    myRatioCorrelatorData = h5py.File(myDataLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
    myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/Matrices_Ratios/'%myEns +  '%s/'%myResamplingScheme)

    if plotCorrs: 
        pcorr.PlotMultiHadronCorrelators(myRatioCorrelatorData, myTypeRs, myVersion, myT0, myPlotLocation, reBin)

    if plotEffMass: 
        peff.PlotMultiHadronsEffectiveMasses(myRatioCorrelatorData, myResamplingScheme, myVersion, myT0, myPlotLocation,reBin)
        
    if plotFits:        
        myFitsLocation = vf.DIRECTORY_EXISTS(myDataLocation + 'Fits_Ratios/')
        myFitCorrelator = h5py.File(myFitsLocation + 'Matrix_correlators_ratios_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion, 'a')
        
        pfit.PlotMultiHadronsFits(myFitCorrelator, myTypeCorrelation, myNrExponentials,  multiTMinsFitPlots, myT0, myVersion, myPlotLocation, reBin)
        
        myFitCorrelator.close()

    if joinPlots:
        print('Now all the plots are in one file')
        
    myRatioCorrelatorData.close()

### --------- PRINTS WHERE IT IS SAVED ---------
        
print('-'*(len(myPlotLocation)+1))
print('Correlator analysis saved : \n' + myPlotLocation)
print('_'*(len(myPlotLocation)+1))

