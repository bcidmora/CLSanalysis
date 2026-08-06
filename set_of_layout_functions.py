import os
from pathlib import Path
import set_of_plot_functions as vfp
from dataclasses import dataclass
from typing import Optional, Tuple


### Here are some classes defined in order to make the selection of analysis tasks easier
@dataclass
class GEVPRuns:
    t0min: Optional[int] = None
    t0max: Optional[int] = None
    sorting: str = "eigenvals" #'vecs_fix' #'eigenvals' # 'vecs_fix' # 'vecs_fix_norm' # 'vecs_var' # 'vecs_var_norm'
    rs_sorting: Optional[str] = None
    td: Optional[int] = None
    

@dataclass
class OperatorRuns:
    t0min: Optional[int] = None
    t0max: Optional[int] = None
    sorting: str = "eigenvals" #'vecs_fix' #'eigenvals' # 'vecs_fix' # 'vecs_fix_norm' # 'vecs_var' # 'vecs_var_norm'
    rs_sorting: Optional[str] = None
    ops_list: list = None


@dataclass
class FitRuns:
    type_fit: str = '1' #'1' #'2' #'g'
    t0: int = 4
    type_correlation: str = 'Correlated' # 'Uncorrelated'
    one_tmin: bool = False # False
    one_t0: bool = True # False


@dataclass
class BinRuns:
    max_bin: int = 20
    fit_range: Tuple[int, int] = (22, 36)
    chosen_bin: int = 5
    
    
@dataclass
class DispRuns:
    mode: str = "both"  # abs / rel / both
    

@dataclass
class NrIrreps:
    first_irrep : int
    last_irrep: int
    steps: Optional[int]
    nr_irreps : int
    

@dataclass
class NrConfigs:
    first_config : int
    last_config: int
    steps: Optional[int]
    nr_configs : int
    

@dataclass
class Runs:
    ensemble: str
    correlator: str
    rs_type: str
    isospin: Optional[str]
    version: Optional[str]

    rebin: bool
    rb: int
    
    ### What to do 
    corrs: bool
    eigenvals: bool
    effmass: bool
    fits: bool
    disp: bool
    corrdiff: bool
    ops: bool
    binning: bool

    ### Subflags deriving from the main tasks
    fit: FitRuns
    gevp: GEVPRuns
    binning_run: BinRuns
    disp_run: DispRuns
    ops_run: OperatorRuns    
    
    ### These variables have a default setting
    kbt: int 
    dist_eff_mass: int
    
    ### These are intended to differentiate between the main set of correlators and the reduced set when obtianing the eff masses and performing the fits
    ops_flag: bool
    gevp_flag: bool
    
    ### Getting the irreps if not all are wanted
    the_irreps: NrIrreps
    the_configs: NrConfigs
    

import argparse
### Comments:
# this function allows you to put all inputs in the terminal as variables. Easier for the amount of variables we have
def parse_args():
    parser = argparse.ArgumentParser(description="Correlator analysis")

    ### These are mandatory
    parser.add_argument("-e", "--ensemble", required=True) # lowercase ensemble name
    parser.add_argument("-c", "--correlator", required=True, choices=["s", "m", "mr"]) # single-, multi-hadron or ratio of corrs.
    parser.add_argument("-rs", "--rs-type", required=True, choices=["jk", "bt"]) # resampling schemes    

    ### These are optional
    parser.add_argument("-i", "--isospin", default=None, choices=["s", "d", "t", "q", "dk", "dd", "dd*"]) # isosinglet, isodoublet, isotriplet, isoquartet, isodoublet k-meson, isodoublet d-meson, isodoublet d-star-meson
    parser.add_argument("-v", "--version", type=str, default='test') # extension of the file's name
    
    ### Binning info
    parser.add_argument("--rebin", action="store_true")
    parser.add_argument("-rb", "--rb", type=int, default=10)

    ### What to do
    parser.add_argument("--corrs", action="store_true") # Correlator analysis
    parser.add_argument("--eigenvals", action="store_true") # GEVP
    parser.add_argument("--effmass", action="store_true") # Effective masses
    parser.add_argument("--fits", action="store_true") # Fitting procedure
    parser.add_argument("--disp-rel", action="store_true") # Dispersion relations
    parser.add_argument("--corr-diff", action="store_true") # Correlated differences
    parser.add_argument("--ops", action="store_true") # Operator analysis
    parser.add_argument("--binning", action="store_true") # Binning analysis
    
    ### These are for the GEVP and ops analysis
    parser.add_argument("--t0min", type=int)
    parser.add_argument("--t0max", type=int)    
    
    ### these are for the dispersion relations
    parser.add_argument("--disp-mode", default="both", choices=["abs", "rel", "both"])
    
    ### These are for the effective masses
    parser.add_argument("--ops-flag", action="store_true")
    parser.add_argument("--gevp-flag", action="store_true")
    
    parser.add_argument("--dist-eff-mass", type=int, default=1)
    
    ### Fit params 
    parser.add_argument("--fit-type", choices=["1", "2", "g"], default="1")
    parser.add_argument("--fit-correlation", choices=["Correlated", "Uncorrelated"], default="Correlated")
    parser.add_argument("--fit-t0", default=4)
    parser.add_argument("--fit-one-tmin", choices=[True, False], default=False)
    
    ### These is for the Bootstrap
    parser.add_argument("-kbt", "--k-bootstrap", type=int, default=500) # resampling schemes
    
    ### How many Irreps to do
    parser.add_argument("-fi", "--start-irrep", type=int)
    parser.add_argument("-li","--last-irrep", type=int)
    parser.add_argument("-ir","--nr-irreps", type=int)
    
    ### How many configs to do
    parser.add_argument("-fc", "--start-config", type=int)
    parser.add_argument("-lc","--last-config", type=int)
    parser.add_argument("-nc","--nr-configs", type=int)
    
    return parser.parse_args()


### Commments:
# This function checks if all the needed parameters were included before running anything
def VALIDATE_RUNS(r: Runs):
    if r.correlator in ("m", "mr") and r.isospin is None:
         raise ValueError("Multihadron analysis requires -i or --isospin.")
     
    if r.eigenvals or r.ops:
        if r.gevp.t0min is None or r.gevp.t0max is None:
            raise ValueError("GEVP/OPS requires --t0min and --t0max")

    if r.binning and r.binning_run is None:
        raise ValueError("Binning requires binning config")

    if r.disp and r.disp_run is None:
        raise ValueError("Plotting dispersion mode missing --disp-mode")
    
    

### Comments:
# This routine selects which runs to do, for example the corrs, the effective masses, etc.
def WhichRuns(args, the_ensemble_data):
    
    ### Fit parameters
    fit_run = FitRuns(type_fit = args.fit_type if args.fits else "1", type_correlation = args.fit_correlation, t0=args.fit_t0, one_tmin = args.fit_one_tmin)
    
    ### GEVP Parameters
    gevp_run = GEVPRuns(t0min=args.t0min if (args.eigenvals or args.ops) else None,t0max=args.t0max if (args.eigenvals or args.ops) else None,)
    ops_run = OperatorRuns(t0min=args.t0min if (args.eigenvals or args.ops) else None,t0max=args.t0max if (args.eigenvals or args.ops) else None, ops_list = the_ensemble_data.get("operatorsChoice", None))
    
    ### Other analyses
    bin_run = BinRuns()
    disp_run = DispRuns(mode=args.disp_mode if args.disp_rel else "both")

    ### Bootstrap parameters
    kbt = args.k_bootstrap if args.rs_type == "bt" else 500
    
    ### Which irreps to run
    the_irreps = NrIrreps(first_irrep = args.start_irrep if args.start_irrep else None, last_irrep = args.last_irrep if args.last_irrep else None, nr_irreps =  args.nr_irreps if args.nr_irreps else None, steps = 1)
    the_configs = NrConfigs(first_config = args.start_config if args.start_config else 0, last_config = args.last_config if args.last_config else -1, nr_configs = args.nr_configs if args.nr_configs else None, steps = 1)

    return Runs(
        ensemble=args.ensemble.upper(),
        correlator=args.correlator.lower(),
        rs_type=args.rs_type,
        isospin=args.isospin.lower() if args.isospin else None,
        version=args.version,

        rebin=args.rebin,
        rb=args.rb,

        corrs=args.corrs,
        eigenvals=args.eigenvals,
        effmass=args.effmass,
        fits=args.fits,
        disp=args.disp_rel,
        corrdiff=args.corr_diff,
        ops=args.ops,
        binning=args.binning,

        fit=fit_run,
        gevp=gevp_run,
        ops_run=ops_run,
        binning_run=bin_run,
        disp_run=disp_run,
        
        kbt=kbt, 
        dist_eff_mass=args.dist_eff_mass,
        
        gevp_flag=args.gevp_flag,
        ops_flag=args.ops_flag,
        
        the_irreps = the_irreps,
        the_configs = the_configs,
        )
    
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
# This function prints the information of the ensemble and the type of correlators one have.
def INFO_PRINTING(the_corr_type, the_ensemble):
    if the_corr_type=='s':
        print('.............................................................................')
        print('                         SINGLE HADRON CORRELATORS')
        print(f'                               ENSEMBLE {the_ensemble}')
        print('.............................................................................')
    elif the_corr_type=='m':
        print('.............................................................................')
        print('                         MULTIHADRON CORRELATORS')
        print(f'                               ENSEMBLE {the_ensemble}')
        print('.............................................................................')
    elif the_corr_type=='mr':
        print('.............................................................................')
        print('                       MULTIHADRON RATIO CORRELATORS')
        print(f'                               ENSEMBLE {the_ensemble}')
        print('.............................................................................')
        
    

### Comments:
# This function is intended to provide a meny with the choices of hadrons to do the dispersion relation plots
def HADRON_MENU(the_irreps):
    the_rest_irreps = [the_irrep for the_irrep in the_irreps if the_irrep.startswith("PSQ0_")]
    the_hadrons = sorted({vfp.IrrepInfo(the_irrep).Hadron for the_irrep in the_rest_irreps})
    
    if len(the_hadrons)==1:
        print(f"Available hadrons: {the_hadrons}")
        return the_hadrons
    else:
        print("Available hadrons:")
        for ii, hh in enumerate(the_hadrons):
            print(f"[{ii+1}] {hh}")
        print(f"[{len(the_hadrons)+1}] All")
            
        the_choice = input("Select hadrons (separated by commas): ")
        the_indices = [int(the_x.strip()) for the_x in the_choice.split(",")]
        
        if (len(the_hadrons) + 1) in the_indices:
            return the_hadrons
        else:
            return [the_hadrons[ii-1] for ii in the_indices]


### Comments:
def SINGLE_INFO_PRINTING(the_irreps, bin_analysis, corr_analysis):
    info_list = []
    if bin_analysis:
        print('Choose one irrep for the bin-size analysis with PSQ = 0.')
    if corr_analysis:
        print('Choose the irreps you want to include in the correlator analysis.')
    else: print('Choose the irreps you want to include in the effective mass- or fit analysis.')
    for i in range(len(the_irreps)):
        info_list.append(f'{the_irreps[i]} = {i}')
    print(*info_list, sep='\n')


### Comments:
# prints the name of the irreps and the index of how they appear in the file as well as the operator combinations and the corresponding indices
def IRREP_OP_INFO_PRINTING(the_new_archivo):

    the_irreps = list(the_new_archivo.keys())
    info_list = []
    print('Irreps and the operator combinations you can choose:')
    for i in range(len(the_irreps)):
        irrep_info = []
        op_list = list(the_new_archivo[the_irreps[i]].keys())
        irrep_info.append(the_irreps[i] + ' = ' + str(i))
        for ii in range(len(op_list)):
            irrep_info.append('      ' + op_list[ii] + ' = ' + str(ii))
        info_list.append(irrep_info)
    columns_per_row = 4
    column_width = 60  # Adjust for spacing, here taken 60 because then long operator lists can be shown correctly
    for i in range(0, len(info_list), columns_per_row):
        chunk = info_list[i:i + columns_per_row]
        # Print headlines
        for sublist in chunk:
            print(f"{sublist[0]:<{column_width}}", end='')
        print()
        # Find the tallest sublist
        max_height = max(len(sublist) for sublist in chunk) - 1
        # Print remaining elements
        for row in range(max_height):
            for sublist in chunk:
                if row + 1 < len(sublist):
                    print(f"{sublist[row + 1]:<{column_width}}", end='')
                else:
                    print(" " * column_width, end='')
            print()
        # Print separation between blocks
        print()  # Blank line for spacing
        

### Comments:
# this functions creates a data set if it doesn exist, or it replaces if it exists
def REPLACE_DATASET(the_group, the_name, the_data):
    if the_name in the_group:
        del the_group[the_name]
    the_group.create_dataset(the_name, data=the_data)
        
        
        
### Comments:
# Prints the index and the name of the irrep in the file and reminds one to choose a rest frame irrep
def IRREP_BIN_SIZE_INFO_PRINTING(the_irreps):
    info_list = []
    print('Choose one irrep for the bin-size analysis with PSQ = 0:')
    for i in range(len(the_irreps)):
        info_list.append(f'{the_irreps[i]} = {i}')
    print(*info_list, sep='\n')
    
    
    
    
    
    
    

if __name__=="__main__":
   print('Nothing to run here, unless you want to change something.')
   

