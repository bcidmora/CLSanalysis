import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_layout_functions as vfl
import set_of_analysis_functions as vfa
import set_of_plot_functions as vfp


def PlotSingleHadronsFits(the_single_fit_data, the_type_fit, the_nr_exps, the_tmins, the_version, the_location, the_rebin, the_irreps, **kwargs):    
    
    ### These are the irreps in this file
    s_irreps = list(the_single_fit_data.keys())
    
    the_zoom_flag = kwargs.get('zoom_fit')
    the_chi_plots_flag = kwargs.get('chi_plots')
    the_total_chi_flag = kwargs.get('total_chi')
    the_delta_chi_flag = kwargs.get('delta_chi')
    the_all_fits_flag = kwargs.get('all_fits_comparison')
    
    # ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(the_irreps)
    
    s_irreps = s_irreps[the_first_irrep:the_last_irrep]
    
    ### Loop over the irreps found in the fits file.
    for the_irrep in s_irreps:
        
        print("---------------------------------------------------------------------------")
        print(f"Irrep: {the_irrep}")
        
        ### This is the info of the fits for this irrep
        dis_set =  np.asarray(the_single_fit_data[f'{the_irrep}/{the_nr_exps}exp/Tmin/{the_type_fit}/Mean'])
        
        ### Mean values of the fits
        the_fit_data = dis_set[2]
        
        ### Statistical errors of those central values
        the_fit_sigmas = dis_set[3]
        
        ### The Chi^2
        the_chi_corr = dis_set[4]
        
        ### The error of these Chi^2's
        the_chi_sigmas = dis_set[5]
        
        ### This is the range of min time slices which the fits were performed
        the_nt = [int(x) for x in dis_set[0]]
        
        ### Max time slice used for the fit
        the_nt_max = int(dis_set[1][0])
        
        ### This is only for the plots
        the_nt_ticks = np.arange(the_nt[0], the_nt[-1], int(len(the_nt)/5))
        
        ### This is the tmin chosen for this fit. It can be changed in "file_ens.py"
        the_chosen_tmin = the_nt.index(the_tmins[the_irreps.index(the_irrep)])
        
        ### This is just to write the errors properly in the plot
        the_mean_fit_string = str(np.round(the_fit_data[the_chosen_tmin], 5))
        the_error_string = vfp.WRITTING_ERRORS_PLOTS(the_fit_sigmas[the_chosen_tmin],5)
        the_sigmas_fit_string = the_error_string[0]
            
        if the_error_string[1]==False:
            the_mean_fit_string = str(f'{np.round(the_fit_data[the_chosen_tmin], the_error_string[2]):.{the_error_string[2]}f}')
        
        ### The SH operator that appears in the plot
        the_op = list(the_single_fit_data[f'{the_irrep}/Operators'])[0]
        theOperatorNamePlot = vfp.OPERATORS_SH(the_op.decode('utf-8'))
        
        da_irrep = vfp.IrrepInfo(the_irrep)
        MomentumIrrep = da_irrep.TotalMomPlot
        nrMomentumIrrep = da_irrep.Momentum
        NameIrrepPlot = da_irrep.NamePlot 
    
        OperatorNamePlot = vfp.SH_OPERATORS_RELABEL(theOperatorNamePlot, NameIrrepPlot, nrMomentumIrrep )
        
        print('Fits vs tmin plot in progress...')
        # fit_fig = plt.figure(figsize=(5.5,4.5))
        fit_fig = plt.figure()        
        the_label = r'${[t_{\mathrm{min}}}, t_{\mathrm{max}} ]= [%sa,$'%str(int(the_nt[the_chosen_tmin])) +r'$%sa]$'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/\mathrm{d.o.f} = %s$'%np.round(the_chi_corr[the_chosen_tmin],3) + '\n' + r'$E_{\mathrm{fit}} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
        
        vfp.PLOT_FITS(the_nt, the_fit_data, the_fit_sigmas, the_chosen_tmin, the_label, r'${t_{\mathrm{min}}/\,a}$', r'$a \; E_{\mathrm{lab}}$', f'{OperatorNamePlot} ({the_nr_exps}-exp)', the_nt_ticks)
        plt.show()
        fit_fig.savefig(f'{the_location}Tmin_Fits_{the_irrep[:4]}_{the_irrep[-1]}_{the_nr_exps}exp{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        
        if the_zoom_flag: 
            print('Zoom Fits vs tmin plot in progress...')
            zoom_fit_fig = plt.figure()  #plt.figure(figsize=(5.5,4.5))
            vfp.PLOT_FITS(the_nt, the_fit_data, the_fit_sigmas, the_chosen_tmin, the_label, r'${t_{\mathrm{min}}}$', r'$a_{t} \; E_{\mathrm{lab}}$', fr'{OperatorNamePlot} ({the_nr_exps}-exp)', the_nt_ticks, zoom=True, the_ll =3, the_ul=4)
            plt.show()
            zoom_fit_fig.savefig(f'{the_location}Tmin_Fits_Zoom_{the_irrep[:4]}_{the_irrep[-1]}_{the_nr_exps}exp{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        
        if the_chi_plots_flag:
            print('Chi2 vs tmin plot in progress...')
            chi_fig = plt.figure()  #plt.figure(figsize=(5.5,4.5))
            vfp.PLOT_FITS(the_nt, the_chi_corr, the_chi_sigmas, the_chosen_tmin, the_label, r'${t_{\mathrm{min}}}\,/\,a$', r'$\chi^{2}/\mathrm{d.o.f} $', f'{OperatorNamePlot} ({the_nr_exps}-exp)', the_nt_ticks)
            plt.show()
            chi_fig.savefig(f'{the_location}Tmin_Chisqr_{the_irrep[:4]}_{the_irrep[-1]}_{the_nr_exps}exp{the_rebin}_{the_version}.pdf', bbox_inches='tight')
            
        if the_total_chi_flag:
            print('Total Chi2 vs tmin plot in progress...')
            total_chi = vfa.TOTAL_CHI(the_chi_corr, the_nt, dis_set[1], 2)
            total_chi_fig = plt.figure() 
            vfp.PLOT_CHI_FITS(the_nt, total_chi, the_chosen_tmin, the_label, r'${t_{\mathrm{min}}}$', r'$\chi^{2}$', fr'{OperatorNamePlot} ({the_nr_exps}-exp)', the_nt_ticks)
            total_chi_fig.savefig(f'{the_location}Tmin_TotalChisqr_{the_irrep[:4]}_{the_irrep[-1]}_{the_nr_exps}exp{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        
        if the_delta_chi_flag:
            delta_chis = vfa.DELTA_CHI(total_chi)
            delta_chi_fig = plt.figure()
            vfp.PLOT_CHI_FITS(the_nt[:-1], delta_chis, the_chosen_tmin, r'$\Delta\chi^{2} = \chi^{2}({t}) - \chi^{2}({t}+1)$', r'${t_{\mathrm{min}}}$', r'$\Delta\chi^{2}$', fr'{OperatorNamePlot} ({the_nr_exps}-exp)', the_nt_ticks)
            
            delta_chi_fig.savefig(f'{the_location}Tmin_DeltaChisqr_{the_irrep[:4]}_{the_irrep[-1]}_{the_nr_exps}exp{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        
        if the_all_fits_flag:
            the_fits_avail = [kk for kk in the_single_fit_data[the_irrep].keys() if 'exp' in kk]
            print(f"All the fits available are: {the_fits_avail}")
            
            ### The SH operator that appears in the plot
            the_op = list(the_single_fit_data[f'{the_irrep}/Operators'])[0]
            theOperatorNamePlot = vfp.OPERATORS_SH(the_op.decode('utf-8'))
                    
            da_irrep = vfp.IrrepInfo(the_irrep)
            MomentumIrrep = da_irrep.TotalMomPlot
            nrMomentumIrrep = da_irrep.Momentum
            NameIrrepPlot = da_irrep.NamePlot 
        
            OperatorNamePlot = vfp.SH_OPERATORS_RELABEL(theOperatorNamePlot, NameIrrepPlot, nrMomentumIrrep )
            
            a_counter = 0
            for the_fit in the_fits_avail:
                ### This is the info of the fits for this irrep
                dis_set =  np.asarray(the_single_fit_data[f'{the_irrep}/{the_fit}/Tmin/{the_type_fit}/Mean'])
                
                ### All important data
                the_fit_data, the_fit_sigmas, the_chi_corr = dis_set[2], dis_set[3], dis_set[4]
                
                print('All fits vs tmin plot in progress...')
                if '1' in the_fit: 
                    the_label = "Single exponential"
                    the_color_1 = vfp.the_colors[0]
                    the_marker_1 = vfp.the_markers_list[0]
                elif '2' in the_fit: 
                    the_label = "Double exponential"
                    the_only_data = [the_nt, the_fit_data, the_fit_sigmas]
                    the_color_1 = vfp.the_colors[1]
                    the_marker_1 = vfp.the_markers_list[1]
                elif 'g' in the_fit: the_label = "Geometric series"

                plt.errorbar(the_nt, the_fit_data, yerr= the_fit_sigmas, marker=vfp.the_markers_list[a_counter], ls='None', ms=5, markeredgewidth=1.5, lw=1.5, elinewidth=1.5, zorder=3, capsize=4.5, label = the_label, color = vfp.the_colors[a_counter])
                a_counter+=1
            
            plt.errorbar([the_only_data[0][the_chosen_tmin]], [the_only_data[1][the_chosen_tmin]], yerr= [the_only_data[2][the_chosen_tmin]], marker=the_marker_1, ls='None', ms=6, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=6, color = the_color_1, markerfacecolor='white')
            
            plt.xlabel(r'$t_{\mathrm{min}}\,/\,a$',fontsize=36)
            plt.ylabel(r'$a\, E_{\mathrm{lab}}$',fontsize=36)
            plt.title(OperatorNamePlot, fontsize=26)
            plt.xticks(the_nt_ticks,fontsize=18)
            plt.yticks(fontsize=18)
            plt.legend(fontsize=16, handletextpad=0.01)
            plt.tight_layout()
            plt.show()
            fig_all.savefig(f'{the_location}Different_Fits_{the_irrep[:4]}_{the_irrep[-1]}_{the_rebin}_{the_version}.pdf', bbox_inches='tight')
        
        
def PlotMultiHadronsFits(the_multi_hadrons_fit_data, the_quantum_number, the_type_fit, the_nr_exps, the_type_rs, the_tmins, the_t0, the_version, the_location, the_rebin, the_irreps, **kwargs):  
    
    m_irreps = list(the_multi_hadrons_fit_data.keys())

    the_zoom_flag = kwargs.get('zoom_fit')
    the_chi_plots_flag = kwargs.get('chi_plots')
    the_total_chi_flag = kwargs.get('total_chi')
    the_delta_chi_flag = kwargs.get('delta_chi')
    the_ops_analysis_flag = kwargs.get('ops_analysis')
    
   ### How many irreps do you want to study        
    the_nr_irreps = kwargs.get('nr_irreps')
    the_first = kwargs.get('first_irrep')
    the_last = kwargs.get('last_irrep')

    if the_nr_irreps is not None:
        the_first_irrep = 0
        the_last_irrep = int(the_nr_irreps)
    else:
        the_first_irrep = int(the_first) - 1 if the_first is not None else 0
        the_last_irrep = int(the_last) if the_last is not None else len(m_irreps)
    
    ### The name of the irreps
    m_irreps = m_irreps[the_first_irrep:the_last_irrep]

    ### Loop over the irreps of this file
    for the_irrep in m_irreps:
        # print("in here now")
        print("---------------------------------------------------------------------------")
        print(f"Irrep: {the_irrep}")
        
        ### Searching if the fit was done and the gevp plots must be included
        if f'{the_nr_exps}exp' in list(the_multi_hadrons_fit_data[the_irrep].keys()): 
            
            ### Retrieving the data
            the_data = the_multi_hadrons_fit_data[f'{the_irrep}/{the_nr_exps}exp/t0_{the_t0}/Tmin/{the_type_fit}/Mean']
            
            ### Loop over the eigenvalues of this irrep
            for bb in range(len(list(the_data.keys()))):
                
                ### bb-th Eigenvalue
                dis_set = np.asarray(the_data[f'lambda_{bb}'])
                
                ### The central values of the diagonalized correlator
                the_fit_data = dis_set[2]
                
                ### The statistical errors
                the_fit_sigmas = dis_set[3]
                
                ### The chi^{2} of the fit
                the_chi_corr = dis_set[4]
                
                ### The statistical error of the chi^{2}
                the_chi_sigmas = dis_set[5]
                
                ## The minimum time slices that the fit was performed
                the_nt = [int(x) for x in dis_set[0]]
                
                ### The max time slice that used for the fit
                the_nt_max = int(dis_set[1][0])
                
                ### These are the ticks that appear in the plot
                the_nt_ticks = np.arange(the_nt[0], the_nt[-1], 2)

                ### Information about the irre
                da_irrep = vfp.IrrepInfo(the_irrep)
                MomentumIrrep = da_irrep.TotalMomPlot
                NameIrrepPlot = da_irrep.NamePlot

                ### Getting the position of the chosen tmin for the fit
                the_chosen_tmin = the_nt.index(the_tmins[the_irreps.index(the_irrep)][bb])
                
                the_mean_fit_string = str(np.round(the_fit_data[the_chosen_tmin], 5))
                the_error_string = vfp.WRITTING_ERRORS_PLOTS(the_fit_sigmas[the_chosen_tmin],5)
                the_sigmas_fit_string = the_error_string[0]
                    
                if the_error_string[1]==False:
                    the_mean_fit_string = str(f'{np.round(the_fit_data[the_chosen_tmin], the_error_string[2]):.{the_error_string[2]}f}')
                
                the_label = r'$[t_{\mathrm{min}}, t_{\mathrm{max}}] = [%sa, $'%str(int(the_nt[the_chosen_tmin])) + r'$%sa] $'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/\mathrm{d.o.f} = %s$'%np.round(the_chi_corr[the_chosen_tmin],3) + '\n' + r'$E_{\mathrm{fit}} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
                
                print('Fit vs Tmin plot in progress...')
                fit_fig = plt.figure()                
                vfp.PLOT_FITS(the_nt, the_fit_data, the_fit_sigmas, the_chosen_tmin, the_label, r'$t_{\mathrm{min}}\,/\,a$', r'$a \; E_{\mathrm{lab}}$',  f'{NameIrrepPlot} ({MomentumIrrep}): ' + r' $\to \;\lambda_{%s}$'%str(bb)  + r' ($t_{0} = %s$)'%str(the_t0), the_nt_ticks)
                plt.show()
                fit_fig.savefig(f'{the_location}Tmin_Fits_{the_quantum_number}_{the_irrep}_{bb}_t0_{the_t0}_{the_nr_exps}exp{the_rebin}_{the_version}.pdf', bbox_inches='tight')
                
                if the_zoom_flag:
                    print('Zoom Fit vs Tmin plot in progress...')
                    fit_fig = plt.figure()
                    vf.PLOT_FITS(the_nt, the_fit_data, the_fit_sigmas, the_chosen_tmin, the_label, r'$t_{\mathrm{min}}$', r'$a_{t} \;\Delta E$',  f'{NameIrrepPlot} ({MomentumIrrep}): ' + r' $\to \;\lambda_{%s}$'%str(bb)  + r' ($t_{0} = %s$)'%str(the_t0), the_nt_ticks, zoom=True, the_ll=2, the_ul =5)
                    plt.show()
                    fit_fig.savefig(f'{the_location}Tmin_Fits_Zoom_{the_quantum_number}_{the_irrep}_{bb}_t0_{the_t0}_{the_nr_exps}exp{the_rebin}_{the_version}.pdf', bbox_inches='tight')
                
                if the_chi_plots_flag:
                    print('Chi2 vs tmin plot in progress...')
                    chi_fit_fig = plt.figure()                
                    vf.PLOT_FITS(the_nt, the_chi_corr, the_chi_sigmas, the_chosen_tmin, the_label, r'$t_{\mathrm{min}}\,/\,a$', r'$\chi^{2}/\mathrm{d.o.f} $',  f'{NameIrrepPlot} ({MomentumIrrep}): ' + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%str(the_t0), the_nt_ticks)
                    plt.show()
                    chi_fit_fig.savefig(f'{the_location}Tmin_Chisqr_{the_quantum_number}_{the_irrep}_{bb}_t0_{the_t0}_{the_nr_exps}exp{the_rebin}_{the_version}.pdf', bbox_inches='tight')    
                
                if the_total_chi_flag:
                    print('Total Chi2 vs tmin plot in progress...')
                    total_chi = vfa.TOTAL_CHI(the_chi_corr, the_nt, dis_set[1], 2)
                    total_chi_fig = plt.figure()
                    plt.plot(the_nt, total_chi, marker='o', ls='None', ms=4, markeredgewidth=1.75, lw=1.75, zorder=3)
                    plt.plot([the_nt[the_chosen_tmin]], [total_chi[the_chosen_tmin]], marker='o', ls='None', ms=4, markeredgewidth=1.75, markerfacecolor='white', lw=1.75, zorder=3, label= the_label)
                    plt.xlabel(r'$t_{\mathrm{min}}$')
                    plt.ylabel(r'$\chi^{2}$')
                    plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %s$)'%str(the_t0))
                    plt.xticks(the_nt_ticks)
                    plt.legend()
                    plt.tight_layout()
                    plt.show()
                    total_chi_fig.savefig(f'{the_location}Tmin_TotalChisqr_{the_quantum_number}_{the_irrep}_{bb}_t0_{the_t0}_{the_nr_exps}exp{the_rebin}_{the_version}.pdf', bbox_inches='tight')                
                
                if the_delta_chi_flag:
                    print('Delta Chi2 vs tmin plot in progress...')
                    
                    delta_chis = vf.DELTA_CHI(total_chi)
                    
                    delta_chi_fig = plt.figure()
                    plt.plot(the_nt[:-1], delta_chis, marker='o', ls='None', ms=4, markeredgewidth=1.75, lw=1.75, zorder=3, label=r'$\Delta\chi^{2} = \chi^{2}(t) - \chi^{2}(t+1)$')
                    plt.plot([the_nt[the_chosen_tmin]], [delta_chis[the_chosen_tmin]], marker='o', ls='None', ms=4, markeredgewidth=1.75, markerfacecolor='white', lw=1.75, zorder=3, label=the_label)
                    plt.xlabel(r'$t_{\mathrm{min}}$')
                    plt.ylabel(r'$\Delta \chi^{2} $')
                    plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb) + r' ($t_{0} = %s$)'%str(the_t0))
                    plt.legend()
                    plt.xticks(the_nt_ticks)
                    plt.tight_layout()
                    plt.show()
                    delta_chi_fig.savefig(f'{the_location}Tmin_DeltaChisqr_{the_quantum_number}_{the_irrep}_{bb}_t0_{the_t0}_{the_nr_exps}exp{the_rebin}_{the_version}.pdf', bbox_inches='tight')
                
                
                the_t0_data = list(the_multi_hadrons_fit_data[f'{the_irrep}/{the_nr_exps}exp/'].keys())
                
                the_t0s = sorted([int(item[3:]) for item in the_t0_data])
                
                if len(the_t0_data)>2:
                    print('t0s vs energy plot in progress...')
                    the_t0s_axis = []
                    the_t0s_axis_error = []
                    for this_t0 in the_t0s:
                        the_t0s_axis.append(the_multi_hadrons_fit_data[f'{the_irrep}/{the_nr_exps}exp/t0_{this_t0}/Tmin/{the_type_fit}/Mean/lambda_{bb}'][2][the_chosen_tmin])
                        the_t0s_axis_error.append(the_multi_hadrons_fit_data[f'{the_irrep}/{the_nr_exps}exp/t0_{this_t0}/Tmin/{the_type_fit}/Mean/lambda_{bb}'][3][the_chosen_tmin])
                    
                    t0s_fig = plt.figure()
                    plt.errorbar(the_t0s, the_t0s_axis, yerr = the_t0s_axis_error, marker='o', ls='None', ms=4, markeredgewidth=1.75, lw=1.75, elinewidth=1.75, zorder=3, capsize=2.85)
                    
                    plt.errorbar([the_t0s[the_t0-the_t0s[0]]], [the_t0s_axis[the_t0-the_t0s[0]]], yerr = [the_t0s_axis_error[the_t0-the_t0s[0]]], capsize=2.75, marker='o', ls='None', ms=4, markeredgewidth=1.75, markerfacecolor='white', lw=1.75, zorder=3, label=the_label)
                    plt.legend()
                    plt.xlabel(r'$t_{0}$')
                    plt.ylabel(r'$a_{t} \;\Delta E$')
                    plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb))
                    plt.tight_layout()
                    #plt.show()
                    t0s_fig.savefig(f'{the_location}T0s_{the_quantum_number}_{the_irrep}_{bb}_{the_nr_exps}exp{the_rebin}_{the_version}.pdf', bbox_inches='tight')
                    
                        

def PlotRatioMultiHadronsFits(the_multi_hadrons_fit_data, the_quantum_number, the_type_fit, the_nr_exps, the_type_rs, the_tmins, the_t0, the_version, the_location, the_rebin, the_irreps, **kwargs):  
    
    m_irreps = list(the_multi_hadrons_fit_data.keys())
    # print("in here now")

    the_chi_plots_flag = kwargs.get('chi_plots')
    
    ### If not all the irreps are wanted t be plotted
    if kwargs.get('nr_irreps')!=None:
        the_first_irrep = 0
        the_last_irrep = int(kwargs.get('nr_irreps'))
    else:
        the_first_irrep = 0
        the_last_irrep = len(m_irreps)
    ### This one checks for an irrep in particular
    if kwargs.get('first_irrep')!=None and kwargs.get('last_irrep')!=None:
        the_first_irrep = int(kwargs.get('first_irrep'))
        the_last_irrep = int(kwargs.get('last_irrep'))
    elif kwargs.get('last_irrep')!=None and kwargs.get('first_irrep')==None:
        the_first_irrep = 0
        the_last_irrep = int(kwargs.get('last_irrep'))
    elif kwargs.get('first_irrep')!=None and kwargs.get('last_irrep')==None:
        the_first_irrep = int(kwargs.get('first_irrep'))
        the_last_irrep = len(m_irreps)
    
    m_irreps = m_irreps[the_first_irrep:the_last_irrep]
    

    ### Loop over the irreps of this file
    for the_irrep in m_irreps:
        # print("in here now")
        print("---------------------------------------------------------------------------")
        print("Irrep: ", the_irrep)
        
        ### Searching if the fit was done and the gevp plots must be included
        if '%sexp'%the_nr_exps in list(the_multi_hadrons_fit_data[the_irrep].keys()) and kwargs.get('gevp')==True: 
            
            ### Retrieving the data
            the_data = the_multi_hadrons_fit_data[the_irrep + '/%sexp/'%the_nr_exps + 't0_%s/Tmin/'%the_t0 + '%s/Mean'%the_type_fit]
            
            ### Loop over the eigenvalues of this irrep
            for bb in range(len(list(the_data.keys()))):
                
                ### bb-th Eigenvalue
                dis_set = np.array(the_data.get('lambda_%s'%bb))
                
                ### The central values of the diagonalized correlator
                the_fit_data = dis_set[2]
                
                ### The statistical errors
                the_fit_sigmas = dis_set[3]
                
                ### The chi^{2} of the fit
                the_chi_corr = dis_set[4]
                
                ### The statistical error of the chi^{2}
                the_chi_sigmas = dis_set[5]
                
                ## The minimum time slices that the fit was performed
                the_nt = [int(x) for x in dis_set[0]]
                
                ### The max time slice that used for the fit
                the_nt_max = int(dis_set[1][0])
                
                ### These are the ticks that appear in the plot
                the_nt_ticks = np.arange(the_nt[0], the_nt[-1], 2)

                ### Information about the irre
                da_irrep = vf.IrrepInfo(the_irrep)
                MomentumIrrep = da_irrep.TotalMomPlot
                NameIrrepPlot = da_irrep.NamePlot

                ### Getting the position of the chosen tmin for the fit
                the_chosen_tmin = the_nt.index(the_tmins[the_irreps.index(the_irrep)][bb])
                
                ### THIS IS OLD
                # the_chosen_tmin = the_tmins[the_irreps.index(the_irrep)][bb]-the_nt[0]
                
                the_mean_fit_string = str(np.round(the_fit_data[the_chosen_tmin], 5))
                the_error_string = vf.WRITTING_ERRORS_PLOTS(the_fit_sigmas[the_chosen_tmin],5)
                the_sigmas_fit_string = the_error_string[0]
                    
                if the_error_string[1]==False:
                    the_mean_fit_string = str(f'{np.round(the_fit_data[the_chosen_tmin], the_error_string[2]):.{the_error_string[2]}f}')
            
                # the_label = r'$t_{\mathrm{min}} = %s$'%str(int(the_nt[the_chosen_tmin])) + '\n' + r'$t_{\mathrm{max}}$ = %s'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/\mathrm{d.o.f} = %s$'%np.round(the_chi_corr[the_chosen_tmin],3) + '\n' + r'$E_{\mathrm{fit}} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
                
                the_label = r'$[t_{\mathrm{min}}, t_{\mathrm{max}}] = [%sa, $'%str(int(the_nt[the_chosen_tmin])) + r'$%sa] $'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/\mathrm{d.o.f} = %s$'%np.round(the_chi_corr[the_chosen_tmin],3) + '\n' + r'$E_{\mathrm{fit}} = %s$'%(the_mean_fit_string + the_sigmas_fit_string)
                
                print('Fit vs Tmin plot in progress...')
                fit_fig = plt.figure()                
                # vf.PLOT_FITS(the_nt, the_fit_data, the_fit_sigmas, the_chosen_tmin, the_label, r'$t_{\mathrm{min}}$', r'$a_{t} \;\Delta E$',  NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb)  + r' ($t_{0} = %s$)'%str(the_t0), the_nt_ticks)
                vf.PLOT_FITS(the_nt, the_fit_data, the_fit_sigmas, the_chosen_tmin, the_label, r'$t_{\mathrm{min}}\,/\,a$', r'$a \; E_{\mathrm{lab}}$',  NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%str(bb)  + r' ($t_{0} = %s$)'%str(the_t0), the_nt_ticks)
                plt.show()
                fit_fig.savefig(the_location + 'Tmin_Fits' + the_quantum_number + the_irrep + '_%s'%str(bb) + '_t0_%s'%str(the_t0)  + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version, bbox_inches='tight')
                
                if the_chi_plots_flag:
                    print('Chi2 vs tmin plot in progress...')
                    chi_fit_fig = plt.figure()                
                    vf.PLOT_FITS(the_nt, the_chi_corr, the_chi_sigmas, the_chosen_tmin, the_label, r'$t_{\mathrm{min}}\,/\,a$', r'$\chi^{2}/\mathrm{d.o.f} $',  NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%str(the_t0), the_nt_ticks)
                    plt.show()
                    chi_fit_fig.savefig(the_location + 'Tmin_Chisqr' + the_quantum_number + the_irrep + '_%s'%str(bb) + '_t0_%s'%str(the_t0) + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version, bbox_inches='tight')
            
                    
        
        
                



### ------------------------------- END FUNCTIONS ----------------------------------------------------



### --------------------------------------------------------------------------------------------------




### ------------------------------- START EXECUTING --------------------------------------------------




if __name__=="__main__":
    
    myEns = str(sys.argv[1]).upper()
    myWhichCorrelator = str(sys.argv[2]).lower()
    myTypeRs = str(sys.argv[3]).lower()
    myRebinOn = str(sys.argv[4]).lower()
    
    myTypeFit = 'Correlated'
    myNrExponentials =  '1'
    myRb = 2
    myVersion = 'test'
    myT0 = 4 
    
    myDataLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/CorrelatorData/%s/'%myEns)
    
    if myEns == 'N451': from files_n451 import singleTMinsFitPlots, multiTMinsFitPlots
    elif myEns == 'N201': from files_n201 import singleTMinsFitPlots, multiTMinsFitPlots 
    elif myEns == 'D200': from files_d200 import singleTMinsFitPlots, multiTMinsFitPlots
    elif myEns == 'X451': from files_x451 import singleTMinsFitPlots, multiTMinsFitPlots
    
    
    vf.INFO_PRINTING(myWhichCorrelator, myEns)
    
    if myRebinOn=='rb': 
        reBin = '_bin'+str(myRb)
    else:
        reBin = '' 
        
    if myTypeRs=='jk':
        myResamplingScheme='Jackknife'
    elif myTypeRs=='bt':
        myResamplingScheme='Bootstrap'    

    if myWhichCorrelator=='s':
        mySingleCorrelatorData = h5py.File(myDataLocation + 'Fits_SingleHadrons/Single_correlators_' + myTypeRs + reBin + '_fits_v%s.h5'%myVersion,'r')
        
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/Plots/%s/SingleHadrons/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotSingleHadronsFits(mySingleCorrelatorData, myTypeFit, myNrExponentials,  singleTMinsFitPlots, myVersion, myPlotLocation, reBin)
        
        mySingleCorrelatorData.close()
    
    elif myWhichCorrelator=='m':
        myMatrixCorrelatorData = h5py.File(myDataLocation + 'Fits_Matrices/Matrix_correlators_' + myTypeRs + reBin +'_fits_v%s.h5'%myVersion,'r')
        
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/Plots/%s/Matrices/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotMultiHadronsFits(myMatrixCorrelatorData, myTypeFit, myNrExponentials, multiTMinsFitPlots, myT0, myVersion, myPlotLocation, reBin)
        
        myMatrixCorrelatorData.close()
        
    elif myWhichCorrelator=='mr':
        myRatioCorrelatorData = h5py.File(myDataLocation + 'Fits_Ratios/Matrix_correlators_ratios_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
        
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/Documents/Chris Files/Plots/%s/Matrices_Ratios/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotMultiHadronsFits(myRatioCorrelatorData, myTypeFit, myNrExponentials, multiTMinsFitPlots, myT0, myVersion, myPlotLocation, reBin)
        
        myRatioCorrelatorData.close()
         
