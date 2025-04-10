import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
import sys
import set_of_functions as vf


def PlotSingleHadronsFits(the_single_fit_data, the_type_fit, the_nr_exps, the_tmins, the_version, the_location, the_rebin):    
    s_irreps = list(the_single_fit_data.keys())

    for irrep in range(len(s_irreps)):
        dis_set =  np.array(the_single_fit_data[s_irreps[irrep] + '/%sexp'%the_nr_exps + '/Tmin/%s/Mean'%the_type_fit])
        
        the_mean_corr = dis_set[2]
        the_sigmas_corr = dis_set[3]
        the_chi_corr = dis_set[4]
        the_chi_sigmas = dis_set[5]
        
        the_nt = [int(x) for x in dis_set[0]]
        the_nt_max = int(dis_set[1][0])
        the_nt_ticks = np.arange(the_nt[0], the_nt[-1], 2)
                
        da_irrep = vf.IrrepInfo(s_irreps[irrep])
        MomentumIrrep = da_irrep.TotalMomPlot
        NameIrrepPlot = da_irrep.NamePlot
        NameIrrep = da_irrep.Name
        
        print('Fits vs tmin plot in progress...')
        fit_fig = plt.figure()
        plt.errorbar(the_nt, the_mean_corr, yerr = the_sigmas_corr, marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5)
        
        plt.errorbar([the_nt[the_tmins[irrep]-the_nt[0]]], [the_mean_corr[the_tmins[irrep]-the_nt[0]]], yerr = [the_sigmas_corr[the_tmins[irrep]-the_nt[0]]], marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, markerfacecolor = 'white' , capsize=2.5, label = r'$t_{min} = %s$'%str(int(the_nt[the_tmins[irrep]-the_nt[0]])) + '\n' +  r'$t_{max} = %s$'%str(int(the_nt_max))  + '\n' + r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_tmins[irrep]-the_nt[0]],3) + '\n' + r'$E_{fit} = %s$'%np.round(the_mean_corr[the_tmins[irrep]-the_nt[0]], 5) + '\n' + r'$\sigma_{fit} = %s$'%np.round(the_sigmas_corr[the_tmins[irrep]-the_nt[0]], 5))
        plt.legend()
        plt.xlabel(r'$t_{min}$')
        plt.ylabel(r'$a_{t} \;\Delta E_{lab}$')
        plt.title( s_irreps[irrep][-1] + ' (%s): '%MomentumIrrep + r' $\to$ %s '%NameIrrepPlot + '(%s-exp)'%the_nr_exps)
        plt.xticks(the_nt_ticks)
        plt.tight_layout()
        #plt.show()
        fit_fig.savefig(the_location + 'Tmin_Fits_' + s_irreps[irrep][:4] +'_%s'%s_irreps[irrep][-1] + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version)
        
        
        print('Zoom Fits vs tmin plot in progress...')
        zoom_fit_fig = plt.figure()
        plt.errorbar(the_nt[5:int(len(the_nt)/2)], the_mean_corr[5:int(len(the_nt)/2)], yerr = the_sigmas_corr[5:int(len(the_nt)/2)], marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5)
        plt.errorbar([the_nt[the_tmins[irrep]- the_nt[0] ]], [the_mean_corr[the_tmins[irrep]- the_nt[0] ]], yerr = [the_sigmas_corr[the_tmins[irrep]] - the_nt[0] ], marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, markerfacecolor = 'white' , capsize=2.5, label = r'$t_{min} = %s$'%str(int(the_nt[the_tmins[irrep]-the_nt[0]])) + '\n' +  r'$t_{max} = %s$'%str(int(the_nt_max))  + '\n' + r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_tmins[irrep]-the_nt[0]],3) + '\n' + r'$E_{fit} = %s$'%np.round(the_mean_corr[the_tmins[irrep]-the_nt[0]], 5))
        plt.legend()
        plt.xlabel(r'$t_{min}$')
        plt.ylabel(r'$a_{t} \;\Delta E_{lab}$')
        plt.title( s_irreps[irrep][-1] + ' (%s): '%MomentumIrrep + r' $\to$ %s'%NameIrrepPlot + '(%s-exp)'%the_nr_exps)
        plt.xticks(the_nt_ticks[3:6])
        plt.tight_layout()
        #plt.show()
        zoom_fit_fig.savefig(the_location + 'Tmin_Fits_Zoom_' + s_irreps[irrep][:4] +'_%s'%s_irreps[irrep][-1] + '_%sexp'%the_nr_exps + the_rebin +'_v%s.pdf'%the_version)
        
        
        
        print('Chi2 vs tmin plot in progress...')
        chi_fig = plt.figure()
        the_chi_nt_ticks = np.arange(the_nt[int(len(the_nt)/3)], the_nt[int(len(the_nt)*.82)], 2)
        plt.errorbar(the_nt[int(len(the_nt)/3):int(len(the_nt)*.82)], the_chi_corr[int(len(the_nt)/3):int(len(the_nt)*.82)], yerr = the_chi_sigmas[int(len(the_nt)/3):int(len(the_nt)*.82)], marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, zorder=3, capsize=2.5)
        plt.errorbar([the_nt[the_tmins[irrep] - the_nt[0] ]], the_chi_corr[the_tmins[irrep] - the_nt[0] ], yerr = the_chi_sigmas[the_tmins[irrep] - the_nt[0] ], marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, markerfacecolor = 'white' , capsize=2.5, label = r'$t_{min} = %s$'%str(int(the_nt[the_tmins[irrep]-the_nt[0]])) + '\n' +  r'$t_{max} = %s$'%str(int(the_nt_max))  + '\n' + r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_tmins[irrep]-the_nt[0]],3) + '\n' + r'$E_{fit} = %s$'%np.round(the_mean_corr[the_tmins[irrep]-the_nt[0]], 5))
        plt.legend()
        plt.xlabel(r'$t_{min}$')
        plt.ylabel(r'$\chi^{2}/d.o.f $')
        plt.title( s_irreps[irrep][-1] + ' (%s): '%MomentumIrrep + r' $\to$ %s'%NameIrrepPlot + '(%s-exp)'%the_nr_exps)
        plt.xticks(the_chi_nt_ticks)
        plt.tight_layout()
        #plt.show()
        chi_fig.savefig(the_location + 'Tmin_Chisqr_' + s_irreps[irrep][:4] +'_%s'%s_irreps[irrep][-1] + '_%sexp'%the_nr_exps + the_rebin+ '_v%s.pdf'%the_version)
        
        
        print('Total Chi2 vs tmin plot in progress...')
        total_chi = vf.TOTAL_CHI(the_chi_corr, the_nt, dis_set[1], 2)
        total_chi_fig = plt.figure()
        plt.plot(the_nt[int(len(the_nt)/3):int(len(the_nt)*.82)], total_chi[int(len(the_nt)/3):int(len(the_nt)*.82)], marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, zorder=3)
        plt.plot([the_nt[the_tmins[irrep] - the_nt[0] ]], total_chi[the_tmins[irrep] - the_nt[0] ],  marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, zorder=3, markerfacecolor = 'white', label = r'$t_{min} = %s$'%str(int(the_nt[the_tmins[irrep]-the_nt[0]])) + '\n' +  r'$t_{max} = %s$'%str(int(the_nt_max))  + '\n' + r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_tmins[irrep]-the_nt[0]],3) + '\n' + r'$E_{fit} = %s$'%np.round(the_mean_corr[the_tmins[irrep]-the_nt[0]], 5))
        plt.legend()
        plt.xlabel(r'$t_{min}$')
        plt.ylabel(r'$\chi^{2}$')
        plt.title( s_irreps[irrep][-1] + ' (%s): '%MomentumIrrep + r' $\to$ %s'%NameIrrepPlot + '(%s-exp)'%the_nr_exps)
        plt.xticks(the_chi_nt_ticks)
        plt.tight_layout()
        total_chi_fig.savefig(the_location + 'Tmin_TotalChisqr_' + s_irreps[irrep][:4] +'_%s'%s_irreps[irrep][-1] + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version)
        
        
        delta_chis = vf.DELTA_CHI(total_chi)
        delta_chi_fig = plt.figure()
        plt.plot(the_nt[int(len(the_nt)/3):int(len(the_nt)*.82)], delta_chis[int(len(the_nt)/3):int(len(the_nt)*.82)], marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, zorder=3, label=r'$\Delta\chi^{2} = \chi^{2}(t) - \chi^{2}(t+1)$')
        plt.plot([the_nt[the_tmins[irrep] - the_nt[0] ]], [delta_chis[the_tmins[irrep] - the_nt[0] ]], marker='o', ls='None', ms=4, markeredgewidth=1.1, markerfacecolor='white', lw=0.85, zorder=3, label = r'$t_{min} = %s$'%str(int(the_nt[the_tmins[irrep]-the_nt[0]])) + '\n' +  r'$t_{max} = %s$'%str(int(the_nt_max))  + '\n' + r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_tmins[irrep]-the_nt[0]],3) + '\n' + r'$E_{fit} = %s$'%np.round(the_mean_corr[the_tmins[irrep]-the_nt[0]], 5))        
        plt.legend()
        plt.xlabel(r'$t_{min}$')
        plt.ylabel(r'$\Delta\chi^{2}$')
        plt.title( s_irreps[irrep][-1] + ' (%s): '%MomentumIrrep + r' $\to$ %s'%NameIrrepPlot + '(%s-exp)'%the_nr_exps)
        plt.xticks(the_chi_nt_ticks)
        plt.tight_layout()
        #plt.show()
        delta_chi_fig.savefig(the_location + 'Tmin_DeltaChisqr_' + s_irreps[irrep][:4] +'_%s'%s_irreps[irrep][-1] + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version)
        
        
        
def PlotMultiHadronsFits(the_multi_hadrons_fit_data, the_type_fit, the_nr_exps, the_type_rs, the_tmins, the_t0, the_version, the_location, the_rebin):  
    m_irreps = list(the_multi_hadrons_fit_data.keys())

    for irrep in range(len(m_irreps)):
        the_data = the_multi_hadrons_fit_data[m_irreps[irrep] + '/%sexp/'%the_nr_exps + 't0_%s/Tmin/'%the_t0 + '%s/Mean'%the_type_fit]
        for bb in range(len(list(the_data.keys()))):
            dis_set = np.array(the_data.get('lambda_%s'%bb))
            
            the_mean_corr = dis_set[2]
            the_sigmas_corr = dis_set[3]
            the_chi_corr = dis_set[4]
            the_chi_sigmas = dis_set[5]
            #for abc in range(len(the_sigmas_corr)):
                #if np.isnan(the_sigmas_corr[abc]): the_sigmas_corr[abc]=0.00001
            the_nt = [int(x) for x in dis_set[0]]
            the_nt_max = int(dis_set[1][0])
            the_nt_ticks = np.arange(the_nt[0], the_nt[-1], 2)

            da_irrep = vf.IrrepInfo(m_irreps[irrep])
            MomentumIrrep = da_irrep.TotalMomPlot
            NameIrrepPlot = da_irrep.NamePlot
            NameIrrep = da_irrep.Name
            
            print('Fit vs Tmin plot in progress...')
            fit_fig = plt.figure()
            plt.errorbar(the_nt, the_mean_corr, yerr = the_sigmas_corr, marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5)
            #plt.plot(the_nt, the_mean_corr, marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, zorder=3)
            
            plt.errorbar([the_nt[the_tmins[irrep][bb]-the_nt[0]]], [the_mean_corr[the_tmins[irrep][bb]-the_nt[0]]], yerr = [the_sigmas_corr[the_tmins[irrep][bb]-the_nt[0]]], capsize=2.5, marker='o', ls='None', ms=4, markeredgewidth=1.1, markerfacecolor='white', lw=0.85, zorder=3, label=r'$t_{min} = %s$'%str(int(the_nt[the_tmins[irrep][bb]-the_nt[0]])) + '\n' + r'$t_{max}$ = %s'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_tmins[irrep][bb]-the_nt[0]],2) + '\n' + r'$E_{fit} = %s$'%np.round(the_mean_corr[the_tmins[irrep][bb]-the_nt[0]], 5))
            plt.legend()
            plt.xlabel(r'$t_{min}$')
            plt.ylabel(r'$a_{t} \;\Delta E$')
            plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb  + r' ($t_{0} = %s$)'%the_t0)
            plt.xticks(the_nt_ticks[2:-1])
            plt.tight_layout()
            #plt.show()
            fit_fig.savefig(the_location + 'Tmin_Fits_' + m_irreps[irrep] + '_%s'%bb  + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version)
            
            
            print('Zoom Fit vs Tmin plot in progress...')
            fit_fig = plt.figure()
            plt.errorbar(the_nt[int(len(the_nt)/3):int(len(the_nt)*.7)], the_mean_corr[int(len(the_nt)/3):int(len(the_nt)*.7)], yerr = the_sigmas_corr[int(len(the_nt)/3):int(len(the_nt)*.7)], marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5)
            
            plt.errorbar([the_nt[the_tmins[irrep][bb]-the_nt[0]]], [the_mean_corr[the_tmins[irrep][bb]-the_nt[0]]], yerr = [the_sigmas_corr[the_tmins[irrep][bb]-the_nt[0]]], capsize=2.5, marker='o', ls='None', ms=4, markeredgewidth=1.1, markerfacecolor='white', lw=0.85, zorder=3, label=r'$t_{min} = %s$'%str(int(the_nt[the_tmins[irrep][bb]-the_nt[0]])) + '\n' + r'$t_{max}$ = %s'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_tmins[irrep][bb]-the_nt[0]],2) + '\n' + r'$E_{fit} = %s$'%np.round(the_mean_corr[the_tmins[irrep][bb]-the_nt[0]], 5))
            plt.legend()
            plt.xlabel(r'$t_{min}$')
            plt.ylabel(r'$a_{t} \;\Delta E$')
            plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb  + r' ($t_{0} = %s$)'%the_t0)
            plt.xticks(the_nt_ticks[2:-1])
            plt.tight_layout()
            #plt.show()
            fit_fig.savefig(the_location + 'Tmin_Fits_Zoom_' + m_irreps[irrep] + '_%s'%bb + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version)
            
            
            print('Chi2 vs tmin plot in progress...')
            zoom_fit_fig = plt.figure()
            the_nt_chi_ticks = np.arange(int(len(the_nt)/3),int(len(the_nt)*.82), 2)
            plt.plot(the_nt[int(len(the_nt)/3):int(len(the_nt)*.82)], the_chi_corr[int(len(the_nt)/3):int(len(the_nt)*.82)], marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, zorder=3)
            
            plt.plot(the_nt[int(the_tmins[irrep][bb])-the_nt[0]], the_chi_corr[the_tmins[irrep][bb]-the_nt[0]], marker='o', ls='None', ms=4, markeredgewidth=1.1, markerfacecolor='white', lw=0.85, zorder=3, label=r'$t_{min} = %s$'%str(int(the_nt[the_tmins[irrep][bb]-the_nt[0]])) + '\n' + r'$t_{max}$ = %s'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_tmins[irrep][bb]-the_nt[0]],2) + '\n' + r'$E_{fit} = %s$'%np.round(the_mean_corr[the_tmins[irrep][bb]-the_nt[0]], 5))
            plt.xlabel(r'$t_{min}$')
            plt.legend()
            plt.ylabel(r'$\chi^{2}/d.o.f $')
            plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0)
            plt.xticks(the_nt_chi_ticks)
            plt.tight_layout()
            zoom_fit_fig.savefig(the_location + 'Tmin_Chisqr_' + m_irreps[irrep] + '_%s'%bb+ '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version)
  
  
            print('Total Chi2 vs tmin plot in progress...')
            total_chi = vf.TOTAL_CHI(the_chi_corr, the_nt, dis_set[1], 2)
            total_chi_fig = plt.figure()
            plt.plot(the_nt[int(len(the_nt)/3):int(len(the_nt)*.82)], total_chi[int(len(the_nt)/3):int(len(the_nt)*.82)], marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, zorder=3)
            plt.plot([the_nt[int(the_tmins[irrep][bb])-the_nt[0]]], [total_chi[the_tmins[irrep][bb]-the_nt[0]]], marker='o', ls='None', ms=4, markeredgewidth=1.1, markerfacecolor='white', lw=0.85, zorder=3, label=r'$t_{min} = %s$'%str(int(the_nt[the_tmins[irrep][bb]-the_nt[0]])) + '\n' + r'$t_{max}$ = %s'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_tmins[irrep][bb]-the_nt[0]],2) + '\n' + r'$E_{fit} = %s$'%np.round(the_mean_corr[the_tmins[irrep][bb]-the_nt[0]], 5))
            plt.xlabel(r'$t_{min}$')
            plt.legend()
            plt.ylabel(r'$\chi^{2}$')
            plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0)
            plt.xticks(the_nt_chi_ticks)
            plt.tight_layout()
            total_chi_fig.savefig(the_location + 'Tmin_TotalChisqr_' + m_irreps[irrep] +'_%s'%bb + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version)
            
            
            if len(the_nt)>1:
                print('Delta Chi2 vs tmin plot in progress...')
                delta_chis = vf.DELTA_CHI(total_chi)
                
                delta_chi_fig = plt.figure()
                plt.plot(the_nt[int(len(the_nt)/3):int(len(the_nt)*.82)], delta_chis[int(len(the_nt)/3):int(len(the_nt)*.82)], marker='o', ls='None', ms=4, markeredgewidth=1.1, lw=0.85, zorder=3, label=r'$\Delta\chi^{2} = \chi^{2}(t) - \chi^{2}(t+1)$')

                plt.plot([the_nt[the_tmins[irrep][bb]-the_nt[0]]], [delta_chis[int(the_tmins[irrep][bb])-the_nt[0]]], marker='o', ls='None', ms=4, markeredgewidth=1.1, markerfacecolor='white', lw=0.85, zorder=3, label=r'$t_{min} = %s$'%str(int(the_nt[the_tmins[irrep][bb]-the_nt[0]])) + '\n' + r'$t_{max}$ = %s'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_tmins[irrep][bb]-the_nt[0]],2) + '\n' + r'$E_{fit} = %s$'%np.round(the_mean_corr[the_tmins[irrep][bb]-the_nt[0]], 5))
                
                plt.xlabel(r'$t_{min}$')
                plt.legend()
                plt.ylabel(r'$\Delta \chi^{2} $')
                plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb + r' ($t_{0} = %s$)'%the_t0)
                plt.xticks(the_nt_chi_ticks)
                plt.tight_layout()
                #plt.show()
                delta_chi_fig.savefig(the_location + 'Tmin_DeltaChisqr_' + m_irreps[irrep] + '_%s'%bb + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version)
            
            
            the_t0_data = list(the_multi_hadrons_fit_data[m_irreps[irrep] + '/%sexp/'%the_nr_exps].keys())
            the_t0s = sorted([int(item[3:]) for item in the_t0_data])

            if len(the_t0_data)>2:
                print('t0s vs energy plot in progress...')
                the_t0s_axis = []
                the_t0s_axis_error = []
                for this_t0 in the_t0s:
                    the_t0s_axis.append(the_multi_hadrons_fit_data[m_irreps[irrep] + '/%sexp/'%the_nr_exps + 't0_%s/Tmin/'%this_t0 + '%s/Mean'%the_type_fit+ '/lambda_%s'%bb][2][the_tmins[irrep][bb]-the_nt[0]])
                    the_t0s_axis_error.append(the_multi_hadrons_fit_data[m_irreps[irrep] + '/%sexp/'%the_nr_exps + 't0_%s/Tmin/'%this_t0 + '%s/Mean'%the_type_fit+ '/lambda_%s'%bb][3][the_tmins[irrep][bb]-the_nt[0]])
                
                t0s_fig = plt.figure()
                plt.errorbar(the_t0s, the_t0s_axis, yerr = the_t0s_axis_error, marker='o', ls='None', ms=2.5, markeredgewidth=1.1, lw=0.85, elinewidth=0.85, zorder=3, capsize=2.5)
                
                plt.errorbar([the_t0s[the_t0-the_t0s[0]]], [the_t0s_axis[the_t0-the_t0s[0]]], yerr = [the_t0s_axis_error[the_t0-the_t0s[0]]], capsize=2.5, marker='o', ls='None', ms=4, markeredgewidth=1.1, markerfacecolor='white', lw=0.85, zorder=3, label=r'$t_{min} = %s$'%str(int(the_nt[the_tmins[irrep][bb]-the_nt[0]])) + '\n' + r'$t_{max}$ = %s'%str(int(the_nt_max)) + '\n' + r'$\chi^{2}/d.o.f = %s$'%np.round(the_chi_corr[the_tmins[irrep][bb]-the_nt[0]],2) + '\n' + r'$E_{fit} = %s$'%np.round(the_mean_corr[the_tmins[irrep][bb]-the_nt[0]], 5))
                plt.legend()
                plt.xlabel(r'$t_{0}$')
                plt.ylabel(r'$a_{t} \;\Delta E$')
                plt.title( NameIrrepPlot+ ' (%s): '%MomentumIrrep + r' $\to \;\lambda_{%s}$'%bb)
                #plt.xticks(the_nt_ticks[2:-1])
                plt.tight_layout()
                ##plt.show()
                t0s_fig.savefig(the_location + 'T0s_' + m_irreps[irrep] + '_%s'%bb + '_%sexp'%the_nr_exps + the_rebin + '_v%s.pdf'%the_version)
            



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
    
    myDataLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH(SAME_THAN_CORRS_SCRIPT_OUTPUT)$/%s/'%myEns)
    
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
        
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/SingleHadrons/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotSingleHadronsFits(mySingleCorrelatorData, myTypeFit, myNrExponentials,  singleTMinsFitPlots, myVersion, myPlotLocation, reBin)
        
        mySingleCorrelatorData.close()
    
    elif myWhichCorrelator=='m':
        myMatrixCorrelatorData = h5py.File(myDataLocation + 'Fits_Matrices/Matrix_correlators_' + myTypeRs + reBin +'_fits_v%s.h5'%myVersion,'r')
        
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/Matrices/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotMultiHadronsFits(myMatrixCorrelatorData, myTypeFit, myNrExponentials, multiTMinsFitPlots, myT0, myVersion, myPlotLocation, reBin)
        
        myMatrixCorrelatorData.close()
        
    elif myWhichCorrelator=='mr':
        myRatioCorrelatorData = h5py.File(myDataLocation + 'Fits_Ratios/Matrix_correlators_ratios_' + myTypeRs + reBin +'_v%s.h5'%myVersion,'r')
        
        myPlotLocation = vf.DIRECTORY_EXISTS(os.path.expanduser('~')+'/$YOUR_OUTPUT_PATH_FOR_PLOTS$/Plots/%s/Matrices_Ratios/'%myEns +  '%s/'%myResamplingScheme)
        
        PlotMultiHadronsFits(myRatioCorrelatorData, myTypeFit, myNrExponentials, multiTMinsFitPlots, myT0, myVersion, myPlotLocation, reBin)
        
        myRatioCorrelatorData.close()
         
