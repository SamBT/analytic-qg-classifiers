import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../utils/")
from train_tools import pre_process, train_qg_pfn
from sklearn.metrics import roc_auc_score, roc_curve

def plot_rocs(pfn_roc,optimal_roc,mult_roc,pfn_auc,optimal_auc,mult_auc):    
    plt.rcParams['figure.figsize'] = (8,8)
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.autolayout'] = True

    pfn_fp, pfn_tp, threshs = pfn_roc
    opt_fp, opt_tp, opt_thresh = optimal_roc
    mult_fp, mult_tp, mult_thresh = mult_roc
    plt.plot(pfn_tp, 1-pfn_fp, '-', label='PFN ({0:.3f})'.format(pfn_auc))
    plt.plot(opt_tp,1-opt_fp,"-",label="Predicted ({0:.3f})".format(optimal_auc))
    plt.plot(mult_tp,1-mult_fp,"-",label="Multiplicity ({0:.3f})".format(mult_auc))
    
    #Draw 'random' line
    plt.plot([0,1],[1,0],linestyle='dashed',color='gray')

    # axes labels
    plt.xlabel('Quark Jet Efficiency',fontsize=14)
    plt.ylabel('Gluon Jet Rejection',fontsize=14)

    # axes limits
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # make legend and show plot
    plt.legend(loc='lower left', frameon=False,fontsize=16)
    plt.text(0.02,0.2,r"$C_F$ = $C_A$ = 3, Kernel Order = $-1$",fontsize=18)
    plt.show()
    
def plot_many(rocs,aucs,labels,ax,ncol=1,fontsize=16):    
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.autolayout'] = True
    plt.rcParams['text.usetex'] = True
    
    for i, roc in enumerate(rocs):
        fp, tp, threshs = roc
        ax.plot(tp, 1-fp, '-', label=labels[i]+' ({0:.3f})'.format(aucs[i]))
    
    #Draw 'random' line
    ax.plot([0,1],[1,0],linestyle='dashed',color='gray')

    # axes labels
    plt.xlabel('Quark Jet Efficiency',fontsize=fontsize)
    plt.ylabel('Gluon Jet Rejection',fontsize=fontsize)

    # axes limits
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xticks(fontsize=fontsize-2)
    plt.yticks(fontsize=fontsize-2)

    # make legend and show plot
    plt.legend(loc='lower left', frameon=False,fontsize=fontsize,ncol=ncol)
    #plt.show()
    
def stamp(left_x, top_y,radius=0.4,fontsize=16,
          line_0=None,
          line_1=r'\textsc{Pythia} 8.303 + \textsc{Dire}',
          line_2=r'$e^+e^- \to H \to q\bar{q}/gg$',
          line_3=r'$\sqrt{s} = 250$ GeV, $C_F = C_A = 3$',
          line_4=r'Parton Level, R = 0.4',
          delta_y=0.035,
          textops_update={},
          ax=None,
          **kwargs):
    
    line_4 = r'Hadron Level, R = {0:.1f}'.format(radius)
    textops = {'horizontalalignment': 'left', 'verticalalignment': 'center', 'fontsize': fontsize}
    textops.update(textops_update)
    if ax is None: ax = plt.gca()
    textops['transform'] = ax.transAxes
    kwargs.update({'line_0':line_0,'line_1': line_1, 'line_2': line_2, 'line_3': line_3, 'line_4':line_4})
    for i in range(len(kwargs)):
        y = top_y - i*delta_y
        t = kwargs.get('line_' + str(i))
        if t is not None:
            ax.text(left_x, y, t, **textops)
            
def duffPS_stamp(left_x, top_y,radius=0.4,fontsize=16,
          line_0=None,
          line_1=r'Simple Parton Shower',
          line_2=r'LL Splittings, $C_F = C_A = 3$',
          delta_y=0.035,
          textops_update={},
          ax=None,
          **kwargs):
    
    textops = {'horizontalalignment': 'left', 'verticalalignment': 'center', 'fontsize': fontsize}
    textops.update(textops_update)
    if ax is None: ax = plt.gca()
    textops['transform'] = ax.transAxes
    kwargs.update({'line_0':line_0,'line_1': line_1, 'line_2': line_2})
    for i in range(len(kwargs)):
        y = top_y - i*delta_y
        t = kwargs.get('line_' + str(i))
        if t is not None:
            ax.text(left_x, y, t, **textops)