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
    
def plot_many(rocs,aucs,labels):    
    plt.rcParams['figure.figsize'] = (8,8)
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.autolayout'] = True
    
    for i, roc in enumerate(rocs):
        fp, tp, threshs = roc
        plt.plot(tp, 1-fp, '-', label=labels[i]+' ({0:.3f})'.format(aucs[i]))
    
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
    plt.legend(loc='center left', frameon=False,fontsize=16)
    plt.text(0.02,0.02,r"$C_F$ = $C_A$ = 3, Kernel Order = $-1$",fontsize=18)
    plt.show()