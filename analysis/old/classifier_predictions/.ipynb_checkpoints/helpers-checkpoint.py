import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../utils/")
from train_tools import pre_process, train_qg_pfn
from sklearn.metrics import roc_auc_score, roc_curve

def plot_rocs(pfn_out,optimal_roc,mult_roc):
    pfn_roc = pfn_out[1][0]

    plt.rcParams['figure.figsize'] = (8,8)
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.autolayout'] = True

    pfn_fp, pfn_tp, threshs = pfn_roc
    opt_fp, opt_tp, opt_thresh = optimal_roc
    mult_fp, mult_tp, mult_thresh = mult_roc
    plt.plot(pfn_tp, 1-pfn_fp, '-', label='PFN ROC')
    plt.plot(opt_tp,1-opt_fp,"-",label="Optimal ROC")
    plt.plot(mult_tp,1-mult_fp,"-",label="Multiplicity")
    
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
    plt.legend(loc='lower left', frameon=False,fontsize=14)
    plt.text(0.45,0.03,r"$C_F$ = $C_A$ = 3, Kernel Order = $-1$",fontsize=16)
    plt.show()
    
def nk(z,pwr):
    out = np.array([[f**pwr if f > 0 else 0 for f in j] for j in z])
    return out