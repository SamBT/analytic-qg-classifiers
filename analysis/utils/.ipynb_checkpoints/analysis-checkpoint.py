import numpy as np
import train_tools as tools
import plotting
import os
from sklearn.metrics import roc_auc_score, roc_curve

def get_jet_data(s_quark,s_gluon,prefix='/home/sambt/pythia83-samples/optimal-classifiers/kernel-1/',
                suffix='/total.root',nmax=-1,njets=-1,eta=-1):
    fq = prefix+s_quark+suffix
    fg = prefix+s_gluon+suffix
    X,Y = tools.pre_process(fq,fg,nev_max=nmax,njets=njets,eta=eta,quiet=False)
    z = X[:,:,0].copy()
    dEta = X[:,:,1].copy()
    dPhi = X[:,:,2].copy()
    dR = np.sqrt(dEta**2 + dPhi**2)
    mults = np.count_nonzero(z,axis=1)

    return z,dEta,dPhi,dR,mults,Y

def run_pipeline(s_quark,s_gluon,opt_func,prefix='/home/sambt/pythia83-samples/optimal-classifiers/kernel-1/',
                suffix='/total.root',nmax=-1,njets=-1,eta=-1,verbose=0,summary=False,quiet=False,predicted="Predicted",
                F_dropouts=0.2,epochs=3,modpath='/home/sambt/analytic-qg-classifiers/analysis/final-studies/saved_models/',savename=None):
    fq = prefix+s_quark+suffix
    fg = prefix+s_gluon+suffix
    X,Y = tools.pre_process(fq,fg,nev_max=nmax,njets=njets,eta=eta,quiet=False)
    z = X[:,:,0].copy()
    mults = np.count_nonzero(z,axis=1)
    optimal = opt_func(z)
    
    pfn, pfn_roc, pfn_auc = tools.train_qg_pfn(X,Y,n_epoch=epochs,summary=summary,verbose=verbose,F_dropouts=F_dropouts)
    pfn_z, pfn_z_roc, pfn_z_auc = tools.train_qg_pfn_no_angular(X,Y,n_epoch=epochs,summary=summary,verbose=verbose,F_dropouts=F_dropouts)
    pfn_noz, pfn_noz_roc, pfn_noz_auc = tools.train_qg_pfn_only_angular(X,Y,n_epoch=epochs,summary=summary,verbose=verbose,F_dropouts=F_dropouts)
    efn, efn_roc, efn_auc = tools.train_qg_efn(X,Y,n_epoch=epochs,summary=summary,verbose=verbose,F_dropouts=F_dropouts)
    
    optimal_auc = roc_auc_score(Y[:,1],-optimal)
    if optimal_auc < 0.5:
        optimal_auc = roc_auc_score(Y[:,1],optimal)
        optimal_roc = roc_curve(Y[:,1],optimal)
    else:
        optimal_roc = roc_curve(Y[:,1],-optimal)
    
    mult_auc = roc_auc_score(Y[:,1],mults)
    if mult_auc < 0.5:
        mult_auc = roc_auc_score(Y[:,1],-mults)
        mult_roc = roc_curve(Y[:,1],-mults)
    else:
        mult_roc = roc_curve(Y[:,1],mults)
    
    rocs = [pfn_roc,pfn_z_roc,pfn_noz_roc,efn_roc,optimal_roc,mult_roc]
    aucs = [pfn_auc,pfn_z_auc,pfn_noz_auc,efn_auc,optimal_auc,mult_auc]
    labels = ["PFN",r"PFN[$z$]",r"PFN[$\eta$,$\phi$]","EFN",predicted,"Multiplicity"]
    
    del X,Y,z,mults,optimal,pfn,pfn_z,pfn_noz,efn

    return rocs,aucs,labels

def duffPS_run_pipeline(fq,fg,opt_func,prefix='/home/sambt/pythia83-samples/optimal-classifiers/kernel-1/',
                suffix='/total.root',nmax=-1,njets=-1,eta=-1,verbose=1,summary=False,quiet=False,predicted="Predicted",
                F_dropouts=0.2,epochs=3,modpath='/home/sambt/analytic-qg-classifiers/analysis/final-studies/saved_models/',savename=None):
    X,Y = tools.duffPS_pre_process(fq,fg,nmax=nmax)
    z = X[:,:,0].copy()
    mults = np.count_nonzero(z,axis=1)
    optimal = opt_func(z)
    
    epochs = 5
    
    pfn, pfn_roc, pfn_auc = tools.duffPS_train_pfn(X,Y,n_epoch=epochs,summary=summary,verbose=verbose,F_dropouts=F_dropouts)
    pfn_z, pfn_z_roc, pfn_z_auc = tools.train_qg_pfn_no_angular(X,Y,n_epoch=epochs,summary=summary,verbose=verbose,F_dropouts=F_dropouts)
    
    optimal_auc = roc_auc_score(Y[:,1],-optimal)
    if optimal_auc < 0.5:
        optimal_auc = roc_auc_score(Y[:,1],optimal)
        optimal_roc = roc_curve(Y[:,1],optimal)
    else:
        optimal_roc = roc_curve(Y[:,1],-optimal)
    
    mult_auc = roc_auc_score(Y[:,1],mults)
    if mult_auc < 0.5:
        mult_auc = roc_auc_score(Y[:,1],-mults)
        mult_roc = roc_curve(Y[:,1],-mults)
    else:
        mult_roc = roc_curve(Y[:,1],mults)
    
    rocs = [pfn_roc,pfn_z_roc,optimal_roc,mult_roc]
    aucs = [pfn_auc,pfn_z_auc,optimal_auc,mult_auc]
    labels = [r"PFN[$z,\theta$]",r"PFN[$z$]",predicted,"Multiplicity"]
    
    del X,Y,z,mults,optimal,pfn

    return rocs,aucs,labels

def run_pfn(s_quark,s_gluon,opt_func,prefix='/home/sambt/pythia83-samples/optimal-classifiers/kernel-1/',
                nmax=-1,njets=-1,eta=-1,verbose=0,summary=False,F_dropouts=0.2,suffix='/total.root'):
    
    fq = prefix+s_quark+suffix
    fg = prefix+s_gluon+suffix
    X,Y = tools.pre_process(fq,fg,nev_max=nmax,njets=njets,eta=eta)
    z = X[:,:,0].copy()
    mults = np.count_nonzero(z,axis=1)
    optimal = opt_func(z)
    
    pfn, pfn_roc, pfn_auc = tools.train_qg_pfn(X,Y,n_epoch=3,summary=summary,verbose=verbose,F_dropouts=F_dropouts)
    
    optimal_auc = roc_auc_score(Y[:,1],-optimal)
    if optimal_auc < 0.5:
        optimal_auc = roc_auc_score(Y[:,1],optimal)
        optimal_roc = roc_curve(Y[:,1],optimal)
    else:
        optimal_roc = roc_curve(Y[:,1],-optimal)
    mult_roc = roc_curve(Y[:,1],mults)
    mult_auc = roc_auc_score(Y[:,1],mults)
    rocs = [pfn_roc,optimal_roc,mult_roc]
    aucs = [pfn_auc,optimal_auc,mult_auc]
    labels = ["PFN","Predicted","Multiplicity"]

    del X,Y,z,mults,optimal,pfn
    
    return rocs,aucs,labels

def run_custom_obs(s_quark,s_gluon,func,prefix='/home/sambt/pythia83-samples/optimal-classifiers/kernel-1/',
                suffix='/total.root',nmax=-1,njets=-1,eta=-1,quiet=False):
    
    fq = prefix+s_quark+suffix
    fg = prefix+s_gluon+suffix
    X,Y = tools.pre_process(fq,fg,nev_max=nmax,njets=njets,eta=eta,quiet=quiet)
    z = X[:,:,0].copy()
    
    if isinstance(func,list):
        optimal_rocs,optimal_aucs = [],[]
        for f in func:
            optimal = f(z)
            optimal_auc = roc_auc_score(Y[:,1],-optimal)
            if optimal_auc < 0.5:
                optimal_auc = roc_auc_score(Y[:,1],optimal)
                optimal_roc = roc_curve(Y[:,1],optimal)
            else:
                optimal_roc = roc_curve(Y[:,1],-optimal)
            optimal_rocs.append(optimal_roc)
            optimal_aucs.append(optimal_auc)
        del X,Y,z
        return optimal_rocs,optimal_aucs
    else:
        optimal = func(z)
        optimal_auc = roc_auc_score(Y[:,1],-optimal)
        if optimal_auc < 0.5:
            optimal_auc = roc_auc_score(Y[:,1],optimal)
            optimal_roc = roc_curve(Y[:,1],optimal)
        else:
            optimal_roc = roc_curve(Y[:,1],-optimal)
        del X,Y,z
        return optimal_roc,optimal_auc
