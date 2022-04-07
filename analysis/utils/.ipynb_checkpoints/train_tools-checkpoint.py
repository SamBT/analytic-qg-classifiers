from __future__ import absolute_import, division, print_function
import uproot
import numpy as np
# energyflow imports
import energyflow as ef
from energyflow.archs import *
from energyflow.utils import data_split, remap_pids, to_categorical

from sklearn.metrics import roc_auc_score, roc_curve

#import time library for testing purposes
import time

def pre_process(fname_q,fname_g,nev_max=-1,njets=-1,eta=-1,quiet=False):
    t_start = time.time()
    
    f_q = uproot.open(fname_q)['EventTree']
    f_g = uproot.open(fname_g)['EventTree']
    
    if not quiet: 
        print("Loaded files at "+str(time.time() - t_start))
    
    #Always train on leading jet from q/g samples
    gjet_pt = f_g.array("lead_constit_pt")
    gjet_eta = f_g.array("lead_constit_eta")
    gjet_phi = f_g.array("lead_constit_phi")
    gjet_njets = f_g.array("nJets")
    gjet_lead_etas = f_g.array("jet_eta")
    gjet_lead_etas = np.array([et[0] if len(et) > 0 else 99999 for et in gjet_lead_etas])

    qjet_pt = f_q.array("lead_constit_pt")
    qjet_eta = f_q.array("lead_constit_eta")
    qjet_phi = f_q.array("lead_constit_phi")
    qjet_njets = f_q.array("nJets")
    qjet_lead_etas = f_q.array("jet_eta")
    qjet_lead_etas = np.array([et[0] if len(et) > 0 else 99999 for et in qjet_lead_etas])
    
    if not quiet: 
        print("Read in arrays at "+str(time.time()-t_start))

    #remove events where there is no leading quark or gluon jet for some reason
    g_mask = np.count_nonzero(gjet_pt,axis=1) > 0
    if njets > 0:
        g_mask = np.logical_and(g_mask,gjet_njets == njets)
    if eta > 0:
        g_mask = np.logical_and(g_mask,np.abs(gjet_lead_etas) < eta)
    gjet_pt = gjet_pt[g_mask]
    gjet_eta = gjet_eta[g_mask]
    gjet_phi = gjet_phi[g_mask]

    q_mask = np.count_nonzero(qjet_pt,axis=1) > 0
    if njets > 0:
        q_mask = np.logical_and(q_mask,qjet_njets == njets)
    if eta > 0:
        q_mask = np.logical_and(q_mask,np.abs(qjet_lead_etas) < eta)
    qjet_pt = qjet_pt[q_mask]
    qjet_eta = qjet_eta[q_mask]
    qjet_phi = qjet_phi[q_mask]
    
    if not quiet:
        print("Cleaned events at "+str(time.time()-t_start))
    
    if nev_max != -1:
        gjet_pt = gjet_pt[:nev_max]
        gjet_eta = gjet_eta[:nev_max]
        gjet_phi = gjet_phi[:nev_max]
        qjet_pt = qjet_pt[:nev_max]
        qjet_eta = qjet_eta[:nev_max]
        qjet_phi = qjet_phi[:nev_max]
        
    #get size of constituent pt/eta/phi arrays; set to 100 in MC generation (more than needed) and padded with zeros
    #q_max_mult = np.max(np.array([np.count_nonzero(k) for k in qjet_pt]))
    #g_max_mult = np.max(np.array([np.count_nonzero(k) for k in gjet_pt]))
    q_max_mult = np.max(np.count_nonzero(qjet_pt,axis=1))
    g_max_mult = np.max(np.count_nonzero(gjet_pt,axis=1))
    pad_size = np.max([q_max_mult,g_max_mult])
    nev_gg = np.size(gjet_pt,axis=0)
    nev_qq = np.size(qjet_pt,axis=0)

    quarks = np.dstack((qjet_pt[:,:pad_size],qjet_eta[:,:pad_size],qjet_phi[:,:pad_size]))
    gluons = np.dstack((gjet_pt[:,:pad_size],gjet_eta[:,:pad_size],gjet_phi[:,:pad_size]))

    #quarks = np.array([[[qjet_pt[i,j],qjet_eta[i,j],qjet_phi[i,j]] for j in range(pad_size)] for i in range(nev_qq)])
    #gluons = np.array([[[gjet_pt[i,j],gjet_eta[i,j],gjet_phi[i,j]] for j in range(pad_size)] for i in range(nev_gg)])
    
    
    if not quiet:
        print("Made quark/gluon input arrays at "+str(time.time()-t_start))

    #make vectors with truth labels, combine q & g samples, shuffle
    quark_labs = np.ones(np.size(quarks,axis=0))
    glu_labs = np.zeros(np.size(gluons,axis=0))

    X = np.concatenate((quarks,gluons))
    y = np.concatenate((quark_labs,glu_labs))
    
    nev = np.size(X,axis=0)

    #convert quark/gluon labels to categorical
    Y = to_categorical(y,num_classes=2)

    #pre-process by centering jets and normalizing pts, using masked arrays
    msk = ~(np.abs(X) > 0)
    Xm = np.ma.array(X,mask=msk)
    pts = Xm[:,:,0]
    yphi_avg = np.ma.average(Xm[:,:,1:3],axis=1,weights=np.ma.repeat(pts,2,axis=1).reshape(pts.shape[0],pts.shape[1],2))
    Xm[:,:,1:3] -= np.ma.repeat(yphi_avg,Xm[:,:,1:3].shape[1],axis=0).reshape(yphi_avg.shape[0],Xm[:,:,1:3].shape[1],yphi_avg.shape[1])
    ptsum = np.sum(pts,axis=1)
    Xm[:,:,0] /= np.ma.repeat(ptsum,Xm[:,:,0].shape[1]).reshape(ptsum.shape[0],Xm[:,:,0].shape[1])
    X = Xm.filled(fill_value=0.0)

    # preprocess by centering jets and normalizing pts
    """for x in X:
        mask = x[:,0] > 0
        yphi_avg = np.average(x[mask,1:3], weights=x[mask,0], axis=0)
        x[mask,1:3] -= yphi_avg
        x[mask,0] /= x[:,0].sum()"""
    
    print('Finished preprocessing at '+str(time.time()-t_start))
    
    return (X,Y)

def duffPS_parse_jets(file,pad=200,nmax=-1):
    jet_flav = []
    jet_z = []
    jet_theta = []
    jet_phi = []
    jet_death = []
    
    nev = 0
    with open(file) as f:
        flav = []
        z = []
        theta = []
        phi = []
        death = []
        for l in f:
            if nev == nmax:
                break
            if "BEGIN" in l:
                flav = []
                z = []
                theta = []
                phi = []
                death = []
            elif "END" in l:
                if len(flav) < pad:
                    while len(flav) < pad:
                        flav.append(-99)
                        z.append(0.0)
                        theta.append(0.0)
                        phi.append(0.0)
                        death.append(0.0)
                jet_flav.append(np.array(flav))
                jet_z.append(np.array(z))
                jet_theta.append(np.array(theta))
                jet_phi.append(np.array(phi))
                jet_death.append(np.array(death))
                nev += 1
            else:
                vals = l.split(" ")
                flav.append(float(vals[0]))
                z.append(float(vals[1]))
                theta.append(float(vals[2]))
                phi.append(float(vals[3]))
                death.append(float(vals[4]))
    return np.array(jet_flav), np.array(jet_z), np.array(jet_theta), np.array(jet_phi), np.array(jet_death)

def duffPS_pre_process(fq,fg,nmax=-1):
    qf, qz, qth, qph, qd = duffPS_parse_jets(fq,nmax=nmax)
    gf, gz, gth, gph, gd = duffPS_parse_jets(fg,nmax=nmax)
    
    print('qz shape is {0}'.format(qz.shape))
    
    quarks = np.dstack((qz,qth))
    gluons = np.dstack((gz,gth))
    
    quark_labs = np.ones(np.size(quarks,axis=0))
    glu_labs = np.zeros(np.size(gluons,axis=0))

    X = np.concatenate((quarks,gluons))
    y = np.concatenate((quark_labs,glu_labs))
    
    nev = np.size(X,axis=0)

    #convert quark/gluon labels to categorical
    Y = to_categorical(y,num_classes=2)
    
    return (X,Y)

def duffPS_train_pfn(X,Y,summary=True,n_epoch=5,verbose=1,filepath=None,F_dropouts=0.2):
    
    t_start = time.time()
    
    #network parameters
    train, test, val = 0.7, 0.15, 0.15
    Phi_sizes, F_sizes = (100, 100, 128), (100, 100, 100)
    num_epoch = n_epoch
    batch_size = 500
    
    # do train/val/test split 
    (X_train, X_val, X_test,
     Y_train, Y_val, Y_test) = data_split(X, Y, train=train, val=val, test=test, shuffle=True)

    print('Training PFN[z,theta]')

    # build architecture
    pfn = PFN(input_dim=X.shape[-1], Phi_sizes=Phi_sizes, F_sizes=F_sizes,F_dropouts=F_dropouts,summary=summary,filepath=filepath)

    # train model
    pfn.fit(X_train, Y_train,
              epochs=num_epoch,
              batch_size=batch_size,
              validation_data=(X_val, Y_val),
              verbose=verbose)
    
    # get predictions on test data
    preds = pfn.predict(X_test, batch_size=1000)

    #Make ROC curve
    roc = roc_curve(Y_test[:,1], preds[:,1])
    auc = roc_auc_score(Y_test[:,1], preds[:,1])
    print('PFN[eta,phi] AUC:', auc)
    
    return pfn,roc,auc

def train_qg_pfn(X,Y,summary=True,n_epoch=5,verbose=1,filepath=None,F_dropouts=0.2):
    
    t_start = time.time()
    
    #network parameters
    train, test, val = 0.7, 0.15, 0.15
    Phi_sizes, F_sizes = (100, 100, 128), (100, 100, 100)
    num_epoch = n_epoch
    batch_size = 500
    
    # do train/val/test split 
    (X_train, X_val, X_test,
     Y_train, Y_val, Y_test) = data_split(X, Y, train=train, val=val, test=test, shuffle=True)

    print('Training PFN')

    # build architecture
    pfn = PFN(input_dim=X.shape[-1], Phi_sizes=Phi_sizes, F_sizes=F_sizes,F_dropouts=F_dropouts,summary=summary,filepath=filepath)

    # train model
    pfn.fit(X_train, Y_train,
              epochs=num_epoch,
              batch_size=batch_size,
              validation_data=(X_val, Y_val),
              verbose=verbose)
    
    # get predictions on test data
    preds = pfn.predict(X_test, batch_size=1000)

    #Make ROC curve
    roc = roc_curve(Y_test[:,1], preds[:,1])
    auc = roc_auc_score(Y_test[:,1], preds[:,1])
    print('PFN AUC:', auc)

    #Make ROCs for mass, multiplicity
    #masses = np.asarray([ef.ms_from_p4s(ef.p4s_from_ptyphims(x).sum(axis=0)) for x in X])
    #mults = np.asarray([np.count_nonzero(x[:,0]) for x in X])
    #mass_roc = roc_curve(Y[:,1], -masses)
    #mass_auc = roc_auc_score(Y[:,1],-masses)
    #mult_roc = roc_curve(Y[:,1], -mults)
    #mult_auc = roc_auc_score(Y[:,1],-mults)
    
    return pfn,roc,auc

def train_qg_pfn_no_angular(X,Y,summary=True,n_epoch=5,verbose=1,filepath=None,F_dropouts=0.2):
    
    t_start = time.time()
    
    #isolate JUST the z_i information
    X = X[:,:,:1]
    np.reshape(X,(X.shape[0],X.shape[1],1))
    
    #network parameters
    train, test, val = 0.7, 0.15, 0.15
    Phi_sizes, F_sizes = (100, 100, 128), (100, 100, 100)
    num_epoch = n_epoch
    batch_size = 500
    
    # do train/val/test split 
    (X_train, X_val, X_test,
     Y_train, Y_val, Y_test) = data_split(X, Y, train=train, val=val, test=test, shuffle=True)

    print('Training PFN[z]')

    # build architecture
    pfn = PFN(input_dim=X.shape[-1], Phi_sizes=Phi_sizes, F_sizes=F_sizes,F_dropouts=F_dropouts,summary=summary,filepath=filepath)

    # train model
    pfn.fit(X_train, Y_train,
              epochs=num_epoch,
              batch_size=batch_size,
              validation_data=(X_val, Y_val),
              verbose=verbose)

    # get predictions on test data
    preds = pfn.predict(X_test, batch_size=1000)

    #Make ROC curve
    roc = roc_curve(Y_test[:,1], preds[:,1])
    auc = roc_auc_score(Y_test[:,1], preds[:,1])
    print('PFN[z] AUC:', auc)
    
    return pfn,roc,auc

def train_qg_pfn_only_angular(X,Y,summary=True,n_epoch=5,verbose=1,filepath=None,F_dropouts=0.2):
    
    t_start = time.time()
    
    #isolate JUST the angular information
    X = X[:,:,1:]
    
    #network parameters
    train, test, val = 0.7, 0.15, 0.15
    Phi_sizes, F_sizes = (100, 100, 128), (100, 100, 100)
    num_epoch = n_epoch
    batch_size = 500
    
    # do train/val/test split 
    (X_train, X_val, X_test,
     Y_train, Y_val, Y_test) = data_split(X, Y, train=train, val=val, test=test, shuffle=True)

    print('Training PFN[eta,phi]')

    # build architecture
    pfn = PFN(input_dim=X.shape[-1], Phi_sizes=Phi_sizes, F_sizes=F_sizes,F_dropouts=F_dropouts,summary=summary,filepath=filepath)

    # train model
    pfn.fit(X_train, Y_train,
              epochs=num_epoch,
              batch_size=batch_size,
              validation_data=(X_val, Y_val),
              verbose=verbose)
    
    # get predictions on test data
    preds = pfn.predict(X_test, batch_size=1000)

    #Make ROC curve
    roc = roc_curve(Y_test[:,1], preds[:,1])
    auc = roc_auc_score(Y_test[:,1], preds[:,1])
    print('PFN[eta,phi] AUC:', auc)
    
    return pfn,roc,auc

def train_qg_efn(X,Y,summary=True,n_epoch=5,verbose=1,filepath=None,F_dropouts=0.2):
    
    t_start = time.time()
    
    #network parameters
    train, test, val = 0.7, 0.15, 0.15
    Phi_sizes, F_sizes = (100, 100, 128), (100, 100, 100)
    num_epoch = n_epoch
    batch_size = 500
    
    # do train/val/test split 
    (z_train, z_val, z_test, 
     p_train, p_val, p_test,
     Y_train, Y_val, Y_test) = data_split(X[:,:,0], X[:,:,1:], Y, val=val, test=test)

    print('Training EFN')

    # build architecture
    efn = EFN(input_dim=2, Phi_sizes=Phi_sizes, F_sizes=F_sizes,F_dropouts=F_dropouts,summary=summary,filepath=filepath)

    # train model
    efn.fit([z_train, p_train], Y_train,
              epochs=num_epoch,
              batch_size=batch_size,
              validation_data=([z_val, p_val], Y_val),
              verbose=verbose)
    
    # get predictions on test data
    preds = efn.predict([z_test, p_test], batch_size=1000)

    #Make ROC curve
    roc = roc_curve(Y_test[:,1], preds[:,1])
    auc = roc_auc_score(Y_test[:,1], preds[:,1])
    print('EFN AUC:', auc)
    
    return efn,roc,auc

def mapto(arr,xmin,xmax):
    r = xmax - xmin
    d = arr.max() - arr.min()
    out = r*arr/d
    out = out - out.min() + xmin
    return out
