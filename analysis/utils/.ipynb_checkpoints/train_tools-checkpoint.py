from __future__ import absolute_import, division, print_function
import uproot
import numpy as np
# energyflow imports
import energyflow as ef
from energyflow.archs import *
from energyflow.utils import data_split, remap_pids, to_categorical

from sklearn.metrics import roc_auc_score, roc_curve

def train_qg_pfn(fname_q,fname_g):
    f_q = uproot.open(fname_q)['EventTree']
    f_g = uproot.open(fname_g)['EventTree']
    
    #Always train on leading jet from q/g samples
    gjet_pt = f_g.array("plead_constit_pt")
    gjet_eta = f_g.array("plead_constit_eta")
    gjet_phi = f_g.array("plead_constit_phi")

    qjet_pt = f_q.array("plead_constit_pt")
    qjet_eta = f_q.array("plead_constit_eta")
    qjet_phi = f_q.array("plead_constit_phi")

    #remove events where there is no leading quark or gluon jet for some reason
    g_mask = np.count_nonzero(gjet_pt,axis=1) > 0
    gjet_pt = gjet_pt[g_mask]
    gjet_eta = gjet_eta[g_mask]
    gjet_phi = gjet_phi[g_mask]

    q_mask = np.count_nonzero(qjet_pt,axis=1) > 0
    qjet_pt = qjet_pt[q_mask]
    qjet_eta = qjet_eta[q_mask]
    qjet_phi = qjet_phi[q_mask]

    #get size of constituent pt/eta/phi arrays; set to 100 in MC generation (more than needed) and padded with zeros
    pad_size = np.size(qjet_pt,axis=1)
    nev_gg = np.size(gjet_pt,axis=0)
    nev_qq = np.size(qjet_pt,axis=0)

    quarks = np.array([[[qjet_pt[i,j],qjet_eta[i,j],qjet_phi[i,j]] for j in range(pad_size)] for i in range(nev_qq)])
    gluons = np.array([[[gjet_pt[i,j],gjet_eta[i,j],gjet_phi[i,j]] for j in range(pad_size)] for i in range(nev_gg)])

    #make vectors with truth labels, combine q & g samples, shuffle
    quark_labs = np.ones(np.size(quarks,axis=0))
    glu_labs = np.zeros(np.size(gluons,axis=0))

    X = np.concatenate((quarks,gluons))
    y = np.concatenate((quark_labs,glu_labs))

    shuf = np.arange(np.size(X,axis=0))
    np.random.shuffle(shuf)

    X = X[shuf]
    y = y[shuf]

    #network parameters
    train, test, val = 75000, 10000, 15000
    Phi_sizes, F_sizes = (100, 100, 128), (100, 100, 100)
    num_epoch = 5
    batch_size = 500

    #convert quark/gluon labels to categorical
    Y = to_categorical(y,num_classes=2)

    # preprocess by centering jets and normalizing pts
    for x in X:
        mask = x[:,0] > 0
        yphi_avg = np.average(x[mask,1:3], weights=x[mask,0], axis=0)
        x[mask,1:3] -= yphi_avg
        x[mask,0] /= x[:,0].sum()

    print('Finished preprocessing')

    # do train/val/test split 
    (X_train, X_val, X_test,
     Y_train, Y_val, Y_test) = data_split(X, Y, train=train, val=val, test=test)

    print('Done train/val/test split')

    print('Model summary:')

    # build architecture
    pfn = PFN(input_dim=X.shape[-1], Phi_sizes=Phi_sizes, F_sizes=F_sizes)

    # train model
    pfn.fit(X_train, Y_train,
              epochs=num_epoch,
              batch_size=batch_size,
              validation_data=(X_val, Y_val),
              verbose=1)

    # get predictions on test data
    preds = pfn.predict(X_test, batch_size=1000)

    #pfn.model.save("saved_models/kern{0}/PFN_CF{1:.2f}_CA{2:.2f}".format(kern,CF,CA))

    #Make ROC curve
    roc = roc_curve(Y_test[:,1], preds[:,1])
    auc = roc_auc_score(Y_test[:,1], preds[:,1])
    print()
    print('PFN AUC:', auc)
    print()

    #Make ROCs for mass, multiplicity
    masses = np.asarray([ef.ms_from_p4s(ef.p4s_from_ptyphims(x).sum(axis=0)) for x in X])
    mults = np.asarray([np.count_nonzero(x[:,0]) for x in X])
    mass_roc = roc_curve(Y[:,1], -masses)
    mass_auc = roc_auc_score(Y[:,1],-masses)
    mult_roc = roc_curve(Y[:,1], -mults)
    mult_auc = roc_auc_score(Y[:,1],-mults)
    
    return (pfn, (roc,auc), (mass_roc,mass_auc), (mult_roc,mult_auc))