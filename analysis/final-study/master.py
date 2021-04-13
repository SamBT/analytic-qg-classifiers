#Global variables, useful imports, etc.

import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
sys.path.append("../utils/")
import train_tools as tools
import plotting
import analysis as an
from sklearn.metrics import roc_auc_score, roc_curve

prefix = '/home/sambt/pythia83-samples/optimal-classifiers/kernel-1/'

f_gluon = prefix+'H2gg-CF3.0CA3.0-sqg0100/total.root'
f_quark = prefix+'H2qq-CF3.0CA3.0-sqg0100/total.root'

quark_variants_exp = ['sqg0100-esq0011','sqg0100-esq1000','sqg0100-esq1011','sqg0100-esq0001']
quark_variants_1overz_exp = ['sq1000-esq0000','sq1000-esq0011','sq1000-esq1000','sq1000-esq1011']