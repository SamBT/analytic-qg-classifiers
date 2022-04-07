#Global variables, useful imports, etc.

import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import train_tools as tools
import plotting
import analysis
from sklearn.metrics import roc_auc_score, roc_curve

plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.autolayout'] = True
plt.rcParams['text.usetex'] = True