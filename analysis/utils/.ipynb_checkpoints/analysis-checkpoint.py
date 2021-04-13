import numpy as np
def nk(z,k):
    out = np.array([np.sum([f**k if f > 0 else 0 for f in j]) for j in z])
    return out