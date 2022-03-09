import numpy as np
from importlib_resources import files

def read():
    frranck_condon = files('src/spectrumutils/data/fulcher').joinpath('franck_condon_factor.txt').read_text()
    # franck_condon_X_to_d = np.loadtxt('data/franck_condon_factor.txt')
