import numpy as np
from importlib_resources import files

def read():
    """
    franck_condon = files('src.spectrumutils.data.fulcher').joinpath('franck_condon_factor.txt').read_text()
    print(franck_condon)
    franck_condon = files('src.spectrumutils.data.fulcher').joinpath('franck_condon_factor_d-to-a.csv').read_text()
    print(franck_condon)
    """
    franck_condon = files('src.spectrumutils.data.fulcher').joinpath('franck_condon_factor.txt')
    with as_file(franck_condon) as f:
        r = np.loadtxt(f)
    print(r)

