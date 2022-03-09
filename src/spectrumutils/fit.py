import matplotlib.pyplot as plt
import numpy as np
import peakutils
from lmfit.models import GaussianModel, ConstantModel
from lmfit.lineshapes import gaussian

def single_peak_fit(x, y, peak_x, width=0.1):
    index = np.where((peak_x-width < x) & (x < peak_x+width))
    xdata = x[index]
    ydata = y[index]

    gaussian_model = GaussianModel(prefix='gaussian_')
    constant_model = ConstantModel(prefix='constant_')

    params = gaussian_model.guess(ydata, x=xdata) + constant_model.guess(ydata, x=xdata)
    model = gaussian_model + constant_model

    result = model.fit(ydata, params, x = xdata)

    plt.plot(xdata, result.data, 'x', label='data')
    plt.plot(xdata, result.best_fit, color='orange', label='fitted')
    plt.axhline(result.best_values['constant_c'], color='orange')
    plt.axvline(peak_x, color='gray', label='x='+str(peak_x))
    plt.legend()

    return {
        'amplitude': result.best_values['gaussian_amplitude'],
        'center': result.best_values['gaussian_center'],
        'sigma': result.best_values['gaussian_sigma'],
        'base': result.best_values['constant_c']
    }

def double_peak_fit():
    print('double')

def triple_peak_fit():
    print('triple')

def quadruple_peak_fit():
    print('quadruple')
