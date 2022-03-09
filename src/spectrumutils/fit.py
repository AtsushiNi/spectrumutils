import matplotlib.pyplot as plt
import numpy as np
import peakutils
from lmfit.models import GaussianModel, ConstantModel
from lmfit.lineshapes import gaussian

def single_peak_fit(x, y, peak_x, width=0.1):
    index = np.where((peak_x-width/2 < x) & (x < peak_x+width/2))
    xdata = x[index]
    ydata = y[index]

    gaussian_model = GaussianModel(prefix='gaussian_')
    constant_model = ConstantModel(prefix='constant_')

    params = gaussian_model.guess(ydata, x=xdata) + constant_model.guess(ydata, x=xdata)
    model = gaussian_model + constant_model

    result = model.fit(ydata, params, x = xdata)

    plt.plot(xdata, result.data, 'x', label='data', color='black')
    plt.plot(xdata, result.best_fit, color='red', label='fitted')
    plt.axhline(result.best_values['constant_c'], color='orange')
    plt.axvline(peak_x, color='gray', label='x='+str(peak_x))
    plt.legend()

    return {
        'amplitude': result.best_values['gaussian_amplitude'],
        'center': result.best_values['gaussian_center'],
        'sigma': result.best_values['gaussian_sigma'],
        'base': result.best_values['constant_c']
    }

def double_peak_fit(x, y, peak_x, other_peak, width=0.15):
    center = (peak_x+other_peak)/2
    index = np.where((center-width/2 < x) & (x < center+width/2))
    xdata = x[index]
    ydata = y[index]

    gaussian1_model = GaussianModel(prefix='gauss1_')
    gaussian2_model = GaussianModel(prefix='gauss2_')
    constant_model = ConstantModel(prefix='constant_')
    model = gaussian1_model + gaussian2_model + constant_model

    params = model.make_params()
    params_value = {
        'gauss1_amplitude': 7e4,
        'gauss1_center': peak_x,
        'gauss1_sigma': 0.01175,
        'gauss2_amplitude': 1e4,
        'gauss2_center': other_peak,
        'gauss2_sigma': 0.010,
        'constant_c': 4300
    }
    for name in model.param_names:
      params[name].set(value = params_value[name])

    result = model.fit(ydata, params, x = xdata)

    gauss1 = gaussian(
        xdata,
        result.best_values['gauss1_amplitude'],
        result.best_values['gauss1_center'],
        result.best_values['gauss1_sigma']
    )
    gauss2 = gaussian(
        xdata,
        result.best_values['gauss2_amplitude'],
        result.best_values['gauss2_center'],
        result.best_values['gauss2_sigma']
    )

    plt.plot(xdata, result.data, 'x', label='data', color='black')
    plt.plot(xdata, gauss1 + result.best_values['constant_c'], color='red', label='fitted')
    plt.plot(xdata, gauss2 + result.best_values['constant_c'], color='orange', label='other peak')
    plt.fill_between(xdata, (gauss1 + result.best_values['constant_c']), result.best_values['constant_c'], hatch='///', facecolor='None', edgecolor='red')
    plt.axhline(result.best_values['constant_c'], color='orange')
    plt.axvline(peak_x, color='gray', label='x='+str(peak_x))
    plt.legend()

    return {
        'amplitude': result.best_values['gauss1_amplitude'],
        'center': result.best_values['gauss1_center'],
        'sigma': result.best_values['gauss1_sigma'],
        'base': result.best_values['constant_c']
    }

def triple_peak_fit(x, y, peak_x, other_peak, width=0.25):
    center = np.average(np.append(other_peak, peak_x))
    index = np.where((center-width/2 < x) & (x < center+width/2))
    xdata = x[index]
    ydata = y[index]

    gaussian1_model = GaussianModel(prefix='gauss1_')
    gaussian2_model = GaussianModel(prefix='gauss2_')
    gaussian3_model = GaussianModel(prefix='gauss3_')
    constant_model = ConstantModel(prefix='constant_')
    model = gaussian1_model + gaussian2_model + gaussian3_model + constant_model

    params = model.make_params()
    params_value = {
        'gauss1_amplitude': 7e4,
        'gauss1_center': peak_x,
        'gauss1_sigma': 0.01175,
        'gauss2_amplitude': 1e4,
        'gauss2_center': other_peak[0],
        'gauss2_sigma': 0.010,
        'gauss3_amplitude': 1e4,
        'gauss3_center': other_peak[1],
        'gauss3_sigma': 0.010,
        'constant_c': 4300
    }
    for name in model.param_names:
      params[name].set(value = params_value[name])

    result = model.fit(ydata, params, x = xdata)

    gauss1 = gaussian(
        xdata,
        result.best_values['gauss1_amplitude'],
        result.best_values['gauss1_center'],
        result.best_values['gauss1_sigma']
    )
    gauss2 = gaussian(
        xdata,
        result.best_values['gauss2_amplitude'],
        result.best_values['gauss2_center'],
        result.best_values['gauss2_sigma']
    )
    gauss3 = gaussian(
        xdata,
        result.best_values['gauss3_amplitude'],
        result.best_values['gauss3_center'],
        result.best_values['gauss3_sigma']
    )

    plt.plot(xdata, result.data, 'x', label='data', color='black')
    plt.plot(xdata, gauss1 + result.best_values['constant_c'], color='red', label='fitted')
    plt.plot(xdata, gauss2 + result.best_values['constant_c'], color='orange', label='other peak1')
    plt.plot(xdata, gauss3 + result.best_values['constant_c'], color='limegreen', label='other peak2')
    plt.fill_between(xdata, (gauss1 + result.best_values['constant_c']), result.best_values['constant_c'], hatch='///', facecolor='None', edgecolor='red')
    plt.axhline(result.best_values['constant_c'], color='orange')
    plt.axvline(peak_x, color='gray', label='x='+str(peak_x))
    plt.legend()

    return {
        'amplitude': result.best_values['gauss1_amplitude'],
        'center': result.best_values['gauss1_center'],
        'sigma': result.best_values['gauss1_sigma'],
        'base': result.best_values['constant_c']
    }

def quadruple_peak_fit():
    print('quadruple')
