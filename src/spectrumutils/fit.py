import matplotlib.pyplot as plt
import numpy as np
import peakutils
from lmfit.models import GaussianModel, ConstantModel
from lmfit.lineshapes import gaussian

# ピークの数を自動判別してフィッティング
def fit(x, y , range):
    index = np.where((range[0] < x) & (x < range[-1]))
    ydata = y[index]
    xdata = x[index]

    peak_index = peakutils.indexes(ydata, thres=0.1)

    peak_num = peak_index.size
    if(peak_num == 0):
        raise Error('could not find any peak')
    peak_x = xdata[peak_index]
    width = xdata[-1] - xdata[0]
    if(peak_num == 1):
        return single_fit(xdata, ydata, peak_x, width=width)
    elif(peak_num == 2):
        return double_fit(xdata, ydata, peak_x, width=width)
    elif(peak_num == 3):
        return triple_fit(xdata, ydata, peak_x, width=width)
    elif(peak_num == 4):
        return quadruple_fit(xdata, ydata, peak_x, width=width)
    else:
        raise Error('wrong peak number')

def single_peak_fit(x, y, peak_x, width=0.1):
    result = single_fit(x, y, peak_x, width=width, target=True)

    return result[0]

def double_peak_fit(x, y, peak_x, other_peak, width=0.15):
    result = double_fit(x, y, np.insert(other_peak, 0, peak_x), width=width, target=True)

    return result[0]

def triple_peak_fit(x, y, peak_x, other_peak, width=0.25):
    result = triple_fit(x, y, np.insert(other_peak, 0, peak_x), width=width, target=True)

    return result[0]

def quadruple_peak_fit(x, y, peak_x, other_peak, width=0.25):
    result = quadruple_fit(x, y, np.insert(other_peak, 0, peak_x), width=width, target=True)

    return result[0]

def single_fit(x, y, peak_x, width=0.1, target=False):
    index = np.where((peak_x-width/2 < x) & (x < peak_x+width/2))
    xdata = x[index]
    ydata = y[index]

    gaussian_model = GaussianModel(prefix='gaussian_')
    constant_model = ConstantModel(prefix='constant_')

    params = gaussian_model.guess(ydata, x=xdata) + constant_model.guess(ydata, x=xdata)
    model = gaussian_model + constant_model

    result = model.fit(ydata, params, x = xdata)

    gauss = gaussian(
        xdata,
        result.best_values['gaussian_amplitude'],
        result.best_values['gaussian_center'],
        result.best_values['gaussian_sigma']
    )

    plt.plot(xdata, result.data, 'x', label='data', color='black')
    plt.plot(xdata, result.best_fit, color='red', label='fitted')
    plt.axhline(result.best_values['constant_c'], color='orange')
    if(target):
        plt.fill_between(xdata, (gauss + result.best_values['constant_c']), result.best_values['constant_c'], hatch='///', facecolor='None', edgecolor='red')
        plt.axvline(peak_x, color='gray', label='x='+str(peak_x))
    plt.legend()

    return [{
        'amplitude': result.best_values['gaussian_amplitude'],
        'center': result.best_values['gaussian_center'],
        'sigma': result.best_values['gaussian_sigma'],
        'base': result.best_values['constant_c']
    }]

def double_fit(x, y, peak_x, width=0.15, target=False):
    center = np.average(peak_x)
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
        'gauss1_center': peak_x[0],
        'gauss1_sigma': 0.01175,
        'gauss2_amplitude': 1e4,
        'gauss2_center': peak_x[1],
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
    plt.plot(xdata, gauss1 + result.best_values['constant_c'], color='red', label='peak1')
    plt.plot(xdata, gauss2 + result.best_values['constant_c'], color='orange', label='peak2')
    plt.axhline(result.best_values['constant_c'], color='orange')
    if(target):
        plt.axvline(peak_x[0], color='gray', label='x='+str(peak_x[0]))
        plt.fill_between(xdata, (gauss1 + result.best_values['constant_c']), result.best_values['constant_c'], hatch='///', facecolor='None', edgecolor='red')
    plt.legend()

    return [{
        'amplitude': result.best_values['gauss1_amplitude'],
        'center': result.best_values['gauss1_center'],
        'sigma': result.best_values['gauss1_sigma'],
        'base': result.best_values['constant_c']
    },
    {
        'amplitude': result.best_values['gauss2_amplitude'],
        'center': result.best_values['gauss2_center'],
        'sigma': result.best_values['gauss2_sigma'],
        'base': result.best_values['constant_c']
    }]

def triple_fit(x, y, peak_x, width=0.25, target=False):
    center = np.average(peak_x)
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
        'gauss1_center': peak_x[0],
        'gauss1_sigma': 0.01175,
        'gauss2_amplitude': 1e4,
        'gauss2_center': peak_x[1],
        'gauss2_sigma': 0.010,
        'gauss3_amplitude': 1e4,
        'gauss3_center': peak_x[2],
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
    plt.plot(xdata, gauss1 + result.best_values['constant_c'], color='red', label='peak1')
    plt.plot(xdata, gauss2 + result.best_values['constant_c'], color='orange', label='peak2')
    plt.plot(xdata, gauss3 + result.best_values['constant_c'], color='limegreen', label='peak3')
    plt.axhline(result.best_values['constant_c'], color='orange')
    if(target):
        plt.axvline(peak_x[0], color='gray', label='x='+str(peak_x[0]))
        plt.fill_between(xdata, (gauss1 + result.best_values['constant_c']), result.best_values['constant_c'], hatch='///', facecolor='None', edgecolor='red')
    plt.legend()

    return [{
        'amplitude': result.best_values['gauss1_amplitude'],
        'center': result.best_values['gauss1_center'],
        'sigma': result.best_values['gauss1_sigma'],
        'base': result.best_values['constant_c']
    },
    {
        'amplitude': result.best_values['gauss2_amplitude'],
        'center': result.best_values['gauss2_center'],
        'sigma': result.best_values['gauss2_sigma'],
        'base': result.best_values['constant_c']
    },
    {
        'amplitude': result.best_values['gauss3_amplitude'],
        'center': result.best_values['gauss3_center'],
        'sigma': result.best_values['gauss3_sigma'],
        'base': result.best_values['constant_c']
    }]

def quadruple_fit(x, y, peak_x, width=0.25, target=False):
    center = np.average(peak_x)
    index = np.where((center-width/2 < x) & (x < center+width/2))
    xdata = x[index]
    ydata = y[index]

    gaussian1_model = GaussianModel(prefix='gauss1_')
    gaussian2_model = GaussianModel(prefix='gauss2_')
    gaussian3_model = GaussianModel(prefix='gauss3_')
    gaussian4_model = GaussianModel(prefix='gauss4_')
    constant_model = ConstantModel(prefix='constant_')
    model = gaussian1_model + gaussian2_model + gaussian3_model + gaussian4_model + constant_model

    params = model.make_params()
    params_value = {
        'gauss1_amplitude': 7e4,
        'gauss1_center': peak_x[0],
        'gauss1_sigma': 0.01175,
        'gauss2_amplitude': 1e4,
        'gauss2_center': peak_x[1],
        'gauss2_sigma': 0.010,
        'gauss3_amplitude': 1e4,
        'gauss3_center': peak_x[2],
        'gauss3_sigma': 0.010,
        'gauss4_amplitude': 1e4,
        'gauss4_center': peak_x[3],
        'gauss4_sigma': 0.010,
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
    gauss4 = gaussian(
        xdata,
        result.best_values['gauss4_amplitude'],
        result.best_values['gauss4_center'],
        result.best_values['gauss4_sigma']
    )

    plt.plot(xdata, result.data, 'x', label='data', color='black')
    plt.plot(xdata, gauss1 + result.best_values['constant_c'], color='red', label='peak1')
    plt.plot(xdata, gauss2 + result.best_values['constant_c'], color='orange', label='peak2')
    plt.plot(xdata, gauss3 + result.best_values['constant_c'], color='limegreen', label='peak3')
    plt.plot(xdata, gauss4 + result.best_values['constant_c'], color='yellow', label='peak4')
    plt.axhline(result.best_values['constant_c'], color='orange')
    if(target):
        plt.fill_between(xdata, (gauss1 + result.best_values['constant_c']), result.best_values['constant_c'], hatch='///', facecolor='None', edgecolor='red')
        plt.axvline(peak_x[0], color='gray', label='x='+str(peak_x[0]))
    plt.legend()

    return [{
        'amplitude': result.best_values['gauss1_amplitude'],
        'center': result.best_values['gauss1_center'],
        'sigma': result.best_values['gauss1_sigma'],
        'base': result.best_values['constant_c']
    },
    {
        'amplitude': result.best_values['gauss2_amplitude'],
        'center': result.best_values['gauss2_center'],
        'sigma': result.best_values['gauss2_sigma'],
        'base': result.best_values['constant_c']
    },
    {
        'amplitude': result.best_values['gauss3_amplitude'],
        'center': result.best_values['gauss3_center'],
        'sigma': result.best_values['gauss3_sigma'],
        'base': result.best_values['constant_c']
    },
    {
        'amplitude': result.best_values['gauss4_amplitude'],
        'center': result.best_values['gauss4_center'],
        'sigma': result.best_values['gauss4_sigma'],
        'base': result.best_values['constant_c']
    }]
