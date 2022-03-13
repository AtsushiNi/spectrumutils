import numpy as np
import matplotlib.pyplot as plt
import copy
import plotly.graph_objects as go

# スペクトルを表示する
def show(wavelength, spectrum):
    pg_data = [
        go.Scatter(x=wavelength, y=spectrum)
    ]
    fig = go.Figure(data=pg_data)
    fig.update_layout(template='plotly_white')
    fig.show()

# 複数本のスペクトルを一本にまとめる
def unify2(wavelengths, spectra):
    new_wavelengths = [a for a in copy.deepcopy(wavelengths)]
    new_spectra = [a for a in copy.deepcopy(spectra)]

    for n in np.arange(wavelengths.shape[0]-1):
        left_wavelength = new_wavelengths[n]
        right_wavelength = new_wavelengths[n+1]
        left_spectra = new_spectra[n]
        right_spectra = new_spectra[n+1]

        for i in np.arange(int(left_wavelength.size/2)):
            if(left_wavelength[-1] < right_wavelength[0]):
                break
            left_wavelength = np.delete(left_wavelength, -1)
            right_wavelength = np.delete(right_wavelength, 0)
            left_spectra = np.delete(left_spectra, -1)
            right_spectra = np.delete(right_spectra, 0)

        new_wavelengths[n] = left_wavelength
        new_wavelengths[n+1] = right_wavelength
        new_spectra[n] = left_spectra
        new_spectra[n+1] = right_spectra

    return np.concatenate(new_wavelengths), np.concatenate(new_spectra)
