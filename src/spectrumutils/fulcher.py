import numpy as np
from importlib_resources import files, as_file
import xarray as xr
import matplotlib.pyplot as plt
import plotly.graph_objects as go

# スペクトルを表示する
# 一つまたは複数のスペクトルを表示
def show(wavelengths, spectra, labels=[""], v=[0]):
    if(len(spectra) > 40): # 一つのスペクトルを表示する場合
        wavelengths = [wavelengths]
        spectra = [spectra]
    max = np.max(np.array(spectra).flatten())

    pg_data = [
        go.Scatter(x=wavelength, y=spectrum, name=label)
        for wavelength, spectrum, label in zip(wavelengths, spectra, labels)
    ]
    fig = go.Figure(data=pg_data)

    wavelength_data = fulcher_wavelength_npy()
    for V in v:
        for vv,lines in enumerate(wavelength_data):
            if vv != V:
                continue
            for line in lines:
                if line != 0:
                    fig.add_shape(
                        type='line',
                        x0=line,y0=0,x1=line,y1=max*1.1,
                        line=dict(width=1, color='gray')
                    )
    fig.update_layout(template='plotly_white')
    fig.show()

# Q1を使ってざっくり波長校正
# borderにはQ1を検出する範囲を指定(v=0,1,2...)
def calibrate_by_Q1(pixel, spectrum, border):
    Q1_pixels = []
    Q1_spectrum = []
    for i in np.arange(len(border)-1):
        index = np.where((border[i] < pixel) & (pixel < border[i+1]))
        pixel_range = pixel[index]
        spectrum_range = spectrum[index]
        Q_index = np.argmax(spectrum_range)
        Q1_pixels.append(pixel_range[Q_index])
        Q1_spectrum.append(spectrum_range[Q_index])

    Q1_pixels = np.array(Q1_pixels)
    Q1_wavelengths = fulcher_wavelength().sel(dv=slice(0,Q1_pixels.size-1), dN=1).values

    fit = np.poly1d(np.polyfit(Q1_pixels, Q1_wavelengths, 1))
    plt.plot(fit(pixel), spectrum, color='black')
    plt.plot(fit(Q1_pixels), Q1_spectrum, 'oC1')

    return fit(pixel)

# 発光線を特定。波長のズレを得る
def detect_lines(wavelength, spectrum, lines, width=0.1):
    plt.figure(figsize=(16,8))
    plt.subplots_adjust(wspace=0.3, hspace=0.5)
    line_num = sum([np.array(a).size for a in lines])
    fulcher_wdata = fulcher_wavelength()
    i=0
    peak_wavelength = []
    target_wavelength = []
    for v, Ns in enumerate(lines):
        for N in Ns:
            target_w = fulcher_wdata.sel(dv=v, dN=N).values
            target_index = np.where((target_w-width<wavelength)&(wavelength<target_w+width))
            display_index = np.where((target_w-0.2<wavelength)&(wavelength<target_w+0.2))
            peak_index = np.argmax(spectrum[target_index])
            peak_w = wavelength[target_index][peak_index]
            peak_s = spectrum[target_index][peak_index]
            i += 1
            plt.subplot(line_num//4+1, 4, i)
            plt.plot(wavelength, spectrum, color='black')
            plt.axvline(target_w, color='red')
            plt.axvline(target_w-width, color='C0')
            plt.axvline(target_w+width, color='C0')
            plt.plot(peak_w, peak_s, 'xC1')
            plt.xlim(target_w-0.2, target_w+0.2)
            plt.ylim(0, np.max(spectrum[target_index])*1.1)
            plt.title('v='+str(v)+', N='+str(N))

            peak_wavelength.append(peak_w)
            target_wavelength.append(target_w)

    return np.array(peak_wavelength), np.array(target_wavelength)

# 各発光線を使って精確に波長校正
def calibrate(wavelength, spectrum, lines, width=0.1):
    return 0

# 発光強度データからボルツマンプロットを作成
def boltzmannplot(data, v, errors=None):
    # 振動準位ごとに
    result = []
    for i, (d, dv) in enumerate(zip(data, v)):
        index = np.nonzero(d)
        N_numbers = index[0]+1
        amplitudes = d[index]

        population = np.zeros(amplitudes.size)
        rot_energy = E_d_rot(dv, N_numbers)

        for j, N in enumerate(N_numbers):
            population[j] = amplitudes[j] * fulcher_wavelength().sel(dv=dv,dN=N) **4 / (2*N+1)/ g_as(N)

        if(errors is None):
            plt.plot(rot_energy, population, '--x')
        else:
            plt.errorbar(rot_energy, population, yerr=error[i])
        plt.yscale('log')
        plt.xlabel('Rotational Energy (eV)')
        plt.ylabel('population (a.u.)')

        result.append({
            'N_numbers': N_numbers,
            'rotation_energies': rot_energy,
            'population': population
        })

    return result

# fulcher-aの波長データ
def fulcher_wavelength():
    franck_condon = files('spectrumutils.data.fulcher').joinpath('fulcher_wavelength.nc')
    with as_file(franck_condon) as f:
        r = xr.open_dataarray(f)
    return r

def fulcher_wavelength_npy(dv=None, dN=None):
    data = files("spectrumutils.data.fulcher").joinpath("fulcher_wavelength.npy")
    with as_file(data) as f:
        r = np.load(f)

    if((dv is None) and (dN is None)):
        return r

    if(dN is None):
        return r[dv]
    if(dv is None):
        return r.T[dN-1]

    return r[dv][dN-1]

# Franck-Condin因子(https://inis.iaea.org/collection/NCLCollectionStore/_Public/37/088/37088524.pdf)
def franck_condon_X_to_d():
    franck_condon = files('spectrumutils.data.fulcher').joinpath('franck_condon_factor.txt')
    with as_file(franck_condon) as f:
        r = np.loadtxt(f)
    return r

def franck_condon_d_to_a():
    franck_condon = files('spectrumutils.data.fulcher').joinpath('franck_condon_factor_d-to-a.csv')
    with as_file(franck_condon) as f:
        r = np.loadtxt(f, delimiter=',')
    return r

# プランク定数[ev/s]
planck_constant = 4.13567e-15
# 光速[m/s]
light_speed = 299792458
# ボルツマン定数
kb = 8.6171e-5 # [eV K^-1]

# NISTより(https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=1000#Diatomic)
Be = [60.853, 30.364, 1.671] # X, d, a準位のBe
ae = [3.062, 1.545, 1.671] # X, d, a準位のαe
De = [0.0471, 0.0191, 0.0216] # X, d, a準位のDe
we = [4401.21, 2371.57, 2664.83] # X, d, a準位のwe
wx = [121.33, 66.27, 71.65] # X, d, a準位のwx
Te = 7 / kb # 電子温度[K] (5eV)

# 核スピンの縮退度(核スピンの統計重率)
def g_as(N):
  return (N % 2) * 2 + 1

###### 回転定数
# 基底準位の回転定数
def B_X(Xv):
  return Be[0] - ae[0]*(Xv+0.5)
# 上準位の回転定数
def B_d(dv):
  return Be[1] - ae[1]*(dv+0.5)

###### 回転エネルギー
# 基底準位の回転エネルギー(eV)
def E_X_rot(Xv, XN):
  return (B_X(Xv)*XN*(XN+1) - De[0]*XN*(XN+1)*XN*(XN+1))*1.23984/1e4
# 上状態の回転エネルギー(eV)
def E_d_rot(dv, dN):
  return (B_d(dv)*dN*(dN+1)-De[1]*dN*(dN+1)*dN*(dN+1))*1.23984/1e4

###### 振動エネルギー
# 基底準位の振動エネルギー
def E_X_vib(v):
  return (we[0]*(v+0.5) - wx[0]*(v+0.5)*(v+0.5))*1.23984/1e4

# 上準位の振動エネルギー
def E_d_vib(v):
  return (we[1]*(v+0.5) - wx[1]*(v+0.5)*(v+0.5))*1.23984/1e4

###### R導出用の諸々
# 振動励起の電子衝突断面積(の多分相対値)
def ccs(vX, vd):
  return np.exp(-((E_d_vib(vd)-E_d_vib(0))-(E_X_vib(vX) - E_X_vib(0)))/ kb / Te)
# 電子速度で平均化した部分断面積　石原さんの修論より
Qr = np.array([0.76,0.122,0.1,0.014])
# 回転構造の分岐比
def branch_ratio(NX, Nd):
  rtp_value = 0
  for r in range(1, 5, 1):
    rtp_value += Qr[r-1]*(2*Nd+1)*(wigner_3j(Nd, r, NX, 1, -1, 0))**2
  return rtp_value
# 禁止遷移は励起係数は0
def kronecker_delta(NX, Nd):
  if(NX % 2 == Nd % 2):
    return 1
  else:
    return 0
