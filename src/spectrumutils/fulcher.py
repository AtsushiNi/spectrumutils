import numpy as np
from importlib_resources import files, as_file
import xarray as xr
import matplotlib.pyplot as plt

# 発光強度データからボルツマンプロットを作成
def boltzmannplot(data, v):
    # 振動準位ごとに
    result = []
    for (d, dv) in zip(data, v):
        index = np.nonzero(d)
        N_numbers = index[0]+1
        amplitudes = d[index]

        population = np.zeros(amplitudes.size)
        rot_energy = E_d_rot(dv, N_numbers)

        for j, N in enumerate(N_numbers):
            population[j] = amplitudes[j] * fulcher_wavelength().sel(dv=dv,dN=N) **4 / (2*N+1)/ g_as(N)

        plt.plot(rot_energy, population, '--x')
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
    franck_condon = files('src.spectrumutils.data.fulcher').joinpath('fulcher_wavelength.nc')
    with as_file(franck_condon) as f:
        r = xr.open_dataarray(f)
    return r

# Franck-Condin因子(https://inis.iaea.org/collection/NCLCollectionStore/_Public/37/088/37088524.pdf)
def franck_condon_X_to_d():
    franck_condon = files('src.spectrumutils.data.fulcher').joinpath('franck_condon_factor.txt')
    with as_file(franck_condon) as f:
        r = np.loadtxt(f)
    return r

def franck_condon_d_to_a():
    franck_condon = files('src.spectrumutils.data.fulcher').joinpath('franck_condon_factor_d-to-a.csv')
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
