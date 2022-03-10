# spectrum-utils
a package for spectroscopy

# Install
from console
~~~
pip install git+https://github.com/AtsushiNi/spectrumutils.git
~~~

from jupyter notebook/lab
~~~
!pip install git+https://github.com/AtsushiNi/spectrumutils.git
~~~

# Usage
## Calculate spectrum amplitude
~~~
from spectrumutils import fit

"""
ガウス関数でフィッテング
ピーク数は1~4個を自動で判別
x, y: フィッティングしたいデータ
range: フィッティングしたい範囲のx座標
thres: どれだけ小さいピークまで検出するか(0~1, default=0.1)
"""
result = fit.fit(x, y, range=[601.1, 602.3], thres=0.1)
~~~

## Create boltzmann plot
~~~
from spectrumutils import fulcher

"""
ボルツマンプロットの作成
発光上準位の占有率(population)の計算
data: fulcher-a発光強度のデータ(計測できなかったものは0埋め)
v: 振動準位
"""
result = fulcher.boltzmannplot([data_v0, data_v1, data_v2], v=[0,1,2])
~~~

# Sample
see the 'sample' folder
