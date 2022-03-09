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
"""
range = [601.1, 602.3]
result = fit.fit(x, y, range=range)
~~~

# Sample
see the 'sample' folder
