[![GitHub license](https://img.shields.io/github/license/Inazuma110/HCBE.jl)](https://github.com/Inazuma110/HCBE.jl/blob/master/LICENSE)
![CI](https://github.com/Inazuma110/HCBE.jl/workflows/CI/badge.svg)

# HCBE

## Overview
HCBE method is one of the hypergraph clustering algorithms.
For more details on the algorithm, please refer to the following paper.

- [エッジ追加に基づくハイパーグラフクラスタリングの特性評価](https://www.jstage.jst.go.jp/article/pjsai/JSAI2020/0/JSAI2020_2P5GS305/_article/-char/ja/)
- [ハイパーグラフクラスタリングにおけるエッジ追加順序の比較](https://www.ipsj.or.jp/event/fit/fit2020/FIT2020_program/data/html/abstract/D-004.html)
- [スター展開を利用したハイパーグラフのソフトクラスタリング](https://jsai.ixsq.nii.ac.jp/ej/index.php?active_action=repository_view_main_item_snippet&page_id=13&block_id=23&pn=1&count=20&order=16&lang=japanese&creator=ito%20shuta)

An English paper summarizing these three papers will be published at iiWAS2020.

## Install
```sh
(REPL)]add https://github.com/Inazuma110/HCBE.jl
```
or
```jl
using Pkg
Pkg.add(url="https://github.com/Inazuma110/HCBE.jl")
```

## Documentation
[Latest documentation](https://inazuma110.github.io/HCBE.jl/dev/)

## See Also
Julia packages providing other clustering methods:

- [Clustering.jl](https://github.com/JuliaStats/Clustering.jl)
- [QuickShiftClustering.jl](https://github.com/rened/QuickShiftClustering.jl)
- [SpectralClustering.jl](https://github.com/lucianolorenti/SpectralClustering.jl)
