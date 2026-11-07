# 計算補足資料

[English](README.md) · [規約](CONVENTIONS.md) · [旧規約との比較](legacy_conventions_ja.md)

`paper06.py recompute --output-dir NEW` は式と図のデータを再導出する。`paper06.py render --data-dir DATA --output-dir NEW` は保存データから図を描く。直下READMEの例のように、repository直下からは `script/` を付けて実行する。両コマンドとも新規出力先を要求し、既存ディレクトリを上書きしない。Python 3.11以降で `requirements.txt` をインストールする。

| Module | 内容 |
| --- | --- |
| `geometry` | 疎な外積形式、class A接続、torsion・曲率・NY・CS・transgression |
| `derive` | 記号恒等式、bulk境界分離、Euler式、拘束、四幾何への特殊化 |
| `shape`, `torus` | shape変分、混合錐、torusの厳密族 |
| `invariants` | 二次・三次縮約、Ricci逆問題 |
| `evolution` | 対角真空ODE、Taylor jet、12点の数値対照 |
| `frames`, `canonical` | 大域frame、端点の正準対応 |
| `components`, `identities` | Koszul接続による成分導出、恒等式検査 |
| `results` | 厳密値の保存、動的制限・場の方程式・曲線データの再計算 |
| `render` | 保存JSON/CSVを読み、PNG/PDFの6図を描画 |

`geometry.class_a()` は形式と全添字を下げた `R[a,b,c,d]` を返す。`derive.derive_geometry()` は幾何と導出したスカラー・端点式、`derive.derive_reduction(geometry_data)` は座標時間の作用、`derive.derive_classification(geometry_data)` は固定空間計量への特殊化を返す。入力記号と成立範囲は規約案内に示す。一般の異方的解や物理スペクトルを分類する関数ではない。

保存データ描画では導出を実行しない。記号再計算ではprofileを制限する前に接続を構成・微分する。本文の完成式へ数値を代入するだけの処理ではない。同梱ファイルと記載したPython依存で実行できる。

## 実行環境と数値対照

検証環境：Python 3.12.10、SymPy 1.14.0、NumPy 2.3.5、SciPy 1.16.3、Matplotlib 3.10.7。固定した依存関係はPython 3.11以降を必要とする。検証時の記号計算backendを再現するには、Pythonの起動前に `SYMPY_GROUND_TYPES=python` を設定する（PowerShellでは `$env:SYMPY_GROUND_TYPES='python'`、POSIX shellでは `export SYMPY_GROUND_TYPES=python`）。

`recompute` はclass A接続とtorsion消去後の独立なshape作用を構成し、新しく導出した曲率から縮約を求め、対角ODEをSciPyのDOP853で積分する。数値対照は指定した12点だけで、q0=1、Nil・SolのH0の両符号、時刻0.01、0.02、0.04、rtol=2e-13、atol=2e-15、max_step=0.002を用いる。これは点での対照であり、区間全体の誤差保証ではない。記号式は厳密であり、浮動小数点の出力は環境によってわずかに異なり得る。

`render` はJSON/CSVを読み、数式文字列の数値化だけを行う。曲率導出・ODE積分は実行しない。本文と付録の論証は実行を前提としない。

三フォルダだけで実行できる。両コマンドとも既存の出力先を拒否する。
