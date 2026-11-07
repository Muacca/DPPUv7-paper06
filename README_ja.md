# Lorentzian Einstein–Cartan Reduction on Four Homogeneous Geometries: Curvature, Boundary Terms, and Reduced Orbits

Muacca — Version 2、2026-11-07。初版：2026-05-16。

[English](README.md) · [英語PDF](DPPUv7-paper06.pdf) · [LaTeX](LaTeX/main.tex) · [版間変更](CHANGELOG.md)

固定class A空間幾何に、一つのscaleと二つの幾何学的torsion振幅を与え、bulk作用、符号付き境界完成、動的曲率量を導く。T3、Nil、Sol、S3のコンパクト参照幾何を比較し、曲率指標χで決まる縮約軌道・局所縮約量と、追加の幾何情報が必要なtracefree方程式・絶対端点規格化を区別する。Nil・Solの縮約軌道は完全真空解とは区別して扱う。

torsion消去後のshape拡張として、torus混合接方向の厳密解族とNil・Solの短時間発展を導く。選定二次・三次縮約に必要なχ・J2・J3の情報、大域frameと端点代表の正準対応を扱う。[参考文献](REFERENCES.md)。

| 本文 | 付録 |
| --- | --- |
| [1 導入](DPPUv7-paper06_sec01.md) | [A 構造定数と成分規約](DPPUv7-paper06_appA.md) |
| [2 幾何と変分領域](DPPUv7-paper06_sec02.md) | [B 曲率恒等式の導出](DPPUv7-paper06_appB.md) |
| [3 Bulk作用と境界完成](DPPUv7-paper06_sec03.md) | [C 変分と接続方程式](DPPUv7-paper06_appC.md) |
| [4 補助場拘束と軌道](DPPUv7-paper06_sec04.md) | [D shape変分と局所真空発展](DPPUv7-paper06_appD.md) |
| [5 共通scaleとshape発展](DPPUv7-paper06_sec05.md) | [E 添付資料](DPPUv7-paper06_appE.md) |
| [6 動的曲率と端点量](DPPUv7-paper06_sec06.md) | [F 改訂注記](DPPUv7-paper06_appF.md) |
| [7 曲率情報と幾何の表示](DPPUv7-paper06_sec07.md) | |
| [8 結論](DPPUv7-paper06_sec08.md) | |

## 再現

添付は `script`、`data`、`LaTeX` の三フォルダで構成する。[計算案内](script/README_ja.md)、[保存データ](data/revised/README.md)、[規約](script/CONVENTIONS.md)、[組版手順](LaTeX/README.md)を参照。Python 3.11以降を使い、出力先は新規ディレクトリを指定する。

```sh
python -m pip install -r script/requirements.txt
python -B script/paper06.py render --data-dir data/revised --output-dir rendered
python -B script/paper06.py recompute --output-dir recomputed
python -B script/paper06.py render --data-dir recomputed --output-dir rendered_recomputed
```

`render` は保存データだけを読む。`recompute` は接続を微分し、曲率・作用・変分・shape族・不変量・短時間数値対照・端点量を再導出する。他のPaperリポジトリへの依存はない。以前の出力との対応は[規約比較](script/legacy_conventions_ja.md)に示す。v1に同梱した汎用の `dppu` パッケージと計算ログは本版に含まない。

## 引用とライセンス

Muacca, *Lorentzian Einstein–Cartan Reduction on Four Homogeneous Geometries: Curvature, Boundary Terms, and Reduced Orbits*, v2 (2026)。既存の[Zenodo record](https://doi.org/10.5281/zenodo.20228881)については、取得した資料の版・日付を確認して引用する。[Repository](https://github.com/Muacca/DPPUv7-paper06) · [MIT License](LICENSE.txt)。連絡先：muacca@dmwp.jp。
