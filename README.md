# Lorentzian Einstein–Cartan Reduction on Four Homogeneous Geometries: Curvature, Boundary Terms, and Reduced Orbits

Muacca — Version 2, revised 2026-11-07. First release: 2026-05-16.

[日本語](README_ja.md) · [English PDF](DPPUv7-paper06.pdf) · [LaTeX](LaTeX/main.tex) · [Changes from v1](CHANGELOG.md)

The paper derives the bulk action, signed endpoint completion and dynamic curvature quantities for a fixed class-A spatial geometry with one scale and two geometric torsion amplitudes. It compares compact representatives of T3, Nil, Sol and S3. The scalar-curvature parameter organizes reduced trajectories and selected local contractions; tracefree field equations and absolute endpoint normalizations retain additional geometric data. Nil and Sol reduced trajectories are distinguished from full vacuum solutions.

The torsionless shape extension gives an exact mixed-tangent cone on the torus and short-time Nil/Sol evolution. Selected quadratic and cubic invariants resolve the information in chi, J2 and J3; global frames and canonical endpoint representatives clarify the limits of structural comparisons.

[References](REFERENCES.md). The Japanese manuscript is organized as follows. The English version is in the PDF and LaTeX source.

| Section | Appendix |
| --- | --- |
| [1 Introduction](DPPUv7-paper06_sec01.md) | [A Structure and conventions](DPPUv7-paper06_appA.md) |
| [2 Geometry](DPPUv7-paper06_sec02.md) | [B Curvature derivations](DPPUv7-paper06_appB.md) |
| [3 Bulk and boundary action](DPPUv7-paper06_sec03.md) | [C Variations](DPPUv7-paper06_appC.md) |
| [4 Auxiliary constraints and trajectories](DPPUv7-paper06_sec04.md) | [D Shape and local vacuum evolution](DPPUv7-paper06_appD.md) |
| [5 Common scale and shape evolution](DPPUv7-paper06_sec05.md) | [E Supplement](DPPUv7-paper06_appE.md) |
| [6 Dynamic curvature and endpoints](DPPUv7-paper06_sec06.md) | [F Revision note](DPPUv7-paper06_appF.md) |
| [7 Curvature information and representations](DPPUv7-paper06_sec07.md) | |
| [8 Conclusions](DPPUv7-paper06_sec08.md) | |

## Reproduce

The supplement consists of `script`, `data` and `LaTeX`. See the [code guide](script/README.md), [saved-data guide](data/revised/README.md), [conventions](script/CONVENTIONS.md) and [build instructions](LaTeX/README.md). Use Python 3.11 or later and new output directories.

```sh
python -m pip install -r script/requirements.txt
python -B script/paper06.py render --data-dir data/revised --output-dir rendered
python -B script/paper06.py recompute --output-dir recomputed
python -B script/paper06.py render --data-dir recomputed --output-dir rendered_recomputed
```

`render` reads saved data only. `recompute` differentiates the connection and derives curvature, the action, variations, shape families, invariant contractions, short-time numerical checks and endpoint quantities. Its results are independent of other paper repositories. The [convention comparison](script/legacy_conventions.md) relates the present quantities to earlier outputs. The generic `dppu` package and calculation logs distributed with v1 are not part of this version.

## Citation and license

Muacca, *Lorentzian Einstein–Cartan Reduction on Four Homogeneous Geometries: Curvature, Boundary Terms, and Reduced Orbits*, v2 (2026). For the existing [Zenodo record](https://doi.org/10.5281/zenodo.20228881), check the edition and date of the artifact obtained. [Repository](https://github.com/Muacca/DPPUv7-paper06) · [MIT License](LICENSE.txt). Contact: muacca@dmwp.jp.
