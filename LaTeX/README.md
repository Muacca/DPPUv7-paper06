# Build the English manuscript

Use a standard LaTeX distribution with the packages listed in `main.tex` (including amsmath, graphicx, tabularx, booktabs, hyperref and xurl). From this directory run:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The second run resolves equation and figure references. The figures are supplied as PDF and PNG pairs. To replace them with newly rendered figures, first use the script guide and compare the results. The manuscript uses stable equation identifiers shared with the Japanese Markdown. `main.pdf` is the build output; the repository's reading copy is `DPPUv7-paper06.pdf` one directory above.

日本語：このフォルダで上のコマンドを2回実行する。図は同梱されており、組版のために再計算する必要はない。式の識別子は日本語Markdownと共通である。

Version 2, revised 2026-11-07: eight sections, Appendices A–F and six figures.
