# Reduced-Sector chi-Universality in Lorentzian EC+NY Minisuperspace: Topology-Robust Admissibility, P-Channel Diagnostics, and Reduced Vacuum Orbits

[日本語 README](README_ja.md)

## Paper Contents

- [Official PDF (English, built from LaTeX)](DPPUv7-paper06.pdf)
- [Draft in Japanese Markdown](DPPUv7-paper06_sec01.md)

## Directory Structure

- `LaTeX/` — LaTeX manuscript and figures
  - `sections/` — TeX files for each section
  - `appendices/` — TeX files for appendices
  - `figures/` — Figures
- `script/` — Data processing and verification scripts
  - For script details, see [script/README.md](script/README.md)
- `data/` — Numerical data and computation logs

## Building the LaTeX PDF

### Basic Build Command

Navigate to the LaTeX directory and run pdflatex.
**At least two compilation passes are required** to resolve cross-references (`\ref`, `\label`, etc.) correctly.

- **First pass**: writes label information to the `.aux` file
- **Second pass**: reads and resolves references from the `.aux` file

If references appear as `??`, run the compilation once more.

```bash
cd LaTeX
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

### Output Files

A successful build produces:

- `main.pdf` — Final PDF file
- `main.aux` — Auxiliary file (cross-reference data)
- `main.log` — Compilation log
- `main.out` — hyperref outline information

### Troubleshooting

#### On Error

Check the log file:
```bash
grep -i "error\|warning" main.log | tail -20
```

#### Clean Build

Delete all auxiliary files and rebuild from scratch:
```bash
cd LaTeX
rm -f *.aux *.log *.out *.synctex.gz *.fdb_latexmk *.fls
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

## Author

Muacca

## License

See the LICENSE file in the repository root.
