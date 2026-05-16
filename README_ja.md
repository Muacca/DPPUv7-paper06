# Reduced-Sector chi-Universality in Lorentzian EC+NY Minisuperspace: Topology-Robust Admissibility, P-Channel Diagnostics, and Reduced Vacuum Orbits

[English README](README.md)

## 論文内容

- [正式な PDF（英語。LaTeX からビルドしたもの）](DPPUv7-paper06.pdf)
- [ドラフトの日本語版 Markdown](DPPUv7-paper06_sec01.md)

## ディレクトリ構造

- `LaTeX/` — LaTeX 原稿と図
  - `sections/` — 各セクションの TeX ファイル
  - `appendices/` — 付録の TeX ファイル
  - `figures/` — 画像
- `script/` — データ処理・検証スクリプト
  - スクリプトの詳細については [script/README_ja.md](script/README_ja.md) を参照ください。
- `data/` — 数値データと計算ログ

## LaTeX のビルド方法

### 基本的なビルドコマンド

LaTeX ディレクトリに移動して pdflatex を実行します。
LaTeX の相互参照（`\ref`、`\label` など）を正しく解決するため、**最低 2 回のコンパイルが必要**です。

- **1 回目**: `.aux` ファイルにラベル情報を書き込む
- **2 回目**: `.aux` ファイルから参照を読み込んで解決

参照が `??` と表示される場合は、もう一度コンパイルしてください。

```powershell
cd LaTeX
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

### 出力ファイル

ビルドが成功すると以下のファイルが生成されます：

- `main.pdf` — 最終的な PDF ファイル
- `main.aux` — 補助ファイル（相互参照情報）
- `main.log` — コンパイルログ
- `main.out` — hyperref アウトライン情報

### トラブルシューティング

#### エラーが出る場合

ログファイルを確認：
```powershell
cat main.log | Select-String -Pattern "Error|Warning" | Select-Object -Last 20
```

#### クリーンビルド

一度すべての補助ファイルを削除してから再ビルド：
```powershell
cd LaTeX
Remove-Item *.aux, *.log, *.out, *.synctex.gz, *.fdb_latexmk, *.fls -ErrorAction SilentlyContinue
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

## 著者

Muacca

## ライセンス

リポジトリ直下の LICENSE ファイルを参照してください。
