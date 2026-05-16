# ログ出力規約

このドキュメントは、`scripts/` 以下の 検証スクリプトにおけるログ出力のアーキテクチャと規約を定義します。

---

## 1. 二層構造

ログシステムは役割の異なる独立した二つの層で構成されます。

| 層 | API | 使用箇所 |
|----|-----|---------|
| **エンジン層** | `ComputationLogger` / `NullLogger` | `dppu/` ライブラリ内部 |
| **スクリプト層** | `setup_log` / `teardown_log` | `scripts/` 以下の検証スクリプト |

この二層は完全に独立しています。  
検証スクリプトは、ユーザー向け出力をスクリプト層で扱います。エンジンを直接構築する場合は、エンジン既定の `NullLogger` に依存するか `NullLogger()` を明示的に渡し、`ComputationLogger` は直接使いません。

図生成スクリプトは意図的な例外です。画像ファイルを `--output` で生成することが主目的であり、実行ログは作成しません。

---

## 2. エンジン層（`ComputationLogger` / `NullLogger`）

実装: `dppu/utils/logger.py`

- `ComputationLogger`: パイプラインの各ステップのタイミングを専用ログファイルに書き込む。
- `NullLogger`: 全てのログ呼び出しを無視するノーオップ。

**規則**: 検証スクリプトが `UnifiedEngine` 等を構築する際は `NullLogger()` に依存します。  
スクリプトの出力はスクリプト層で完結させます。

```python
from dppu.utils.logger import NullLogger

eng = UnifiedEngine(cfg, NullLogger())
```

---

## 3. スクリプト層（`setup_log` / `teardown_log`）

実装: `dppu/utils/tee_logger.py`

`setup_log` は `sys.stdout` を `_Tee` オブジェクトに差し替え、以降の全 `print()` 呼び出しをコンソールとログファイルの両方に書き込みます。  
`teardown_log` は `sys.stdout` を元に戻してログファイルを閉じ、サマリ行を出力します。

### 3.1 `teardown_log` が出力するサマリ

```
（空行）
Computation time: X.Xs
（空行）
*Script: scripts/<サブディレクトリ>/<スクリプト名>.py
```

既に独自でタイミングや `*Script:` 行を出力しているスクリプトは、二重出力を避けるためそれらの `print` を削除してください。

### 3.2 ログファイルの出力先

| スクリプト種別 | `log_dir` 引数 |
|--------------|--------------|
| argparse 検証スクリプト | `args.output_dir` |
| 非 argparse スクリプト | `"output"`（デフォルト）|

指定ディレクトリ（存在しない場合は自動作成）に `<スクリプト名>_YYYYMMDD_HHMMSS.log` の形式で保存されます。

### 3.3 環境変数によるオーバーライド

| 変数名 | 効果 |
|--------|------|
| `DPPU_LOG_DIR` | `log_dir` をこのパスで上書きする |
| `DPPU_LOG_STDOUT` | `0` / `false` にするとコンソール出力を抑制（ログファイルのみ）|

---

## 4. 使用パターン

### パターン A — argparse スクリプト（`main()` あり）

```python
import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from dppu.utils.tee_logger import setup_log, teardown_log

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()

    setup_log(__file__, log_dir=args.output_dir)

    # ... print() による出力 ...

    teardown_log()
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
```

### パターン B — 非 argparse スクリプト（`main()` あり）

```python
if __name__ == "__main__":
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
    from dppu.utils.tee_logger import setup_log, teardown_log
    setup_log(__file__, log_dir='output')
    main()
    teardown_log()
```

`import` を `if __name__` ブロック内に置くことで、スクリプトをモジュールとして import した場合に `sys.path` の変更や dppu の import が実行されないようにします。

### パターン C — モジュールレベルスクリプト（`main()` なし）

```python
import sys, os
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')  # エンコーディング修正を先に実施

import time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from dppu.utils.tee_logger import setup_log, teardown_log
setup_log(__file__, log_dir='output')   # エンコーディング修正の後に呼ぶ

# ... モジュールレベルでの計算と print() 出力 ...

teardown_log()   # ファイル末尾
```

**重要**: `setup_log` は Windows エンコーディング修正（`sys.stdout = io.TextIOWrapper(...)`）の**後**に呼び出します。  
これにより、`_Tee` が修正済みの UTF-8 stdout をラップする順序が保証されます。

---

## 5. 出力規約

- スクリプト内の出力メソッドは `print()` **のみ**を使用します。stdlib の `logging` モジュールは使用しません。
- Paper05 検証スクリプトのセクション区切りには `"=" * 72` を使用します。
- ステップヘッダーは `"--- Step N: ..."` 等の自由形式で記述します。
- `scripts/visualize/` は `[OK] wrote <path>` のような完了行のみを出力してよく、`setup_log` は不要です。

---

## 6. 関連ファイル

| ファイル | 役割 |
|---------|------|
| `dppu/utils/tee_logger.py` | `setup_log`、`teardown_log`、`LogTee` |
| `dppu/utils/logger.py` | `ComputationLogger`、`NullLogger` |
