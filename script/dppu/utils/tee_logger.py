"""
Tee Logger
==========

Utility that writes stdout to both the console and a log file.
The DPPU_LOG_STDOUT environment variable can suppress console output
so that output goes to the log file only.

Usage::

    from dppu.utils.tee_logger import setup_log, teardown_log

    log_path = setup_log(__file__, log_dir="logs")
    print("This output is written to both console and log file.")
    teardown_log()

Or as a context manager::

    from dppu.utils.tee_logger import LogTee

    with LogTee(__file__, log_dir="logs") as log_path:
        print("Logging in progress.")
"""

import sys
import os
import time as _time
from datetime import datetime


class _Tee:
    """Wrapper that writes simultaneously to console and a file."""

    def __init__(self, file_path: str, encoding: str = "utf-8", echo_console: bool = True):
        self._console = sys.stdout
        self._file = open(file_path, "w", encoding=encoding, buffering=1)
        self._echo_console = echo_console

    def write(self, data: str) -> int:
        if self._echo_console:
            try:
                self._console.write(data)
            except UnicodeEncodeError:
                enc = getattr(self._console, "encoding", "ascii") or "ascii"
                self._console.write(data.encode(enc, errors="replace").decode(enc))
        self._file.write(data)
        return len(data)

    def flush(self):
        if self._echo_console:
            self._console.flush()
        self._file.flush()

    def close(self):
        if self._file and not self._file.closed:
            self._file.close()

    # Delegate readline / isatty to the console.
    def isatty(self) -> bool:
        return self._console.isatty() if self._echo_console else False


_active_tee: "_Tee | None" = None
_start_time: "float | None" = None
_script_abs_path: "str | None" = None


def _env_flag(name: str, default: bool) -> bool:
    """Parse a boolean environment flag."""
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() not in {"0", "false", "no", "off"}


def setup_log(script_file: str, log_dir: str = "logs") -> str:
    """
    Open a log file and replace stdout with a Tee.

    Parameters
    ----------
    script_file : str
        Path to the calling script (pass ``__file__``).
        The basename is used as the log file prefix.
    log_dir : str
        Log directory. Resolved relative to the process cwd if not absolute.

    Returns
    -------
    str
        Path of the created log file.
    """
    global _active_tee, _start_time, _script_abs_path

    log_dir_override = os.environ.get("DPPU_LOG_DIR")
    if log_dir_override:
        log_dir = log_dir_override

    # Resolve log directory: if not absolute, treat as relative to cwd.
    if not os.path.isabs(log_dir):
        log_dir = os.path.join(os.getcwd(), log_dir)
    os.makedirs(log_dir, exist_ok=True)

    base = os.path.splitext(os.path.basename(script_file))[0]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(log_dir, f"{base}_{timestamp}.log")

    _script_abs_path = os.path.abspath(script_file)
    _start_time = _time.time()

    echo_console = _env_flag("DPPU_LOG_STDOUT", True)
    _active_tee = _Tee(log_path, echo_console=echo_console)
    sys.stdout = _active_tee
    return log_path


def teardown_log():
    """Restore stdout to its original state and close the log file."""
    global _active_tee, _start_time, _script_abs_path
    if _active_tee is not None:
        elapsed = _time.time() - _start_time if _start_time is not None else 0.0

        # Compute relative script path: relative to the "script/" root
        # (script files live at script/scripts/<subdir>/name.py, lib is script/)
        try:
            lib_dir = os.path.normpath(
                os.path.join(os.path.dirname(_script_abs_path), "..", ".."))
            rel_path = os.path.relpath(_script_abs_path, lib_dir)
            rel_path = rel_path.replace(os.sep, "/")
        except Exception:
            rel_path = os.path.basename(_script_abs_path)

        print()
        print(f"Computation time: {elapsed:.1f}s")
        print()
        print(f"*Script: {rel_path}")

        sys.stdout = _active_tee._console
        _active_tee.close()
        _active_tee = None
        _start_time = None
        _script_abs_path = None


class LogTee:
    """Context manager wrapper around setup_log / teardown_log."""

    def __init__(self, script_file: str, log_dir: str = "logs"):
        self._script_file = script_file
        self._log_dir = log_dir
        self.log_path: str = ""

    def __enter__(self) -> str:
        self.log_path = setup_log(self._script_file, self._log_dir)
        return self.log_path

    def __exit__(self, *_):
        teardown_log()
