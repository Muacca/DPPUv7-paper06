"""
Computation Logging
===================

Provides logging infrastructure for the computation pipeline.

Classes:
- ComputationLogger: Full-featured logger with file output
- NullLogger: No-op logger for silent operation

Usage:
    # With logging
    logger = ComputationLogger("computation.log")
    engine = make_engine(
        TopologyType.S3, mode, variant, logger=logger, enable_squash=True
    )

    # Without logging (silent)
    engine = make_engine(
        TopologyType.S3, mode, variant, logger=NullLogger(), enable_squash=True
    )

Author: Muacca
"""

from datetime import datetime
from pathlib import Path
from typing import Optional


class NullLogger:
    """
    No-operation logger that silently ignores all calls.

    Use this when you want to run computations without any output.
    Implements the same interface as ComputationLogger.
    """

    def step(self, step_id: str, description: str) -> str:
        """Log step start (no-op)."""
        return step_id

    def info(self, message: str, indent: int = 2) -> None:
        """Log info message (no-op)."""
        pass

    def warning(self, message: str) -> None:
        """Log warning message (no-op)."""
        pass

    def error(self, message: str) -> None:
        """Log error message (no-op)."""
        pass

    def success(self, message: str) -> None:
        """Log success message (no-op)."""
        pass

    def finalize(self) -> None:
        """Finalize logging (no-op)."""
        pass


class ComputationLogger:
    """
    Manages step-by-step logging with timestamps and file output.

    Features:
    - Automatic timestamping of steps
    - Duration tracking per step
    - Both console and file output
    - Total runtime tracking

    Example:
        >>> logger = ComputationLogger("my_computation.log")
        >>> logger.step("E4.1", "Setup parameters")
        >>> logger.info("Mode: MX")
        >>> logger.success("Setup complete")
        >>> logger.finalize()
    """

    def __init__(self, log_file: str = "dppu_engine.log"):
        """
        Initialize the logger.

        Args:
            log_file: Path to the log file
        """
        self.log_file = Path(log_file)
        self.current_step: Optional[str] = None
        self.step_start_time: Optional[datetime] = None
        self.total_start_time = datetime.now()

        # Initialize log file
        with open(self.log_file, 'w', encoding='utf-8') as f:
            f.write("=" * 70 + "\n")
            f.write("DPPU Engine Log\n")
            f.write("=" * 70 + "\n")
            f.write(f"Started: {self.total_start_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    def step(self, step_id: str, description: str) -> str:
        """
        Log the start of a new computation step.

        Args:
            step_id: Step identifier (e.g., "E4.1")
            description: Brief description of the step

        Returns:
            The step_id for chaining
        """
        # Log duration of previous step if any
        if self.step_start_time:
            duration = (datetime.now() - self.step_start_time).total_seconds()
            self._log(f"  Duration: {duration:.2f}s\n")

        self.current_step = step_id
        self.step_start_time = datetime.now()
        timestamp = self.step_start_time.strftime("%H:%M:%S")

        desc_str = description.strip() if description else "No description"

        msg = f"[{timestamp}] STEP {step_id}: {desc_str}"
        print("=" * 70)
        print(msg)
        print("=" * 70)
        self._log(msg + "\n")

        return step_id

    def info(self, message: str, indent: int = 2) -> None:
        """
        Log an informational message.

        Args:
            message: The message to log
            indent: Number of spaces to indent
        """
        prefix = " " * indent
        msg = f"{prefix}{message}"
        print(msg)
        self._log(msg + "\n")

    def warning(self, message: str) -> None:
        """
        Log a warning message.

        Args:
            message: The warning message
        """
        msg = f"  [WARNING] {message}"
        print(msg)
        self._log(msg + "\n")

    def error(self, message: str) -> None:
        """
        Log an error message.

        Args:
            message: The error message
        """
        msg = f"  [ERROR] {message}"
        print(msg)
        self._log(msg + "\n")

    def success(self, message: str) -> None:
        """
        Log a success message.

        Args:
            message: The success message
        """
        msg = f"  [SUCCESS] {message}"
        print(msg)
        self._log(msg + "\n")

    def _log(self, message: str) -> None:
        """Write message to log file."""
        with open(self.log_file, 'a', encoding='utf-8') as f:
            f.write(message)

    def finalize(self) -> None:
        """
        Finalize logging and report total runtime.

        Should be called at the end of computation.
        """
        # Log duration of last step
        if self.step_start_time:
            duration = (datetime.now() - self.step_start_time).total_seconds()
            self._log(f"  Duration: {duration:.2f}s\n")

        total_duration = (datetime.now() - self.total_start_time).total_seconds()
        msg = f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        msg += f"Total runtime: {total_duration:.1f}s ({total_duration/60:.1f} min)\n"
        print(msg)
        self._log(msg)
