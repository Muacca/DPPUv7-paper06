"""
Checkpoint Management
=====================

Provides save/restore functionality for computation state.

This allows:
- Resuming interrupted computations
- Debugging specific steps
- Reproducibility

Author: Muacca
"""

import pickle
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Tuple


class CheckpointManager:
    """
    Manages saving and loading of computation checkpoints.

    Checkpoints are saved as pickle files containing:
    - step_id: The step at which the checkpoint was created
    - timestamp: When the checkpoint was created
    - data: The engine's data dictionary
    - metadata: Additional information (mode, variant, etc.)

    Example:
        >>> ckpt = CheckpointManager("checkpoints", enabled=True)
        >>> ckpt.save("E4.7", engine.data, {'mode': 'MX'})
        >>> data, meta = ckpt.load("E4.7")
    """

    def __init__(
        self,
        output_dir: str = "checkpoints",
        enabled: bool = False
    ):
        """
        Initialize the checkpoint manager.

        Args:
            output_dir: Directory to store checkpoint files
            enabled: If False, save/load operations are no-ops
        """
        self.output_dir = Path(output_dir)
        self.enabled = enabled

        if self.enabled:
            self.output_dir.mkdir(parents=True, exist_ok=True)

    def save(
        self,
        step_id: str,
        data: Dict[str, Any],
        metadata: Optional[Dict[str, Any]] = None
    ) -> Optional[Path]:
        """
        Save a checkpoint.

        Args:
            step_id: Step identifier (e.g., "E4.7")
            data: Engine data dictionary to save
            metadata: Additional metadata

        Returns:
            Path to saved checkpoint file, or None if disabled
        """
        if not self.enabled:
            return None

        checkpoint = {
            'step_id': step_id,
            'timestamp': datetime.now().isoformat(),
            'data': data,
            'metadata': metadata or {}
        }

        # Convert step_id to filename (E4.7 -> checkpoint_E4_7.pkl)
        safe_step_id = step_id.replace('.', '_')
        filepath = self.output_dir / f"checkpoint_{safe_step_id}.pkl"

        try:
            with open(filepath, 'wb') as f:
                pickle.dump(checkpoint, f, protocol=pickle.HIGHEST_PROTOCOL)
            return filepath
        except Exception as e:
            print(f"  [ERROR] Failed to save checkpoint: {e}")
            return None

    def load(self, step_id: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """
        Load a checkpoint.

        Args:
            step_id: Step identifier (e.g., "E4.7")

        Returns:
            Tuple of (data, metadata)

        Raises:
            FileNotFoundError: If checkpoint doesn't exist
        """
        safe_step_id = step_id.replace('.', '_')
        filepath = self.output_dir / f"checkpoint_{safe_step_id}.pkl"

        if not filepath.exists():
            raise FileNotFoundError(f"Checkpoint {step_id} not found at {filepath}")

        with open(filepath, 'rb') as f:
            checkpoint = pickle.load(f)

        print(f"  [OK] Loaded checkpoint: {step_id}")
        return checkpoint['data'], checkpoint['metadata']

    def exists(self, step_id: str) -> bool:
        """
        Check if a checkpoint exists.

        Args:
            step_id: Step identifier

        Returns:
            True if checkpoint file exists
        """
        if not self.enabled:
            return False

        safe_step_id = step_id.replace('.', '_')
        filepath = self.output_dir / f"checkpoint_{safe_step_id}.pkl"
        return filepath.exists()

    def list_checkpoints(self) -> list:
        """
        List all available checkpoints.

        Returns:
            List of step_ids with available checkpoints
        """
        if not self.enabled or not self.output_dir.exists():
            return []

        checkpoints = []
        for filepath in self.output_dir.glob("checkpoint_*.pkl"):
            # Extract step_id from filename
            name = filepath.stem  # checkpoint_E4_7
            step_id = name.replace("checkpoint_", "").replace("_", ".")
            checkpoints.append(step_id)

        return sorted(checkpoints)

    def clear(self) -> int:
        """
        Remove all checkpoints.

        Returns:
            Number of checkpoints removed
        """
        if not self.enabled or not self.output_dir.exists():
            return 0

        count = 0
        for filepath in self.output_dir.glob("checkpoint_*.pkl"):
            filepath.unlink()
            count += 1

        return count
