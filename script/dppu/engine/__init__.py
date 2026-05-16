"""
DPPU Engine Layer
=================

Provides computation pipeline infrastructure and low-level geometry.

Modules:
- pipeline: BaseFrameEngine and step orchestration
- checkpoint: CheckpointManager for save/restore
- metric: Frame metric utilities
- levi_civita: Levi-Civita connection (Koszul formula)
- contortion: Contortion tensor from torsion
- ec_connection: Einstein-Cartan connection

Note: ComputationLogger and NullLogger live in utils.logger;
re-exported here for backward compatibility.
"""

from ..utils.logger import ComputationLogger, NullLogger
from .checkpoint import CheckpointManager

__all__ = [
    'ComputationLogger',
    'NullLogger',
    'CheckpointManager',
]
