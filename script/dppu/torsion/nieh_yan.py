"""
Nieh-Yan Variant Enum
=====================

Defines the Nieh-Yan topological term variants.

Physical Background:
    The Nieh-Yan term is a topological invariant in 4D manifolds with torsion:
        N = dT^a ∧ e_a - T^a ∧ T_a ∧ e^b ∧ e_b  (differential form notation)

    In our framework, we decompose it into two parts:
    - TT term: (1/4) ε^{abcd} T^e_{ab} T_{ecd}
    - Ree term: (1/4) ε^{abcd} R_{abcd}

    The full Nieh-Yan is: N_full = N_TT - N_Ree

References:
    - Nieh & Yan (1982): J. Math. Phys. 23, 373
    - Chandia & Zanelli (1997): Phys. Rev. D 55, 7580

Author: Muacca
"""

from enum import Enum


class NyVariant(Enum):
    """
    Nieh-Yan variant selection enumeration.

    Attributes:
        TT: TT-only variant.
            Uses only the TT term: (1/4) ε^{abcd} T^e_{ab} T_{ecd}
            Physical: Pure torsion-torsion coupling

        REE: Ree-only variant.
            Uses only the Ree term: (1/4) ε^{abcd} R_{abcd}
            Physical: Curvature contribution (Euler density related)

        FULL: Full Nieh-Yan variant.
            Uses N = N_TT - N_Ree
            Physical: Complete topological term

    Examples:
        >>> variant = NyVariant.FULL
        >>> variant.value
        'FULL'
        >>> NyVariant['TT']
        <NyVariant.TT: 'TT'>
    """
    TT = "TT"      # TT-only: (1/4) ε^{abcd} T^e_{ab} T_{ecd}
    REE = "REE"    # Ree-only: (1/4) ε^{abcd} R_{abcd}
    FULL = "FULL"  # Full: N_TT - N_Ree

    def __str__(self) -> str:
        return self.value

    @property
    def description(self) -> str:
        """Human-readable description of the variant."""
        descriptions = {
            NyVariant.TT: "TT-only (torsion-torsion)",
            NyVariant.REE: "Ree-only (curvature)",
            NyVariant.FULL: "Full Nieh-Yan (TT - Ree)",
        }
        return descriptions[self]

    @property
    def includes_tt(self) -> bool:
        """True if this variant includes the TT term."""
        return self in (NyVariant.TT, NyVariant.FULL)

    @property
    def includes_ree(self) -> bool:
        """True if this variant includes the Ree term."""
        return self in (NyVariant.REE, NyVariant.FULL)
