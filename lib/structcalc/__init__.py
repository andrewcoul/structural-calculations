"""structcalc: shared calculation engine for structural-calculations Colab notebooks.

US customary units throughout (kip, in, ft, psi/ksi, plf/psf).
Code basis: IBC 2024, ASCE 7-22, AISC 360-22, NDS 2024 / SDPWS 2021,
ACI 318-19.

These notebooks are design aids. A licensed engineer must review every
result before use.
"""

__version__ = "0.1.0"

from .report import Report  # noqa: F401
