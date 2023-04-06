""" Python module for calculating Feldman-Cousins confidence intervals."""

try:
    import FCpy.FCpy as FC
except ModuleNotFoundError:
    import FCpy as FC

__version__ = '0.1.1'
__pkgname__ = 'FCpy'
__author__ = 'Evan Groopman'
__url__ = 'https://github.com/usnistgov/FCpy'
__license__ = 'NIST Software'
__copyright__ = ''
__description__ = __doc__

_changelog = """
Version 0.1.1:
    - Add option to fix minor pathology at fixed n0 with different b.
    - Add optional dask dependency to perform this search faster.

Version 0.1.0:
    - First published version.
    - Primary emphasis on FC_poission() function.
"""