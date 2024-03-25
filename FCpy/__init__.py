""" Python module for calculating Feldman-Cousins confidence intervals."""

try:
    from FCpy.FCpy import FC
except ModuleNotFoundError:
    import FCpy.FC as FC
except ModuleNotFoundError:
    import FC

__version__ = '0.1.3'
__pkgname__ = 'FCpy'
__author__ = 'Evan Groopman'
__url__ = 'https://github.com/usnistgov/FCpy'
__license__ = 'NIST Software'
__copyright__ = ''
__description__ = __doc__

_changelog = """
Version 0.1.3:
    - Add FC_poisson_sensitivity(...) for calculating sensitivity.
    
Version 0.1.2:
    - Add convenience function FC_poisson_confband() for calculating confidence band for a range of n.
    - Convenience function FC_poisson_list() for calculating a list of n and t values.
    - Add optional user-defined limits in FC_poisson() for range of mu to calculate (mumin, mumax).
      This can increase computation speed over the conservative auto bounds.
      Add warnings for CL_high/CL_low hitting the user-defined mumax/mumin.

Version 0.1.1:
    - Add option to fix minor pathology at fixed n0 with different b.
    - Add optional dask dependency to perform this search faster.

Version 0.1.0:
    - First published version.
    - Primary emphasis on FC_poission() function.
"""