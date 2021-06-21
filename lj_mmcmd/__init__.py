"""
lj-mmcmd
A simple implementation of Metropolis Monte Carlo simulations and Molecular dynamics simulations for Lennard Jones particles
"""

# Add imports here
from .scripts import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
