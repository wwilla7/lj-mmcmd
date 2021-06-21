"""
Unit and regression test for the lj_mmcmd package.
"""

# Import package, test suite, and other packages as needed
import lj_mmcmd
import pytest
import sys

def test_lj_mmcmd_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "lj_mmcmd" in sys.modules
