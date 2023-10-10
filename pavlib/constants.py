"""
Program constants
"""

# Version constants
VERSION_MAJOR = 2     # Major change or a large number of minor changes
VERSION_MINOR = 3     # A significant change it PAV or a large number of small incremental changes
VERSION_DEV = 4       # Small changes, usually bug fixes
VERSION_PATCH = None  # Development and test versions, not for releases

# VERSION_PATCH is "None" for releases or set if:
#
# 1) Testing fixes before a release.
#
# 2) In rare cases when an older version of PAV is patched, the patch version is set. Changes to those versions are
# not synced with the main branch, but may contain fixes applied to later and replicated into an earlier branch.
# This is rarely used and only in special circumstances when a very large project is using a specific version but needs
# bug without significant changes to PAV.


def get_version_string():
    """
    Get PAV version as a string.

    :return: PAV version string.
    """

    return f'{VERSION_MAJOR}.{VERSION_MINOR}.{VERSION_DEV}' + (f'.{VERSION_PATCH}' if VERSION_PATCH is not None else '')


########
# Call #
########

# Default merge for INS/DEL/INV
MERGE_PARAM_INSDELINV = 'nr::ro(0.5):szro(0.5,200,2):match'
MERGE_PARAM_SNV = 'nrsnv::exact'

MERGE_PARAM_DEFAULT = {
    'ins': MERGE_PARAM_INSDELINV,
    'del': MERGE_PARAM_INSDELINV,
    'inv': MERGE_PARAM_INSDELINV,
    'snv': MERGE_PARAM_SNV
}


################
# Return Codes #
################

# When scripts/density.py fails to find an inversion but not because of an error, this value is returned to indicate
# that the pipeline should go ahead without crashing. Return codes other than 0 and ERR_INV_FAIL indicate a problem or
# that should be corrected.
ERR_INV_FAIL = 125
