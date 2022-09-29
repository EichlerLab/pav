"""
Program constants
"""

# Version constants
VERSION_MAJOR = 2
VERSION_MINOR = 1
VERSION_DEV = 0


def get_version_string():
    """
    Get PAV version as a string.

    :return: PAV version string.
    """
    return '{}.{}.{}'.format(VERSION_MAJOR, VERSION_MINOR, VERSION_DEV)


############
### Call ###
############

### Haplotype Merging ###

# Default merge for INS/DEL/INV
MERGE_PARAM_INSDELINV = 'nr::ro(0.5):szro(0.5,200,2):match'
MERGE_PARAM_SNV = 'nrsnv::exact'

MERGE_PARAM_DEFAULT = {
    'ins': MERGE_PARAM_INSDELINV,
    'del': MERGE_PARAM_INSDELINV,
    'inv': MERGE_PARAM_INSDELINV,
    'snv': MERGE_PARAM_SNV
}


####################
### Return Codes ###
####################

# When scripts/density.py fails to find an inversion but not because of an error, this value is returned to indicate
# that the pipeline should go ahead without crashing. Return codes other than 0 and ERR_INV_FAIL indicate a problem or
# that should be corrected.
ERR_INV_FAIL = 125
