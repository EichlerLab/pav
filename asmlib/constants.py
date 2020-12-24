"""
Program constants
"""

# Version constants
VERSION_MAJOR = 0
VERSION_MINOR = 9
VERSION_DEV = 1


def get_version_string():
    """
    Get PAV version as a string.

    :return: PAV version string.
    """
    return '{}.{}.{}'.format(VERSION_MAJOR, VERSION_MINOR, VERSION_DEV)
