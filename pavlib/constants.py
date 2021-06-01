"""
Program constants
"""

# Version constants
VERSION_MAJOR = 1
VERSION_MINOR = 1
VERSION_DEV = 0


def get_version_string():
    """
    Get PAV version as a string.

    :return: PAV version string.
    """
    return '{}.{}.{}'.format(VERSION_MAJOR, VERSION_MINOR, VERSION_DEV)
