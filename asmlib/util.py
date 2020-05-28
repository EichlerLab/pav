"""
General utility functions.
"""

def as_bool(val):
    """
    Translate value as a boolean. If `val` is boolean, return `val`. If val is a string, `True` if lower-case string
    is "true", "1", "yes", "t", or "y". `False` if lower-case string is "false", "0", "no", "f", or "n". All other
    values throw a RuntimeError. If `val` is not bool or string, it is converted to a string and the rules above are
    applied.

    :param val: Value to interpret as a boolean.

    :return: `True` or `False` (see above).
    """

    if issubclass(val.__class__, bool):
        return val

    val = str(val).lower()

    if val in {'true', '1', 'yes', 't', 'y'}:
        return True

    if val in {'false', '0', 'no', 'f', 'n'}:
        return False

    raise RuntimeError('Cannot interpret as boolean value: {}'.format(val))
