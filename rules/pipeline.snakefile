"""
Pipeline control functions.
"""

import pavlib

global ASM_TABLE
global config

def get_config(wildcards, key=None, default_val=None, default_none=False):
    """
    Get a config object that might be modified by CONFIG parameters in the assembly table.

    If "key" is None, the full config dictionary is returned. If "key" is defined, then the value of config for
    that key is returned with an optional default value.

    :param wildcards: Rule wildcards.
    :param key: Key of the value to get from config.
    :param default_val: Default value.
    :param default_none: If True and default_val is None, return None instead of throwing an exception.

    :return: Config object. Original global config, if unmodified, or a modified copy of it.
    """

    return pavlib.pipeline.get_config(wildcards, config, ASM_TABLE, key=key, default_val=default_val, default_none=default_none)
