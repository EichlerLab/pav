"""
Pipeline control functions.
"""

import pavlib

global ASM_TABLE
global config

def get_override_config(asm_name):
    return pavlib.config.get_override_config(asm_name, config, ASM_TABLE)

def get_config(key, wildcards):
    """
    Get a config object that might be modified by CONFIG parameters in the assembly table.

    If "key" is None, the full config dictionary is returned. If "key" is defined, then the value of config for
    that key is returned with an optional default value.

    :param wildcards: Rule wildcards.
    :param key: Key of the value to get from config.

    :return: Config object. Original global config, if unmodified, or a modified copy of it.
    """

    return pavlib.config.get_config(key, wildcards, config, ASM_TABLE)
