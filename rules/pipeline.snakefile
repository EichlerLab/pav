"""
Pipeline control functions.
"""

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

    local_config = pavlib.pipeline.get_override_config(config, wildcards.asm_name, ASM_TABLE)

    if key is None:
        return local_config

    if key not in local_config:
        if default_val is None and not default_none:
            raise RuntimeError('Configuration does not include key "{}" for asm_name "{}" with no default specified'.format(key, wildcards.asm_name))

        return default_val

    return local_config[key]
