"""
Alignment parameters and helper functions.
"""

DEFAULT_ALIGNER = {
    1: 'minimap2',
    2: 'minimap2'
}

DEFAULT_ALIGNER_PARAMS = {
    1: {
        'minimap2': '-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5',
        'lra': ''
    },
    2: {
        'minimap2': '-x asm5'
    }
}

KNOWN_ALIGNERS = ['minimap2', 'lra']

def get_aligner(tier, config):
    """
    Get aligner.

    :param tier: Alignment tier.
    :param config: Pipeline config.

    :return: Name of aligner. Will pull from default values if not overridden by config.
    """

    # Get tier
    try:
        tier = int(tier)
    except ValueError:
        raise RuntimeError(f'Parameter "tier" must be an integer: {tier}')

    if tier < 1 or tier > 2:
        raise RuntimeError(f'Parameter "tier" is out of bounds: Expected value between 1 and 2: {tier}')

    # Get aligner
    aligner_config = f'aligner_tier{tier}'

    if tier == 1 and 'aligner' in config:

        # Allow "aligner" as an alias for "aligner_tier1", but only if one is defined or both are equal
        if aligner_config in config and str(config['aligner']).strip() != str(config[aligner_config]).strip():
            raise RuntimeError(f'Conflicting aligner parameters for tier {tier}: Both "aligner", and "{aligner_config}" are defined and are not the same parameter string.')

        aligner_config = 'aligner'

    if aligner_config not in config:
        if tier not in DEFAULT_ALIGNER:
            raise RuntimeError(f'No default aligner for tier {tier}')

        aligner = DEFAULT_ALIGNER[tier]
    else:
        aligner = config[aligner_config]

        if aligner not in KNOWN_ALIGNERS:
            raise RuntimeError(f'Unknown aligner: {aligner}: Known aligners are {", ".join(KNOWN_ALIGNERS)}')

    return aligner


def get_aligner_params(tier, config, aligner=None):
    """
    Get parameters for an alignment tier.

    :param tier: Alignment tier.
    :param config: Pipeline config.
    :param aligner: Alignment program name (minimap2, lra).

    :return: A string of parameters for the aligner. Will pull from default values if not overridden by config.
    """

    # Get tier
    try:
        tier = int(tier)
    except ValueError:
        raise RuntimeError(f'Parameter "tier" must be an integer: {tier}')

    if tier < 1 or tier > 2:
        raise RuntimeError(f'Parameter "tier" is out of bounds: Expected value between 1 and 2: {tier}')

    # Get aligner
    if aligner is None:
        aligner = get_aligner(tier, config)

    # Get aligner config parameter
    aligner_config = f'align_params_tier{tier}'

    if tier == 1 and aligner == 'minimap2' and 'minimap2_params' in config:

        # Backwards compatibility: Allow "minimap2_params" (or "lra_params" to be consistent) for compatibility with older PAV versions
        if aligner_config in config and str(config[f'{aligner}_params']).strip() != str(config[aligner_config]).strip():
            raise RuntimeError(f'Conflicting alignment parameters for tier {tier}: Both "{aligner}_params", and "{aligner_config}" are defined and are not the same parameter string.')

        aligner_config = f'{aligner}_params'

    if aligner_config in config:
        print(f'Getting alignments from config: {aligner_config}')
        align_params = config[aligner_config]
    else:
        print(f'Getting alignments from default')
        align_params = DEFAULT_ALIGNER_PARAMS.get(tier, dict()).get(aligner, None)

        if align_params is None:
            raise RuntimeError(f'No default alignment parameters for tier {tier} and aligner {aligner}')

    return align_params
