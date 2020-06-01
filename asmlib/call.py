"""
Variant caller functions.
"""


def get_gt(row, hap, map_tree):
    """
    Get variant genotype based on haplotype and mappability.

    :param row: Variant call row (post h1/h2 merge).
    :param hap: Haplotype being called.
    :param map_tree: Tree of mappable regions.

    :return: '1' if the variant call is in this haplotype, '0' if it was not called but is in mappable regions,
        and '.' if it is not in a mappable region.
    """

    if hap in row['HAP_SRC'].split(','):
        return '1'

    for interval in map_tree[row['#CHROM']].overlap(row['POS'], row['END']):

        if (interval.begin <= row['POS']) and (interval.end >= row['END']):
            return '0'

    return '.'
