# Pipeline definitions

# Variant call filters.
FILTER_REASON = {
    'PASS': 'Variant passed filters',
    'TIG_FILTER': 'Contig filter region (tig_filter_pattern)',
    'COMPOUND': 'Inside larger variant',
    'COMPOUND_INV': 'Inside a larger inversion',
    'INV_MIN': 'Less than min INV size ({})',
    'INV_MAX': 'Exceeds max INV size ({})',
    'TRIM_REC': 'Alignment trimming removed whole alignment record',
    'TRIM_FLANK': 'Alignment trimming removed variant region'
}

# Variant call merge batches
MERGE_BATCH_COUNT = 20
