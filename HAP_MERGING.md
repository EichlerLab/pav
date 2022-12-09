# Diploid callset merging

For phased assemblies, the final diploid callset is found by merging two independent callsets for one haplotype.

PAV uses SV-Pop's merging facility.

Default parameters are:
* Insertions and deletions (SV and indel): `nr::ro(0.5):szro(0.5,200,2):match`
* SNVs: `nrsnv::exact`

The format of this string is the merging engine ("nr" or "nrsnv"), then two colons, then merging stages separated by
single colons. An optional "match" at the end tells it to match variants by sequence content.

* ro: Reciprocal overlap (RO).
  1. Overlap proportion (0.5 is 50% RO).
* szro: Size-reciprocal-overlap.
  1. Size proportion (like RO if variants were right on top of each other)
  1. Maximum distance (minimum of start positon difference or end position difference)
  1. Maximum distance as a proportion of the variant size (e.g. 2 means offset may not be more than 2x variant size)

"nrsnv" enforces REF and ALT matches, and "exact" means exact position only. Offset distances are allowed (see SV-Pop)

The parameters are tunable using 'merge_TYPE' where TYPE is the type of the variant (acceptable values are ins, del,
inv, insdel, insdelinv, or snv).

The default parameters do a better job not merging multi-allelic sites into the same variant call, but do inflate the
number of variants in tandem repeats. The parameters previously used (early versions of PAV and HGSVC) can be
used by setting `merge_insdelinv` in `config.json` to "nr::ro(0.5):szro(0.5,200,2)" (no "match" directive).
