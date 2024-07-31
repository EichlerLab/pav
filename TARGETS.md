# PAV targets

## Alignments

* align_all: Run alignments and stop
* align_export_all: Export CRAM, BAM, or SAM (.sam.gz)
  * --config options: trim, tier, export_fmt, asm_name, hap
    * Each is a comma-separated list of values. Defaults to all trim levels (trim), tiers (tier), CRAM (export_fmt),
      and all assemblies and haplotypes based on the configuration (asm_name and hap)
    * export_fmt may be "cram", "bam", or "sam". "sam.gz" is an alias for "sam".
  * See next target for details 
* Export alignment: "results/{asm_name}/align/export/pav_align_tier-{tier}_trim-{trim}_{hap}.{ext}"
  * ext may be "cram", "bam", or "sam.gz" ("sam" is not recognized)
  * All files are indexed (SAM is bgzipped and tabix indexed)

