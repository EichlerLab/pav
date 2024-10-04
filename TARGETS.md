# PAV targets


## Init

* data_init: Generate all files that are needed to run PAV across samples. This rule is intended to be used if base
  tables and resources should be setup before parallelizing over individual samples in separate PAV runs. This rule
  is not needed if PAV is run in one command for all samples, and in this case, PAV will run the rules to generate these
  resources as needed. If PAV runs are split by sample, then executing this rule will prevent conflicts when processing
  individual samples in their own PAV runs.

## Alignments

* align_all: Run alignments and stop
* align_export_all: Export CRAM, BAM, or SAM (.sam.gz)
  * --config options: trim, export_fmt, asm_name, hap
    * Each is a comma-separated list of values. Defaults to all trim levels (trim) and CRAM (export_fmt),
      and all assemblies and haplotypes based on the configuration (asm_name and hap)
    * export_fmt may be "cram", "bam", or "sam". "sam.gz" is an alias for "sam".
  * See next target for details 
* Export a single alignment: "results/{asm_name}/align/export/pav_align_trim-{trim}_{hap}.{ext}"
  * ext may be "cram", "bam", or "sam.gz" ("sam" is not recognized)
  * All files are indexed (SAM is bgzipped and tabix indexed)
