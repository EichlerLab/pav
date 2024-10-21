# PAV changes

## 2.4.5
* Corrected bug in VCF output for inversions and VCF output format for other variant types. 

## 2.4.4
* Bad container build, re-released as 2.4.5

## 2.4.3
* Added "no_link_qry" configuration parameter to force PAV to re-write query (assembly) FASTA files if this parameter
  is set to True. By default (no_link_qry = False), an acceptable input FASTA file is linked (single compressed FASTA
  as input).

## 2.4.2
* Bad container build, re-released as 2.4.3

## 2.4.1.2
* VCF writes failed variants with a FILTER
  * INFO/COMPOUND lists larger variants that caused smaller ones underneath it to be removed.
* Added "ignore_env_file" configuration parameter to be used by the Docker entrypoint script to avoid disrputing the container environment if a "setenv.sh" file is present (container should use it's own environment).

## 2.4.1
* Final VCFs write to vcf/merged
* Avoid FASTA copy to temp (symlink)

## 2.4.0.2 (development version)
* Filter variants inside inverted core (uniquely inverted loci, inner INV coordinates). Previous versions filtered all variants inside inversions (outer INV coordinates)
* Skip input FASTA re-indexing if indexes already exist

## 2.4.0.1 (development version)
* Arbitrary number of named haplotypes and new assembly table format (old format supported).
* No support for "asm_pattern" in config (was removed before this version but was left in documentation).
* Removed "chrom_cluster" from configuration options.
* Split alignment trimming into distinct "tig" and "ref" phases (BED for each).
* Keeping dropped variants and adding FILTER column with PASS and non-PASS values. Added COMPOUND field (and filter name
  called COMPOUND) to annotate variants that are dropped because they are inside another variant (i.e. small variant in
  a large DEL, COMPOUND field points to the accepted variant covering the dropped one). Makes it easier to recover
  variants if a large false-SV is removed.
* Added ".fna" and ".fnq" to recognized file types.
* Moving to Snakemake profile. Removed threading parameters from config (map_threads, inv_threads_lg, inv_threads)
  Configure in the profile instead of parameters.
* Added joblib to container dependencies, removed Snakemake version restriction
* VCF default rule updates

## 2.3.4
* Adding "COMPOUND" column for filtered variants annotating which variants it was nested inside of (e.g. small variant
  inside a large deletion). Can be used to recover small variants if a larger variant is removed from the callset.
* Fixed bug that dropped large insertions.
* Fixed a bug generating tracks for references with numeric chromosomes.

## 2.3.3
* Fixed minor alignment trimming bug uncovered by the large deletion fix in 2.3.2 causes large variant detection fail.
* Fixed bug in coordinates causing inversion searches at the edge of chromosomes to fail with an invalid index (-1).

## 2.3.2
* Fixed large deletion bug.
* Removed "fake home" from container run script, causes problems with some systems. Some tools in the container will use real home (e.g. .cache)

## 2.3.1
* Fixed typo in Singularity run script from home directory caching issues (thanks @aidenlx)
* Fixed bug in variant filtering
* Updated documentation for native installs related to deadlocks (OPENBLAS_NUM_THREADS)

## 2.3.0
* Fixed home directory caching issues.

## 2.2.9
* Fixed issues with home directory caching (related to bug fixes in 2.2.7)

## 2.2.8
* Fixed rules for generating UCSC tracks for PAV calls and alignments. 

## 2.2.7
* Fixed bug creating the callable region BED if large SV BED files are empty.
* Eliminating the need for --writable-tmpfs (Singularity), fails for some kernels

## 2.2.6
* SNV ID format is now: CHROM-POS-END-REFALT
  * Eliminated dash between REF and ALT (single base each), consistent with SV-Pop library and makes IDs more paresable.
* Updated SV-Pop submodule.

## 2.2.5
* Added support for systems without /dev/shm
  * Set "inv_threads" and "inv_threads_lg" to "1" in config.json. Slower, but will not try to multiprocess.
* Fixed bug in variant calling ("Passing a set as an indexer is not supported.", #34)

## 2.2.4
* Updating for LRA index file name change from ".mmi" to ".mms".

## 2.2.3
* Minor SV-Pop bug incorrectly de-duplicates variant IDs.

## 2.2.2
* Minor SV-Pop bug incorrectly parsed chromosome names with "." while versioning variants.

## 2.2.1
* Fixed missing library path in scripts/density.py
* Quieted merging output
* Minor documentation changes

## 2.2.0
* Dockerized PAV
* Added a small example (EXAMPLE.md)
* Fixed bug writing VCFs with no variants of a certain type (most often caused by no inversions).
* Fixed bugs when reference has numeric chromosomes (e.g. "1" instead of "chr1")

## 2.1.0
* Added parameter "inv_inner" (boolean, default False) to allow small variants inside inversions
* Added tunable parameter "minimap2_params" (thanks to Edwmard Rice).
* If contig name is a number, convert it to str (thanks to Peng Jia).
* Inversion caller bug fixes:
  * Incorrectly calling large inversions based on misaligned contigs.
  * Small inversions were being missed (scripts/density.py wasn't able to find kanapy, also fixed reason it was failing silently).
* Batching haplotype merges (20 batches) instead of one per chrom (generates too many jobs on unplaced/unlocalized contigs)

## 2.0.1
* If hap and parent are missing in asm_pattern, it's an unphased assembly (no longer an error).
* Fixed INV signature bug when there are no indel clusters flagging INVs.
* Reading small indels that cross Ns (particularly "NA") causes Pandas to read the table back with missing data instead of a "NA" deletion or insertion
* Merging info changes (removed HAP_AC and HAP_AF, not meaningful).
* Updated parameters for redundant callsets (multiple haplotypes in one assembly).
* Added missing code for creating dotplots.
* De-duplicating IDs before merging (for redundant callsets).
* Fixed sys paths in scripts/reconstruct_sam.py.

## 2.0.0
* Haplotype merging uses sequence content by default (not just size/position).
* Fixed tuple assignment bug while finding input files (#12, thanks to Edward Rice).
* Fixed SNV merge strategy typo (#14).

## 1.2.2
* More granular merging control based on variant type.
* VCF SVLEN header had the incorrect data type.

## 1.2.1
* Changed parameter "merge_align" to "merge_match" to better match matching paradigm.

## 1.2.0
* PAV merges using sequence context by default. Sequences are only merged if
  they are 80% identical by alignment parameters. Merging is done with SV-Pop
  parameters "match=0.8,2,-1,-4,-0.25,500000,9" (80% identical, match=2, mismatch=-1,
  gap-open=-4, gap-extend=-0.25, align up to 500,000 bp (Jaccard similarity if larger),
  Jaccard k-mer size 9). Parameter is tunable with option "merge_align", which may be
  "true" (default), "false" (no sequence matching), or a string that replaces everything
  after "match=".

