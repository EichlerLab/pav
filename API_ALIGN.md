# Alignment API and data structures

## Nomenclature

PAV refers to alignments using "query" and "subject". A query is the sequence being aligned (a contig) and the subject
is the sequence it is aligned to (the reference).

## Coordinate systems

Two coordinate systems appear in PAV, BED coordinates and 1-based coordinates.

### BED coordinates

Start (POS) and end (END) positions are in a base-0 system where the first base in the reference or query sequence is
0, the second base is 1, and so on. The start positions give the location of the base where an alignment or event
starts, and the end position is one beyond. This is often referred to as 0-based-half-open coordinates. Computations
with this system are easier to write and less prone to bugs.

Wherever coordinates appear in separate fields, such as POS and END columns in a table, they are this coordinate system.

### 1-based coordinates

Start and end positions point to the base where an event starts and ends with the first base in the sequence is at
position 1. This system is used by UCSC, IGV, and others. It is more intutivie as a string, but harder to handle
programatically.

Whenever coordinates appear as a single string (e.g. "chr1:12345678-12346678"), it is in 1-based coordinates. This
rarely appears in PAV, but some library routines output it.


## Alignment table

The alignment tables in "results/{asm_name}/align" in BED format have these fields:
* #CHROM: Subject (reference) name.
* POS: Aligned subject (reference) start position in BED coordinates.
* END: Aligned subject (reference) end position in BED coordinates.
* INDEX: Index of the alignment record in the original SAM output where the first alignment record is 0. Ordering is
  preserved for alignment records that were removed (i.e. if record 2 was removed, there will be no INDEX with 2 and 3
  is still the index of the record that came after 2). Gives a unique ID to each alignment record.
* QRY_ID: Query name.
* QRY_POS: Query (contig) start position in BED coordinates in the original assembly coordinates (see note below).
* QRY_END: Query (contig) end position in BED coordinates in the original assembly coordinates (see note below).
* QRY_LEN: Length of the full contig (not the aligned part of the contig).
* RG: RG tag from the alignment.
* AO: AO tag from LRA alignments. From LRA documentation (https://github.com/ChaissonLab/LRA), "This number shows the
  order of the aligned segment when a read is split".
* MAPQ: Mapping quality (PHRED scale).
* REV: True if alignment was reverse-complemented during alignment.
* FLAGS: SAM flags in hexadecimal format.
* HAP: Haplotype.
* CIGAR: CIGAR string.
* CALL_BATCH: Batch number for variant calling steps. Alignments are separated by contig IDs into batches of similar
  sizes and each batch can be run in a separate job.
* TRIM_REF_L: Number of subject (reference) bases that were trimmed from the left (upstream) side of alignment during alignment trimming.
* TRIM_REF_R: Number of subject (reference) bases that were trimmed from the right side (downstream) of alignment during alignment trimming.
* TRIM_QRY_L: Number of query (contig) bases that were trimmed from the left side (upstream) of alignment during alignment trimming. Left and right are relative to the aligned reference orientation.
* TRIM_QRY_R: Number of query (contig) bases that were trimmed from the right side (downstream) of alignment during alignment trimming. Left and right are relative to the aligned reference orientation.

TRIM fields are added after the initial stages and are not present in the BED before trimming.

QRY_POS and QRY_END are the aligned positions of the contig in the *original* assembly coordinates regardless of the
contig's orientation with the reference. To extract the sequence for a record, use these coordinates (e.g.
"QRY_ID:(QRY_POS + 1)-QRY_END" in samtools faidx), and to place it in reference orientation, reverse complement the
extracted sequence if REV is True.

