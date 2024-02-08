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
* #CHROM
* POS
* END
* INDEX
* QRY_ID
* QRY_POS
* QRY_END
* QRY_LEN: Length of the full contig (not the aligned part of the contig).
* RG
* AO
* MAPQ
* REV
* FLAGS
* HAP
* CIGAR

