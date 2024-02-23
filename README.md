<h1 align="center"><img width="300px" src="img/logo/PAVLogo_Full.png"/></h1>
<p align="center">Phased Assembly Variant Caller</p>

***
<!-- Templated header from the pbsv github page: https://github.com/PacificBiosciences/pbsv -->

PAV is a tool for discovering variation using assembled genomes aligned to a reference. It supports both phased and
unphased assemblies.

PAV was developed for the Human Genome Structural Variation Consortium (HGSVC)

    Ebert et al., “Haplotype-Resolved Diverse Human Genomes and Integrated Analysis of Structural Variation”, Science, February 25, 2021, eabf7117, https://doi.org/10.1126/science.abf7117.


PAV was originally developed as part of the Eichler lab at UW and is now updated and maintained by the Beck lab at JAX.
Both labs continue to contribute to the HGSVC.

Eichler lab:
https://eichlerlab.gs.washington.edu/

Beck lab:
https://www.jax.org/research-and-faculty/research-labs/the-beck-lab

## Configuring PAV

Change to a clean directory (the ANALYSIS directory) to run PAV. PAV will read `config.json` from this directory and
write output to this directory. If you have a native install, do not run PAV from the PAV install location (the SITE
directory where `Snakefile` and `pavlib` are found).

PAV gets it's configuration from two files:

* `config.json`: Points to the reference genome and can be used to set optional parameters.
* `assemblies.tsv`: A table of input assemblies.

### Base config: config.json

A JSON configuration file, `config.json`, configures PAV. Default options are built-in, and the only required option is
`reference` pointing to a reference FASTA file variants are called against.

Example:

```
{
  "reference": "/path/to/hg38.no_alt.fa.gz"
}
```

Note: The HGSVC reference for long reads can be found here:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/

A no-ALT version of a reference is essential for long read assemblies. A reference containing alternate loci, decoys, or
patches may produce unexpected results and will likely lead to a loss of sensitivity.

#### Haplotype assembly input

Assemblies may be in FASTA, FASTQ, or GFA (hifiasm compatible) formats, and files may be optionally gzipped. File
indexes such as ".fai" files are not needed, PAV will create the indices it needs. It can also take an FOFN
(File Of File Names) pointing to multiple input files and processing them as one.

The assembly table has one line per sample. The `NAME` column contains the assembly name (or sample name). This column
must be present and must not contain duplicates. All sample names must be composed of only alpha-numeric characters,
underscores, and dashes (A-Z, a-z, 0-9, "_", "-").

PAV accepts an arbitrary number of haplotypes per sample. Each haplotype is a column with one of two name formats:
* HAP_* (where * is a haplotype name): A named haplotype. For example, "HAP_h1", "HAP_h2" would be typical columns
  for a phased assembly. Other pseudo-haplotypes may be added, such as "HAP_unphased". Haplotypes must only contain
  alpha-numeric and "-", "+", "." characters. Underscores are not allowed to avoid ambiguity with underscore-separated
  assembly name and haplotypes (e.g. "a_b_c" could be sample "a_b" with haplotype "c", or sample "a" with haplotype
info_header_list.append(('HAP_VARIANTS', '.', 'STR', 'List of variant IDs in the order of "HAP" identifiying the variant merged in for each haplotype'))info_header_list.append(('XXX', '.', 'STR', 'XXX (INFO/HAP order)'))  "b_c").
* HAP* (where * is an integer): Legacy mode supporting tables with "HAP1" and "HAP2" column names. These are translated
  to "HAP_h*" internally.

Each entry in a "HAP_" column is the name of a file (FASTA, FASTQ, GFA, or FOFN). Multiple files can be input by
separating them by semi-colons (i.e. "path/to/file1.fasta;path/to/file2.fasta").

PAV will attempt to replace "{asm_name}" and "{hap}" wildcards in paths with the assembly name (NAME column) and
the haplotype name making it easier to generate paths in the assembly table if input files follow the same file naming
conventions (can use the same path pattern for many samples and haplotype columns). For example, with sample "HG00733"
and haplotype "h2", path pattern "/path/to/assemblies/{asm_name}/{asm_name}_{hap}.fa.gz" becomes
"/path/to/assemblies/HG00733/HG00733_h1.fa.gz". Another wildcard "{sample}" is treated as an alias of "asm_name"
(i.e. the example above could have used either the "{asm_name}" or "{sample}" wildcards to achieve the same result).

PAV will run for all haplotypes with a non-empty "HAP_" column. The variant calling process is done independently for
each haplotype, and variant calls are merged at the end in the order they are found in this table. To include a
haplotype with no input, add an entry to a zero-byte file to make the haplotype appear in the merged variant table and
in VCF genotypes.

#### Assembly-specific configuration

Global configuration parameters (those found in `config.json`) can be set per-assembly, which wil override `config.json`
for specific assemblies. The configuration string can be placed into the optional "CONFIG" column and is a
semicolon-separated list of key-value pairs (i.e. "key1=val1;key2=val2").


#### Assembly filter

A BED file may be input per assembly marking unreliable assembly regions. These regions might be curated using results
from assembly QC tools (such as Flagger or NucFreq). Any variants intersecting these assembly loci (even by 1 bp) will
be marked as filtered with reason "TIG_FILTER". Input BED files may be plain text or gzipped and must be in assembly
coordinates (i.e. #CHROM, POS, and END correspond to locations on the assembly, not the reference).

The configuration table optional "FILTER_" columns point to these files with the haplotype column name following.

Examples:
* If the haplotype column is "HAP_h1", then the corresponding filter column is "FILTER_h1".
* If the haplotype column is "HAP1", then the corresponding filter column is "FILTER_HAP1".

If these columns are absent or empty, no filter is applied.

Like assembly input columns, PAV will attempt to replace "{asm_name}" and "{hap}" in paths with the sample and haplotype
name (also allowing "{sample}" as an alias for "{asm_name}").


### Additional configuration parameters

Additional information about configuration parameters for `config.json` can be found in `CONFIG.md`.

## Running PAV from Docker and Singularity
 
Change to the ANALYSIS directory (where `config.json` is found), then run the container:

Docker:
```
sudo docker run --rm -v ${PWD}:${PWD} --user "$(id -u):$(id -g)" --workdir ${PWD} becklab/pav:latest -c 16
```

Singularity:
```
singularity run --bind "$(pwd):$(pwd)" library://becklab/pav/pav:latest -c 16
```

Notes:
1. **Cores**: Set the maximum number of cores `-c` (or `--cores`) to be used simultaneously.
1. **Directory binding**: You may need to adjust the directory bindings for your machine, but these parameters should work for most.
1. **Version**: You may change "latest" to an explicit PAV version to ensure compatibility among samples.

PAV can process a phased human genome in 4.5 to 5.5 hours with 64 GB of memory and 32 cores with minimap2 alignments.
Actual memory usage is around 52 GB.

## Running native PAV

See `NATIVE_INSTALL.md` for help installing and PAV natively on a machine. This option necessary if Docker and
Singularity are not available or if distribute individual PAV steps over a cluster.

## Running a small example

See `EXAMPLE.md` to setup small example run to test PAV on your system.


## Interpreting output

Most projects will read from the VCF in the root of the run directory, but PAV outputs some other useful information.

The output directory (`results/{asm_name}`) has several subdirectories:

1. align: Information about contig alignments.
    1. Post-trimming alignments.
    1. The BED and FASTA files in this directory could be used to reconstruct a SAM file.
1. bed: Variant calls in formatted BED files
    1. One file for each variant type (sv_ins, sv_del, indel_ins, indel_del, snv_snv)
1. bed/fa: FASTA for inserted and deleted sequences
    1. Unique ID links sequnece to variant call
    1. No FASTA for SNVs (see REF and ALT in variant calls)
1. callable: BED files of callable regions (where contigs aligned) smoothed by 500 bp windows.
1. inv_caller: Intermediate output from the inversion caller.
    1. Flagged loci queried for inversions (not all produce calls).
    1. Contains data useful for visualizing inversions.
1. lg_sv: Intermediate output from large SV calls.

## Haplotype merging

Information about how PAV resolves two haplotypes as one diploid sample can be found in `HAP_MERGING.md`.

## Cite PAV

Ebert et al., “Haplotype-Resolved Diverse Human Genomes and Integrated Analysis of Structural Variation,”
Science, February 25, 2021, eabf7117, https://doi.org/10.1126/science.abf7117 (PMID: 33632895).

PAV was also presented at ASHG 2021:

Audano et al., "PAV: An assembly-based approach for discovering structural variants, indels, and point mutations
in long-read phased genomes," ASHG Annual Meeting, October 20, 2021 (10:45 - 11:00 AM), PrgmNr 1160

## Contact

Please open a case on the Github page for problems.

You may also contact Peter Audano directly (e-mail omitted to reduce SPAM). PAV was developed in the lab of
Dr. Evan Eichler at the University of Washington and is currently maintained in the lab of Dr. Christine beck at The
Jackson Laboratory.
