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

### Assembly Input

PAV can find assemblies in one of two ways:
* `assemblies.tsv`: Use for a small number of samples. Most will use this option.
* `asm_pattern` in `config.json`: If paths to the assemblies are consistent, then this method can be used to
  find the input without creating `assemblies.tsv` with an entry to each. This is useful for processing a large
  collection of assemblies (i.e. consortia-sized data).

Assemblies may be in FASTA, FASTQ, or GFA (hifiasm compatible) formats, and files may be optionally gzipped. File
indexes such as ".fai" files are not needed, PAV will create the indices it needs. It can also take an FOFN
(File Of File Names) pointing to multiple input files and processing them as one.

For phased assemblies, PAV expects two haplotypes separated into two files. It can process unphased assemblies by
leaving the haplotype column blank in `assemblies.tsv` or giving it a path to a 0-byte file (works for both
`assemblies.tsv` and `asm_pattern` input).

#### Finding assemblies: assemblies.tsv

A tab-separated-values (TSV) table is created with one line per assembly:

Create an assemblies TSV file, `assemblies.tsv`, with three columns:
1. NAME: Assembly name
1. HAP1: Path to haplotype 1 FASTA
1. HAP2: Path to haplotype 2 FASTA

The configuration option `assembly_table` (in `config.json`) may be used to set the name of the assemblies TSV file,
and the file may be gzipped.

Most spreadsheet applications including MS Excel and LibreOffice Calc can be used to edit this file. Be sure to save
as the right file type, it must be plain text and end with ".tsv" (not ".tsv.txt" unless the `assembly_table` option
is adjust to look for it).

If you have one assembly per sample (unphased), then leave HAP2 blank for those samples or give it a path to a 0-byte
file.

#### asm_pattern in config.json

Use this method if you have multiple assemblies in a consistent path structure. With this method, it is possible to
process hundreds of assemblies without creating an entry in `assembly_table.tsv` for each one.

Two wildcards will need to be inserted into the path, "asm_name" and "hap". "asm_name" can be any name for the assembly
or sample. "hap" must be "h1" and "h2". The input path goes into `config.json` parameter `asm_pattern`.

For example:

```
{
  "reference": "/path/to/hg38.no_alt.fa.gz",
  "asm_pattern": "/path/to/assemblies/{asm_name}/{hap}.fa.gz"
}
```

In this example, if an assembly with `asm_name` "HG00733" is run, then PAV will expect to find two files:

1. /path/to/assemblies/HG00733_CCS_SS_PG_PRR/h1.fa.gz
1. /path/to/assemblies/HG00733_CCS_SS_PG_PRR/h2.fa.gz

If there is no "hap" or "parent" (see below) wildcard in the path, it is treated as an unphased assembly by PAV and
PAV inputs the assembly into `h1` and never tries to read `h2`.

##### hifiasm-trio

To support hifiasm-trio, you may substitute the "parent" wildcard for "hap" (do not include both "hap" and "parent"
wildcards). In this case, "mat" becomes and alias for "h1", and "pat" becomes an alias for "h2" when searching for
files. The PAV output will still contain "h1" and "h2" for the maternal and paternal haplotypes, respectively.

##### sample wildcard

Optionally, "sample" may be a wildcard in the path. When PAV sees this, it assumes the sample name is the first part of
"asm_name" delimited by underscores. For example, if "asm_name" is "HG00733_CCS_SS_PG_PRR", then "sample" is "HG00733".
This was a feature used mainly for HGSVC where we used a "sample/assembly" directory structure (e.g.
`assemblies/HG00733/HG00733_CCS_SS_PG_PRR_h1.fa.gz` from pattern `assemblies/{sample}/{asm_name}_{hap}.fa.gz`).
This may be useful for consorita with several assemblies per sample.

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
singularity run --bind "$(pwd):$(pwd)" --writable-tmpfs library://becklab/pav/pav:latest -c 16
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
