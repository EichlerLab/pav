# Phased Assembly Variant Caller (PAV)

PAV is a tool for discovering variation using assembled genomes aligned to a reference. It is designed explicitly for
phased assemblies, however, it can be used for squashed assemblies by providing an empty FASTA for the second haplotype.

PAV was developed for the Human Genome Structural Variation Consortium (HGSVC) (Pending publication).

## Installing PAV

### Dependencies

Choose a clean directory to install PAV ("install directory"). This will be a clean directory, and analyses will run in
separate locations.

PAV requires Python 3 with the following Python libraries installed:
1. BioPython
1. intervaltree
1. matplotlib
1. numpy
1. pandas
1. pysam
1. scipy

PAV will also need the contig aligner. This defaults to minimap2, and PAV also supports LRA (Chaisson lab, USC, not
yet unpublished).

Finally, Samtools should be accessible. PAV can create browser tracks of variant calls, and if this feature is used,
UCSC browser tools are needed (bedToBigBed).

### Additional libraries

There are two additional libraries PAV needs installed or linked into the `dep` directory.

    mkdir dep
    git clone https://github.com/EichlerLab/svpop
    ln -s svpop/analib
    git clone git@github.com:paudano/kanapy.git

The SV-Pop libraries contain code for handling variant calls including merging (haplotype callsets to diploid callset).
Kanapy is a small k-mer library used by the inversion caller.


## Configuring PAV

Once PAV is installed, move to a clean directory ("analysis directory") where PAV will be run. The following
configuration steps will be done in this analysis directory.


### Base config: config.json

A JSON configuration file, `config.json`, configures PAV. Default options are built-in, and the only required option is
"reference" pointing to a reference FASTA file variants are called against.

Example:

    {
      "reference": "/path/to/hg38.no_alt.fa.gz"
    }

Note: The HGSVC reference can be found in ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/

A no-ALT version of a reference is essential for PAV. A reference containing alternate loci, decoys, or patches may
produce unexpected results and will likely lead to a loss of sensitivity.

### Assembly Input

There are two ways to find assemblies, "asm_pattern" in `config.json`, or entries in `assemblies.tsv`. Choose one of
these methods.

#### asm_pattern in config.json

Use this method if you have multiple assemblies located in predictable paths.

Two wildcards will need to be inserted into the path, "asm_name" and "hap". "asm_name" can be any name for the assembly
or sample, and PAV results will appear in `results/asm_name` for each assembly processed. "hap" must be "h1" and "h2".
This path goes into `config.json` as "asm_pattern".

For example:

    {
      "reference": "/path/to/hg38.no_alt.fa.gz",
      "asm_pattern": "/path/to/assemblies/{asm_name}/{hap}.fa.gz"
    }

When assembly "HG00733_CCS_SS_PG_PRR" is run, then PAV will expect to find two files:

1. /path/to/assemblies/HG00733_CCS_SS_PG_PRR/h1.fa.gz
1. /path/to/assemblies/HG00733_CCS_SS_PG_PRR/h2.fa.gz

##### hifiasm-trio

To support hifiasm-trio, you may substitute the "hap" wildcard for "parent" (do not include both "hap" and "parent"
wildcards.). In this case, "mat" becomes and alias for "h1", and "pat" becomes an alias for "h2" when searching for
files. The PAV output will still contain "h1" and "h2" for the maternal and paternal haplotypes, respectively.

##### sample wildcard

Optionally, "sample" may be a wildcard in the path. When PAV sees this, it assumes the sample name is the first part of
"asm_name" delimited by underscores. For example, if "asm_name" is "HG00733_CCS_SS_PG_PRR", then "sample" is "HG00733".
This was a feature used mainly for HGSVC where we used a "sample/assembly" directory structure.

#### Finding assemblies: assemblies.tsv

Use this method if you have a small number of assemblies. A table entry will be created for each one.

Create an assemblies TSV file, `assemblies.tsv`. It is a tab-delimited file with three columns:
1. NAME: Assembly name
1. HAP1: Path to haplotype 1 FASTA
1. HAP2: Path to haplotype 2 FASTA

The configuration option "assembly_table" (`config.json`) may be used to set the name of the assemblies TSV file. 

## Running PAV

Examples in this section will assume shell variable `PAV` is set to the PAV install directory (directory with
`Snakemake` in it).

PAV currently writes variant calls as formatted BED files, but a VCF writer is expected to be implemented soon. If
PAV was configured to read from `assemblies.tsv`, then the default rule will run all assemblies.

To run a single sample, request any output file from Snakemake.

For example:

    snakemake -s ${PAV}/Snakefile results/HG00733_CCS_SS_PG_PRR/bed/sv_ins.bed.gz

## Interpreting output

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
1. lg_sv: Intermediate output from the lage

