# PAV Example

This is a small example to help you get started running PAV. Files are included in `files/example` in the PAV
repository.

By convention, PAV requires two directories:

* SITE: Location where PAV code is installed (directory produced by `git clone ...`). Only required if PAV is run
directly (the SITE directory is built into the PAV Docker and Singularity containers).
  * `Snakefile` in the root of the PAV repository will be in this location. 
* ANALYSIS: Location where analysis is carried out.
  * `config.json` will be in this location.
  * Output and temporary files are written to this location.

The PAV run scripts will link the SITE and ANALYSIS directories allowing different versions of PAV to be used for
different projects.

## Install PAV

If you are running PAV outside of a container (Docker or Singularity), then it should be installed into its own
directory (the SITE directory).

Clone PAV:
```
git clone --recursive https://github.com/EichlerLab/pav.git
```

The `pav` directory produced by `git clone` is the SITE directory (contains `Snakefile`).

### Runtime envorinment

At runtime, PAV needs an environment containing Python3, Python libraries, and some command-line programs
(see `README.md`). These can be installed in several ways, such as through a conda environment (activate the environment
before running PAV) or through command-line tools with PATH set accordingly (such as Linux enivornment modules). 

## Download assembly data

This example uses phased assemblies for HG00733 from Ebert 2021 (see citation in README.md) covering
GRCh38 region 22q12 (11.7 Mbp of chr22).

The phased assembly is in two FASTA files, one for each phased haplotype:

* Haplotype 1 (h1): HG00733_22q12_h1.fa.gz (4.1 MB)
* Haplotype 2 (h2): HG00733_22q12_h2.fa.gz (4.0 MB)

Download these from the EBI IGSR FTP:

```
mkdir assemblies

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20221202_PAV_Example/* -P assemblies/
```

This guide assumes `assemblies` is in the ANALYSIS directory for simplicity. Adjust paths in `samples.tsv` if it is
placed somewhere else.

## Download the reference

In this example, `config.json` is configured to use a non-ALT version of the reference, which is the primary assembly.
We do not recommend using an assembly with ALTs, decoys, or patches because they may produce incorrect results with
long read aligners.

Download the reference:

```
mkdir hg38_no_alt

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/* -P hg38_no_alt/
```

This guide assumes `hg38_no_alt` is in the ANALYSIS directory for simplicity. In real analyses, references are usually
not in the ANALYSIS directory (e.g. a shared location on the analysis machine).

### Use another reference

If a different reference is already downloaded, adjust `config.json` (next step) to point to the reference FASTA.
Relative paths are from the PAV analysis directory, absolute paths are relative to the root of the filesystem
(i.e. `/`). The reference FASTA may be uncompressed or gzipped and does not need to be indexed (i.e. the `.fai` file is
not needed).

## Write configuration

Copy the `config.json` and `samples.tsv` configuration files from `files/example` in the PAV repository to the ANALYSIS
directory (e.g. `ANALYSIS/config.json`, do not place in a subdirectory). If a different reference is being used,
change the `reference`.

## Run

PAV may be run in several ways, through built-in run-scripts or through a container (Docker or Singularity). For those
more familiar with Snakemake, directly calling snakemake also works (not covered here).

### Run scripts

There are two run scripts, `rundist` (distribute over a cluster) and `runlocal` (distribute locally). To use them,
change to the ANALYSIS directory (location of `config.json` and `samples.tsv`) and link the script there:

```
ln -s PATH/TO/PAV/runlocal
ln -s PATH/TO/PAV/rundist
```

The script will know the current directory is ANALYSIS directory (where results will be written) and the real location
of the run script is in the SITE directory (where all the PAV code and `Snakefile` are located). This allows you to
maintain multiple versions of PAV, and the link serves as a record pointing to the version used for a specific project.

The first argument to the run script is the number of concurrent jobs to run. No other arguments are requried.

Run PAV with 8 concurrent jobs:

```
./runlocal 8
```

Any additional arguments are given directly to Snakemake and can control Snakemake parameters or request specific
output files.

Run HG00733 and do not remove temporary files (`--nt`) in `temp` directory:

```
./runlocal 8 --nt pav_HG00733.vcf.gz
```

The syntax for `rundist` is the same as `runlocal`, but additional configuration is needed to distribute it over your
cluster (cluster dependent, see `README.md`).
