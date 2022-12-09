# PAV Example

This is a small example to help you get started running PAV. Files are included in `files/example` in the PAV
repository.

This assumes you have installed PAV or are running from a container (Docker or Singularity). Running from a container
is recommended if Docker or Singularity is available.

Native installs are recommended if individual PAV steps should be distributed for performance reasons (e.g. running a
single sample across multiple cluster nodes) or if Docker and Singularity are not available. See native install
instructions in `INSTALL_NATIVE.md`.

This example will assume you are in a clean directory where run configuration files are found and where results will be
written (the ANALYSIS directory).

You may need up to 32 GB of memory to map the alignments to the human reference.

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

This guide assumes `assemblies` is in the ANALYSIS directory for simplicity, although assemblies are often in a
different location in real runs. Adjust paths in `samples.tsv` if it is placed somewhere else.

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
not in the ANALYSIS directory (e.g. a shared location on the analysis machine). The reference path in `config.json` can
be updated to point to an existing reference (see the next step).

## Write configuration

Copy the `config.json` and `samples.tsv` configuration files from `files/example` in the PAV repository to the ANALYSIS
directory (e.g. `ANALYSIS/config.json`, do not place in a subdirectory). If a different reference is being used,
change the `reference`.

## Run

PAV may be run in several ways, through built-in run-scripts or through a container (Docker or Singularity). For those
more familiar with Snakemake, directly calling snakemake also works (not covered here).

### Run container (recommended)

From the ANALYSIS directory, run the container.

Docker:
```
sudo docker run --rm -v ${PWD}:${PWD} --user "$(id -u):$(id -g)" --workdir ${PWD} paudano/pav:2.1.1 -j 4 --nt
```

Singularity:
```
singularity run --bind "$(pwd):$(pwd)" --writable-tmpfs library://paudano/pav/pav:latest -j 4
```

You may need to adjust the directory bindings for your machine, but these parameters should work for most.


### Run scripts (native install)

There are two run scripts, `rundist` (distribute over a cluster, requires additional configuration) and `runlocal`
(run on a single machine). To use them, change to the ANALYSIS directory (location of `config.json` and `samples.tsv`)
and link the script there:

```
ln -s PATH/TO/PAV/runlocal
```

The script will know the current directory is ANALYSIS directory (where results will be written) and the real location
of the run script is in the SITE directory (where all the PAV code and `Snakefile` are located). This allows you to
maintain multiple versions of PAV, and the link serves as a record pointing to the version used for a specific project.

The first argument to the run script is the number of concurrent jobs to run. No other arguments are required.

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
