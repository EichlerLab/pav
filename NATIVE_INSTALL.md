# Native PAV

Follow these instructions to install and run PAV natively on a host. These instructions do not apply to Docker
or Singularity, which encapsulates the PAV deployment inside a container.

A native deployment is useful if you want more control over PAV, such as distributing it over a cluster, which
requires additional cluster configuration. If PAV is run on a single host, containerized PAV (Docker or Singularity) is
faster to setup and easier to run.

## Directories

This documentation will refer to two directories:

* SITE: Location where PAV code is installed (directory produced by `git clone ...`). Only required if PAV is run
directly (the SITE directory is built into the PAV Docker and Singularity containers).
  * `Snakefile` in the root of the PAV repository will be in this location. 
* ANALYSIS: Location where analysis is carried out.
  * `config.json` will be in this location.
  * Output and temporary files are written to this location.

The PAV run scripts will link the SITE and ANALYSIS directories allowing different versions of PAV to be used for
different projects.

## Install PAV

### Pull PAV code

Go to the directory where the SITE directory will be (SITE will be `pwd`/pav) and pull PAV:

```
git clone --recursive https://github.com/EichlerLab/pav.git
```

There is nothing to compile, PAV and submodule libraries are pure Python.

### Dependencies 

PAV requires Python 3 with the following Python libraries installed:
1. BioPython
1. intervaltree
1. matplotlib
1. numpy
1. pandas
1. pysam
1. scipy

Command line tools needed:
1. minimap2 (default aligner)
1. lra (optional alternate aligner)
1. samtools
1. bedToBigBed (optional, from UCSC browser tools)

Optional tools: The aligner defaults to minimap2, but PAV also supports LRA. PAV only requires the aligner you use,
both do not need to be installed. PAV can create browser tracks for variant calls, and if this feature is used, the UCSC
browser tools are needed (bedToBigBed). Pre-compiled binaries can be obtained from
http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

### Setting up the run directory

Link `Snakefile` from the SITE directory to the ANALYSIS directory.


### Run scripts

To execute via run scripts, go to the run directory and link `rundist` and `runlocal` from the PAV install directory

Example:

    ln -s /path/to/PAV/1.0.0/rundist ./
    ln -s /path/to/PAV/1.0.0/runlocal ./

`rundist` will be used to distribute over a cluster, and `runlocal` will be used to run the pipeline in the current
session. Configuring `rundist` will require some knowledge for distributing over a cluster.

Both scripts setup some control variables and pass control over to a user-defined script. `rundist` will search for
`config/rundist.sh`, and `runlocal` will search for `config/runlocal.sh`. It first searches the run directory, then
it searches the PAV install directory. Once found, it calls that script to carry out calling Snakemake. If `runlocal.sh`
is not found, PAV will use a built-in default (`files/run_scripts/runlocal.sh`), but the environment must be setup
with all PAV's dependencies (Snakemake, python with libraries, minimap2, samtools, etc -- See "Dependencies" in this
README).

Generally, you would add your run scripts to the PAV install directory so it could be run for any number of projects.
Some runs may need custom resources not typical for other projects, and for those, you can override the script in the
PAV install directory with one in the run directory.

Examples for what your `config/rundist.sh` and `config/runlocal.sh` might look like are in comments at the bottom of
`rundist` and `runlocal`.

After the symbolic link is created and your `config/rundist.sh` and/or `config/runlocal.sh` are in place, PAV is run:

    ./rundist 20 pav_HG00733_CCS_SS_PG_PRR.vcf.gz

or

    ./runlocal 20 pav_HG00733_CCS_SS_PG_PRR.vcf.gz

The first number (20 in the example) is the number of concurrent jobs. Everything else is passed to Snakemake as a
target, and there may be multiple targets. If there are no targets, the default rule is run, which reads
`assemblies.tsv` and generates the VCF for all. If you are using wildcards in `config.json` instead of `assemblies.tsv`,
the targets are required.

### Snakemake direct

Examples in this section will assume shell variable `PAV` is set to the PAV install directory (directory with
`Snakemake` in it).

If PAV was configured to read from `assemblies.tsv`, then the default rule will run all assemblies.

To run a single sample, request any output file from Snakemake.

For example:

    snakemake -s ${PAV}/Snakefile pav_HG00733_CCS_SS_PG_PRR.vcf.gz

Where "HG00733_CCS_SS_PG_PRR" is the name of the sample or assembly to be run.

## Additianal errors and information

### Deadlock and OPENBLAS_NUM_THREADS

If PAV hangs for no apparent reason, try running with "OPENBLAS_NUM_THREADS=1" set in the environment. This could
be done by adjusting your runlocal.sh or rundist.sh or exporting it before running PAV. There is a reasonable
explanation here:
https://stackoverflow.com/questions/56104472/why-would-setting-export-openblas-num-threads-1-impair-the-performance

This issue was observed in GitHub issue #39. If you see this behavior in the Docker or Singularity container distributed
with PAV, please report a bug.
