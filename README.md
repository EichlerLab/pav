<h1 align="center"><img width="300px" src="img/logo/PAVLogo_Full.png"/></h1>
<p align="center">Phased Assembly Variant Caller</p>

***
<!-- Templated header from the pbsv github page: https://github.com/PacificBiosciences/pbsv -->

PAV is a tool for discovering variation using assembled genomes aligned to a reference. It is designed explicitly for
phased assemblies, however, it can be used for squashed assemblies by providing an empty FASTA for the second haplotype.

PAV was developed for the Human Genome Structural Variation Consortium (HGSVC)

Ebert et al., “Haplotype-Resolved Diverse Human Genomes and Integrated Analysis of Structural Variation,”
Science, February 25, 2021, eabf7117, https://doi.org/10.1126/science.abf7117.

## Installing PAV

### Pull PAV code

Choose a clean directory to install PAV ("install directory"). This will be a clean directory, and analyses will run in
separate locations.

Pull git and submodules:
`git clone --recursive https://github.com/EichlerLab/pav.git`

There is nothing to compile, the submodule libraries are pure Python.

### Dependencies 

PAV requires Python 3 with the following Python libraries installed:
1. BioPython
1. intervaltree
1. matplotlib
1. numpy
1. pandas
1. pysam
1. scipy

Note: We have had problems with pysam 0.16 crashing, which first appears at the `align_get_read_bed` step. We have
had success with pysam 0.15.2.

Command line tools needed:
1. minimap2 (default aligner)
1. lra (optional alternate aligner)
1. samtools
1. bedToBigBed (optional, from UCSC browser tools)

Optional tools: The aligner defaults to minimap2, but PAV also supports LRA (Jingwen Ren and Mark Chaisson, USC, not yet
unpublished). PAV only requires the aligner you use, both do not need to be installed.

PAV can create browser tracks for variant calls, and if this feature is used, the UCSC browser tools are needed
(bedToBigBed). Pre-compiled binaries can be obtained from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/


## Configuring PAV

Once PAV is installed, move to a clean directory ("analysis directory") where PAV will be run. Do not try to run PAV
from the install directory (where `Snakefile` is found). The following configuration steps will be done in the
new analysis directory. This way, PAV can be deployed once and run on any number of projects without duplicating code.


### Base config: config.json

A JSON configuration file, `config.json`, configures PAV. Default options are built-in, and the only required option is
"reference" pointing to a reference FASTA file variants are called against.

Example:

    {
      "reference": "/path/to/hg38.no_alt.fa.gz"
    }

Note: The HGSVC reference can be found in
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/

A no-ALT version of a reference is essential for long read assemblies. A reference containing alternate loci, decoys, or
patches may produce unexpected results and will likely lead to a loss of sensitivity.

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

There are two ways to run PAV.

1. Run scripts: Execute through `rundist` and `runlocal`. This requries a little more configuration, but allows PAV to be quickly
executed for a number of projects. This is the recommended way to run PAV.
1. Snakemake direct: Execute by calling snakemake directly.

Both methods are outlined below.

### Run scripts

To execute Snakemake, go to the run directory and link `rundist` and `runlocal` from the PAV install directory

Example:

    ln -s /path/to/PAV/1.0.0/rundist ./
    ln -s /path/to/PAV/1.0.0/runlocal ./

`rundist` will be used to distribute over a cluster, and `runlocal` will be used to run the pipeline in the current
session. Configuring `rundist` will require some knowledge for distributing over a cluster.

Both scripts setup some control variables and pass control over to a user-defined script. `rundist` will search for
`config/rundist.sh`, and `runlocal` will search for `config/runlocal.sh`. It first searches the run directory, then
it searches the PAV install directory. Once found, it calls that script to carry out calling Snakemake. Generally,
you would add your run scripts to the PAV install directory so it could be run for any number of projects. Some runs
may need custom resources not typical for other projects, and for those, you can override the script in the PAV
install directory with one in the run directory.

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

## Additional configuration options

Additional PAV configuration options for `config.json` in the PAV run directory or command-line
(e.g. `--config config_file="..."`).

The configuration parameter the expected type are given in brackets after the
parameter name. For example, type `int` would expect to take a string representation of a number (e.g. "12") and
translate it to an integer.

Type `bool` is boolean (True/False) values and is not case sensitive:
* True: "true", "1", "yes", "t", and "y"
* False: "false", "0", "no", "f", and "n"
* Other values generate an error.

### Base PAV

* config_file [config.json; string]: The runtime configuration file. This option must be supplied to PAV as a
  command-line parameter to have an effect (PAV reads the config file when it starts). By default, PAV searches for
  `config.json` in the run directory (not the PAV install directory where `Snakefile` is found). Using this option is not
  recommended, but it could be used to run different profiles (e.g. a different number of alignment threads). Changing
  some options after PAV has partially run, such as numbers of batches (see below) may cause variant call dropout
  or failures. 
* env_source [setenv.sh; string]: Source this file before running shell commands. If this file is present, it is sourced by
  prepending to all "shell" directives in Snakemake by calling: `shell.prefix('set -euo pipefail; source setenv.sh; ')`.
  Note that all shell commands evoke BASH strict mode regardless of the `setenv.sh`
  (See http://redsymbol.net/articles/unofficial-bash-strict-mode/).
* assembly_table [assemblies.tsv; string]: Assembly table file name (see "Assembly Input").
* reference [NA; string]: Location of the reference FASTA file, which may be uncompressed or gzip compressed. This
  parameter is required for all PAV runs.

### Alignments

* map_threads [12; int]: Number of threads the aligner will use.
* chrom_cluster [False; bool]: Use only if contigs are clustered by chromosome and if everything before the first "_" is
  the cluster name (e.g. "chr1_tig00001"). This feature was develoed to incoroprate QC information from PGAS
  (Porubski 2020, https://doi.org/10.1126/science.abf7117; Ebert 2021, https://doi.org/10.1126/science.abf7117), which
  clusters contigs by chromosome and haplotype, although the chromosome name is not known during the assembly (it is
  reference-free). If this option is set, PAV will fill the "CLUSTER_MATCH" field (variant call BED file). Each
  variant is called from a contig, and each contig is assinged a cluster. The cluster is assigned to a reference
  chromosome if 85% of the cluster maps to that single chromosome (not counting alignments less than 1 kbp in
  reference space). If the cluster is assigned to a chromosome, then a True/False value is given to this field if the
  chromosome and cluster match. If This feature is disabled (default) or the cluster was not assigned, the
  "CLUSTER_MATCH" field for the variant call will be empty (NA).
* min_trim_tig_len [1000, int]: Ignore contig alignments shorter than this size. This is the number of reference bases
  the contig covers.
* aligner [minimap2; string]: PAV supports minimap2 (aligner: minimap2) and LRA (aligner: lra). LRA is under development
  with the Chaisson lab (USC).


### Diploid callset merging

For phased assemblies where both input FASTA files for one sample are present (neither is empty), the final diploid
callset is found by merging two independent callsets for one haplotype.

The merging process can be found in methods for the HGSVC paper PAV was released in
(Ebert 2021, https://doi.org/10.1126/science.abf7117). SNVs are merged by exact match only (same chromosome position and
alternate base). For all others, two variants are intersected (i.e. homozygous) in a three-step process:

1. Exact match by variant ID.
1. 50% reciprocal overlap
    1. For insertions, the "END" is set to SV position + SV length.
1. Within 200 bp offset and 50% reciprocal size overlap by size only
    1. Reciprocal size overlap is reciprocal overlap if variants were shifted to maximal overlap (e.g. same starting
        position)
    1. Offset is defined by the minimum of the start position offset and end position offset (with insertion end
       position + length).

These parameters are tunable:

* ro_min [0.5; float]: Minimim reciprocal overlap and reciprocal size overlap.
* offset_max [200; float]: Maximum offset.
* merge_threads [12; int]: Number of threads for each merging step.
* merge_by_chrom [True; bool]: If True, do merging by chromosome, then concatenate results. This causes PAV to spaw
  more smaller jobs, which may significantly improve parallel performance when distributed over a cluster. This
  parameter should not impact results.


### Inversions

* inv_region_limit [1200000; int]: When scanning for inversions, the region is expanded until an inversion flanked
  by clean reference sequence is found. Expansion stops when a full inversion is found, no evidence of an inversion
  is seen in the region, or expansion hits this limit (number of reference bases covered).

#### Large inversions

These inversion calls break the alignment into mulitple alignment records. PAV will first try to resolve these
through its k-mer density 

* inv_k_size [31; int]: Use k-mers of this size when resolving inversion breakpoints. Contig k-mers are assigned to
  a reference state (FWD, REV, FWD+REV, NA) by matching k-mers of this size to reference k-mer sets over the
  same region.
* inv_threads_lg [12; int]: Like `inv_threads`, but for large (alignment-truncating) inversion events.
* lg_batch_count [10; int]: Run large variant detection (variants truncating alignments) in this many batches and
  run each batch as one job. Improves the speed of large inversion detection in a cluster environment.

#### Small inversions

Small inversions typically do not cause the aligner to break the contig into multiple records. Instead, it aligns
through the inversion leaving signatures of variation behind, but no inversion call. PAV recognizes these
signatures, such as SV insertions and deletions in close proximity with the same size and resolves the
inversion event. When an inversion is found, the false variants from those signatures are replaced with the
inversion call.

Parameters controlling inversion detection for small inversions:

* inv_threads [4; int]: Scan for inversions using this many simultaneous threads to compute smoothed density
  plots for k-mer states (for breakpoint and structure resolution). This parameter is specifically for resolving
  inversions within aligned contigs (i.e. did not break contig into multiple alignment records).
* inv_min_expand [1; int]: When searching for inversions, expand the search window at least this many times before
  giving up on a region. May be used to increase sensitivity for some inversions, but increasing this
  parameter penalizes performance heavily trying to resolve complex loci where there are no inversions.


##### Inversion signatures

Signatures of small inversions include:
1. Matched SV insertions and deletions in close proximity with the same size (i.e. insertion matched with deletion).
1. Matched indel insertions and deletionsin close proximityd with the same size.
1. Clusters of SNVs.
1. Clusters of indels.

* inv_sig_filter [svindel; string]: Determine which inversion signature patterns PAV should try to resolve into an
  inversion call. See the four recognized signatures above.
    * single_cluster: Accept signatures from only clustered SNVs or clustered indels. This will cause PAV to search
      sites and maximize sensitivity for small inversions. There is a steep performance cost, and possibly a
      sharp increase in false-positive inversions as well. Use this parameter with caution and purpose.
    * svindel: Accept signatures with matched SVs or matched indels.
    * sv: Only accept signatures with matched SVs.
* inv_sig_merge_flank [500; int]: When searching for clustered SNVs or indels, cluster variants within this many bp.
* inv_sig_batch_count [60; int]: Batch signature regions into this many batches for the caller and distribute
  each batch as one job. Improves cluster parallelization, but has no performance effect for
  non-distributed jobs.
* inv_sig_insdel_cluster_flank [2; float]: For each insertion, multiply the SV length by this value and
  search for a matching SV deletion.
* inv_sig_insdel_merge_flank [2000; int]: Merge clusters within this distance.
* inv_sig_cluster_svlen_min [4; int]: Discard indels less than this size when searching for clustered signatures.
* inv_sig_cluster_win [200; int]: Cluster variants within this many bases
* inv_sig_cluster_win_min [500; int]: Clustered window must reach this size (first clustered variant to last
  clustered variant).
* inv_sig_cluster_snv_min [20; int]: Minimum number if SNVs in a window to be considered a SNV-clustered site.
* inv_sig_cluster_indel_min [10; int]: Minimum number of indels in a window to be considered an indel-clustered
  site.

##### Density resolution parameters

In a region that is being scanned for an inversion, PAV knows the reference locus, the matching contig locus,
and the reference orientation (i.e. if the contig was reverse complemented when it was mapped).

PAV first extracts a list of k-mers from the scanned contig region. Over the matching reference region, it gets
a set of k-mers that match the contig k-mers and a set of k-mers that match the reverse-complemented contig
k-mers.

Each contig k-mer is then assigned a state:

* FWD: Contig k-mer matches a reference k-mer in the same orientation.
* REV: Contig k-mer matches a reference k-mer if it is reverse-complemented.
* FWDREV: Contig k-mer and it's reverse complement both match reference k-mers.
  * Common inside inverted repeats flanking NAHR-driven inversions.
* NA: K-mer does not match a reference k-mer in the region.

PAV then removes any k-mers with an NA state and re-numbers them contiguously. This allows the following steps to
be less biased by variation that emerged within the inversion since the inversion event or contcurrently with
the inversion event. For example, an MEI insertion in the inversion will not cause inversion detection to fail. PAV
stores the original k-mer coordinates, which are restored later.

PAV then generates a density function for each state (FWD, REV, FWDREV) and applies it to the contig k-mers.
Density is computed and scaled between 0 and 1 (max value for any k-mer is 1). Using the smoothed density, PAV
determines the internal structure of the inversion by setting breakpoins where the maximum state changes. For
example, a typical inversion flanked by inverted repeats will traverse through states
(FWD > FWDREV > REV > FWDREV > FWD). If PAV sees this pattern, it sets outer breakpoints at the FWD/FWDREV
transition and inner breakpoints at the FWDREV/REV transitions. Note that some inversions generated by other
mechanisms will not have FWDREV states, and the inner/outer breakpoints will be set to the same location where
FWD transitions to either FWDREV or REV on each side. PAV requires max states tho be FWD on each side and
REV occuring at least once inside the locus to call an inversion. Very complex patterns are sometimes seen
inside, such as multiple FWDREV/REV transitions (common in MEI-dense loci).

Once PAV sees the necessary state transitions, it moves back to the original contig coordinates (recall that NA
states were removed for density computations), annotates the inner and outer breakpoints in both reference and
contig space, and emits the inversion call.

Because computing density values is very compute heavy, PAV implements an interpolation scheme to avoid computing
densities for each location in the contig. By default, it computes every 20th k-mer density and interpolates
the values between them. If the maximum state changes or if there is a significant change in the density within
the interpolated region, the actual density is computed for each k-mer in the window.

A parameter, `srs_list` is a list of lists in the JSON configuration file. It can be used to allow PAV to
search for very large inversions while fine-tuning the interpolation distance. Generally, larger inversions can
tolerate more interpolation, however, it is possible that small changes in structure would be missed. We have
observed little gains in our own work trying to tune this way, but it can be helpful if the reference and
sample are significantly diverged.

Each element in `srs_list` is a pair of region length and interpolation distance. See example below. 

* srs_list [None; list of lists]: Set a dynamic smoothing limits based on the size of the region being scanned.
  Larger regions with interpolate over larger distances.


    "srs_list": [
        [0, 20],
        [20000, 40],
        [60000, 80],
        [100000, 120],
        [500000, 400],
        [800000, 600]
    ],

In this example, a region up to 20 kbp will compute denisty on every 20th k-mer, a region from 20 kbp up to 60 kbp
will compute every 40th k-mer, and so on. If the first record does not start with 0, then it is automatically
extended down to 0 (e.g. if the first element is "[100, 20]", it is the same as "[0, 20]"). The last window size
is extended up to infinity (e.g. a 1 Mbp inversion in the example above would compute density for every 600th k-mer).

Recall that PAV will fill actual density values if it sees a state change or a large change in state, so large values
for the second parameter will give diminishing returns.

## Cite PAV

Ebert et al., “Haplotype-Resolved Diverse Human Genomes and Integrated Analysis of Structural Variation,”
Science, February 25, 2021, eabf7117, https://doi.org/10.1126/science.abf7117 (PMID: 33632895).

Audano et al., "PAV: An assembly-based approach for discovering structural variants, indels, and point mutations
in long-read phased genomes," ASHG Annual Meeting, October 20 (10:45 - 11:00 AM), PrgmNr 1160

## Contact

Please open a case on the Github page for problems.

You may also contact Peter Audano directly (e-mail omitted to reduce SPAM). PAV was developed in the lab of
Dr. Evan Eichler at the University of Washington and is currently maintained in the lab of Dr. Christine beck at The
Jackson Laboratory.
