# PAV Configuration Options

## Alignment

* aligner [minimap2]
* min_trim_tig_len [1000]
* redundant_callset [False]: Allow multiple haplotypes in one assembly. This option turns off trimming reference bases
  from alignments unless two alignment records are from the same contig. Can increase false calls, especially for small
  variants, but should not eliminate haplotypic variation in an assembly that is not separated into two FASTA files.
* chrom_cluster [False]
* sample_delimiter [_]: For extracting sample from assembly name ("{sample}" wildcard in input path patterns). May be one of "_", ".", "-", "+", "#".

## Call

* merge_ins [param merge_svindel]: Override default merging parameters for SV/indel insertions (INS).
* merge_del [param merge_svindel]: Override default merging parameters for SV/indel deletions (DEL).
* merge_inv [param merge_svindel]: Override default merging parameters for SV/indel inversions (INV). 
* merge_svindel [nr::exact:ro(0.5):szro(0.5,200):match]: Override default merging parameters for INS, DEL, INV
* merge_snv [nrsnv:exact]: Override default merging parameters for SNVs
* merge_threads [12]
* min_inv [300]
* max_inv [2000000]
* inv_min_svlen [300]

* inv_min [300]
* inv_max [2000000]

## Call - INV
* inv_k_size [31]
* inv_threads [4]
* inv_region_limit [None]
* inv_min_expand [None]
* inv_sig_merge_flank [500]
* inv_sig_batch_count [BATCH_COUNT_DEFAULT]
* inv_sig_filter [svindel]
* inv_sig_insdel_cluster_flank [2]
* inv_sig_insdel_merge_flank [2000]
* inv_sig_cluster_svlen_min [4]
* inv_sig_cluster_win [200]
* inv_sig_cluster_win_min [500]
* inv_sig_cluster_snv_min [20]
* inv_sig_cluster_indel_min [10]

## Call - Large SV
* inv_k_size [31]
* inv_threads_lg [12]
* inv_region_limit [None]
* lg_batch_count [10]

