# PAV Configuration Options

## Alignment

* aligner [minimap2]
* min_trim_tig_len [1000]
* redundant_callset [False]
* chrom_cluster [False]
* sample_delimiter [_]: For extracting sample from assembly name ("{sample}" wildcard in input path patterns). May be one of "_", ".", "-", "+", "#".

## Call

* ro_min [0.5]
* offset_max [200]
* merge_threads [12]
* merge_align [None]
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

