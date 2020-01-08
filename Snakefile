"""
Make variant calls from aligned contigs.
"""


import collections
import gzip
import intervaltree
import numpy as np
import os
import pandas as pd
import pysam
import sys
import re

from Bio import SeqIO
import Bio.bgzf

#
# Definitions
#

REF_FA = '/net/eichler/vol26/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa'

REF_FAI = '/net/eichler/vol26/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa.fai'

WINDOW_SIZE = 500

SMRTSV_DIR = '/net/eichler/vol27/projects/structural_variation/nobackups/tools/smrtsv3/201910/'

SVPOP_DIR = '/net/eichler/vol27/projects/structural_variation/nobackups/tools/svpop/201910'

PRINT_GAPS_LOC = 'scripts/PrintGaps.py'

ASM_ALN_DICT = {
    'GARG_NA12878_PG_PH_PN': '/net/eichler/vol27/projects/phased_hifi_strandseq/nobackups/assemblies/GARG_NA12878_PG_PH_PN/assembly_{hap}.cram',
    'DP_HG00733_PG_PH_PN': '/net/eichler/vol27/projects/phased_hifi_strandseq/nobackups/assemblies/DP_HG00733_PG_PH_PN/assembly_{hap}.cram',
    'DP_HG00733_PG_PH_PRR': '/net/eichler/vol27/projects/phased_hifi_strandseq/nobackups/assemblies/DP_HG00733_PG_PH_PRR/assembly_{hap}.cram',
    'DP_HG00733_CNU_PH_PRR': '/net/eichler/vol27/projects/phased_hifi_strandseq/nobackups/assemblies/DP_HG00733_CNU_PH_PRR/assembly_{hap}.cram',
    'DP1.1_HG00733_PG_PH_PRR': '/net/eichler/vol27/projects/phased_hifi_strandseq/nobackups/assemblies/DP1.1_HG00733_PG_PH_PRR/assembly_{hap}.cram'
}


PRINT_GAPS = os.path.join(SMRTSV_DIR, PRINT_GAPS_LOC)

SINGLE_FILTER_PATTERN = '../consensus_filter/results/{sample}/win_single_merged.bed'


#
# Parameters
#

# Merge aligned contigs within this distance when finding consensus regions.
CONTIG_ALIGN_MERGE_DIST = '500'

# Max variant offset for clustering (RO or size-RO-offset)
OFFSET_MAX = 200

# Max reciprocal overlap or size reciprocal-overlap for clustering (RO or size-RO-offset)
RO_MIN = 0.5


#
# SVPOP imports
#

sys.path.append(SVPOP_DIR)

import analib.variant
import analib.seq


#
# Rules
#

localrules: pg_all

shell.prefix('set -euo pipefail; ')


### Process calls ###

# pg_all
#
# Make all files.
rule pg_all:
    input:
        bed=expand(
            '{asm_name}/bed/{vartype}_{svtype}_h12.bed.gz',
            asm_name=ASM_ALN_DICT.keys(),
            vartype=('sv', 'indel'),
            svtype=('ins', 'del')
        ),
        bed_snv=expand(
            '{asm_name}/bed/snv_snv_h12.bed.gz',
            asm_name=ASM_ALN_DICT.keys(),
        )


### PrintGaps Post-processing ###

# pg_variant_bed
#
# Make variant BED and FASTA. Write all variant calls regardless of consensus loci.
rule pg_variant_bed:
    input:
       bed1='temp/{asm_name}/pg/clustered/{vartype}_{svtype}_h1.bed.gz',
       bed2='temp/{asm_name}/pg/clustered/{vartype}_{svtype}_h2.bed.gz',
       con_bed='{asm_name}/align/depth_1/regions_h12.bed'
    output:
        bed='{asm_name}/bed/{vartype}_{svtype}_h12.bed.gz',
        fa='{asm_name}/bed/fa/{vartype}_{svtype}_h12.fa.gz'
    params:
        sge_opts='-pe serial 6 -l mfree=4G -l h_rt=96:00:00'
    run:
        
        # Get configured merge definition
        if wildcards.vartype == 'snv':
            config_def = 'nrid'
        else:
            config_def = 'nr:szro={}:offset={}:roor'.format(int(RO_MIN * 100), OFFSET_MAX)
        
        print('Merging with def: ' + config_def)
        sys.stdout.flush()
        
        # Merge
        df = analib.svmerge.merge_variants(
            bed_list=[input.bed1, input.bed2],
            sample_names=['h1', 'h2'],
            strategy=config_def,
            threads=6
        )
        
        # Restructure columns
        del(df['DISC_CLASS'])
        
        df.columns = [re.sub('^MERGE_', 'HAP_', val) for val in df.columns]
        
        # Load consensus regions
        consensus_tree = collections.defaultdict(intervaltree.IntervalTree)
        
        df_con = pd.read_csv(input.con_bed, sep='\t', header=None, names=('#CHROM', 'POS', 'END'))
        
        for index, row in df.iterrows():
            consensus_tree[row['#CHROM']][row['POS']:row['END']] = True
        
        # Define a function to annotate consensus regions
        def con_match(row):
            for interval in consensus_tree[row['#CHROM']].overlap(row['POS'], row['END']):
                if (interval.begin <= row['POS']) and (interval.end >= row['END']):
                    return True
                    
            return False
        
        # Annotate consensus regions
        df['HAP_CONSENSUS'] = df.apply(con_match, axis=1)
        
        # Save SEQ as a FASTA
        if 'SEQ' in df.columns:
            with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
                SeqIO.write(analib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')
            
            del(df['SEQ'])
        else:
            with open(output.fa, 'w') as out_file:
                pass
        
        # Save BED
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# pg_cluster_merge
#
# Merge variants from each chromosome.
rule pg_cluster_merge:
    input:
        bed=expand('temp/{{asm_name}}/pg/clustered/by_chrom/{{vartype}}_{{svtype}}_{{hap}}/{chrom}.bed.gz', chrom=analib.ref.get_df_fai(REF_FAI).index)
    output:
        bed='temp/{asm_name}/pg/clustered/{vartype}_{svtype}_{hap}.bed.gz'
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|snv',
        hap='h1|h2'
    params:
        sge_opts='-l mfree=24G -l h_rt=12:00:00'
    run:
        
        # Get file list, exclude 0-byte files
        file_list = [file_name for file_name in input.bed if os.stat(file_name).st_size > 0]
        
        pd.concat(
            [
                pd.read_csv(file_name, sep='\t') for file_name in file_list
            ],
            axis=0
        ).sort_values(
            ['#CHROM', 'POS']
        ).to_csv(
            output.bed, sep='\t', index=False, compression='gzip'
        )

# pg_cluster
#
# Cluster variants.
rule pg_cluster:
    input:
        bed='temp/{asm_name}/pg/raw/{vartype}_{hap}.bed',
        bed_tile='{asm_name}/align/central_tiling_path_{hap}.bed'
    output:
        bed='temp/{asm_name}/pg/clustered/by_chrom/{vartype}_{svtype}_{hap}/{chrom}.bed.gz',
        bed_dropped='{asm_name}/bed/dropped/aux_{vartype}_{svtype}_{hap}_{chrom}.bed.gz'
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|snv',
        hap='h1|h2'
    params:
        sge_opts='-l mfree=8G -l h_rt=72:00:00'
    run:
        
        # Read variants
        print('Reading...')
        sys.stdout.flush()
        
        df = pd.read_csv(input.bed, sep='\t')
        
        df = df.loc[df['#CHROM'] == wildcards.chrom]
    
        df = df.loc[df['SVTYPE'] == wildcards.svtype.upper()].copy()
        
        if df.shape[0] == 0:
            print('Empty variant set')
            
            with open(output.bed, 'w') as out_file:
                pass
            
            df.to_csv(output.bed_dropped, sep='\t', index=False)
            
            return
        
    
        if wildcards.vartype == 'snv':

            df['REF'] = df['REF'].apply(lambda val: val.upper())
            df['ALT'] = df['ALT'].apply(lambda val: val.upper())
            
            df['QUERY_END'] = df['QUERY_POS'] + 1
            
            df['SVLEN'] = 1
        
        # Transform columns (by ref), set ID, sort
        df['END'] = df.apply(lambda row: (row['POS'] + 1) if row['SVTYPE'] == 'INS' else row['END'], axis=1)
        df['ID'] = analib.variant.get_variant_id(df)
        df.sort_values(['#CHROM', 'POS'], inplace=True)
        
        df['HAP'] = wildcards.hap
        
        df = analib.variant.order_variant_columns(df, tail_cols=['HAP', 'QUERY_ID', 'QUERY_POS', 'QUERY_END', 'QUERY_STRAND'])
        
        # Remove duplicate ID/contig pairs (can occur if the same contig is mapped multiple times to a reference location,
        # e.g. a large duplication relative the reference. Happens rarely, but will crash the pipeline.
        
        df.drop_duplicates(['ID', 'QUERY_ID'], inplace=True)
        
        df.reset_index(drop=True, inplace=True)
        
        # Read tiling BED and make tree
        print('Tiling...')
        sys.stdout.flush()
        
        df_tile = pd.read_csv(input.bed_tile, sep='\t')
        
        tiling_tree = collections.defaultdict(intervaltree.IntervalTree)
        
        for index, row in df_tile.iterrows():
            if row['END'] - row['POS'] > 0:
                tiling_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']
        
        def get_central_id(row):
            central_set = tiling_tree[row['#CHROM']][row['POS']:row['END']]
            
            if len(central_set) < 1:
                return False
            
            return row['QUERY_ID'] == list(central_set)[0].data
        
        df['IS_CENTRAL'] = df.apply(get_central_id, axis=1)
        
        # Setup clustering data structures
        print('Setup clustering...')
        
        cluster_key = set()  # Primary variants (from centrally located variants)
        dropped_key = set()  # Aux variants that do not intersect a central variant
        cluster_support = collections.defaultdict(set)
    
        n_processed = 0

        cluster_tree = collections.defaultdict(intervaltree.IntervalTree)
        
        # Separate variants into central and auxiliary
        # * Central: Located in the central tiling path: when contigs overlap, choose a variant
        #   from the contig where it is further from an alignment end.
        # * Auxiliary: Variants in a region with multiple contig alignments, but closer to the end of
        #   a contig than another variant. If these intersect a central variant, annotate the central variant.
        #   For auxiliary variants that do not intersect a central variant, discard.
        df_central = df.loc[df['IS_CENTRAL']]
        df_aux = df.loc[~ df['IS_CENTRAL']]
        
        df_central.drop_duplicates(['ID'], inplace=True)  # Contigs with multiple alignments generate the same variant more than once
        
        # Check for no central variants
        if df_central.shape[0] == 0:
            print('Empty variant set')
            
            with open(output.bed, 'w') as out_file:
                pass
            
            df.to_csv(output.bed_dropped, sep='\t', index=False)
            
            return
        
        # Add primary calls from most centrally-located contig alignments
        for index, row in df_central.iterrows():
            cluster_tree[row['#CHROM']][row['POS']:row['END']] = index
            cluster_key.add(index)
        
        # Cluster aux calls
        print('Clustering AUX...')
        
        for index, row in df_aux.iterrows():
            region_start = row['POS'] - OFFSET_MAX
            region_end = row['END'] + OFFSET_MAX
            
            # Report progress
            if n_processed % 1000 == 0:
                print('\t\t* {} of {}'.format(n_processed + 1, df_aux.shape[0]))
            
            n_processed += 1
            
            # Get cluster
            cluster_set = cluster_tree[row['#CHROM']][region_start:region_end]
            
            # Check for no cluster intersection
            if not cluster_set:
                dropped_key.add(index)
                continue
            
            cluster_id_list = [interval.data for interval in cluster_set]
            
            # Check for exact ID match
            for cluster_id in cluster_id_list:
                if row['ID'] == df.loc[cluster_id, 'ID']:
                    cluster_support[cluster_id].add(index)
                    index = -1  # Flag to stop processing this variant
                    break
            
            # Drop SNVs unless they match by ID
            if wildcards.vartype == 'snv' and index != -1:
                dropped_key.add(index)
                continue
            
            # Matched variant by ID, stop processing
            if index < 0:
                continue
            
            # Intersect this variant with the cluster
            df_source = df.loc[[index]]
            df_target = df.loc[cluster_id_list]
            
            df_target.drop_duplicates('ID', inplace=True)
            
            df_intersect = analib.variant.nearest_by_svlen_overlap(
                df_source, df_target,
                szro_min=RO_MIN,
                offset_max=OFFSET_MAX
            )
            
            if not df_intersect.shape[0] > 0:
                # No match, drop aux variant
                dropped_key.add(index)
                
            else:
                # Add support for existing cluster.
                cluster_support[
                    df_source.loc[df_source['ID'] == df_intersect.iloc[0]['ID']].index[0]
                ].add(index)
        
        # Process clusters into a callset
        print('Post-cluster merging...')
        sys.stdout.flush()
        
        df_merge = df.loc[cluster_key]
        
        df_merge_support = pd.concat(
            [
                pd.Series(
                    [
                        index,
                        1 + len(cluster_support[index]),
                        ','.join(
                            [
                                df_merge.loc[index]['QUERY_ID']
                            ] + (list(
                                df.loc[cluster_support[index], 'QUERY_ID']
                            ) if cluster_support[index] else [])
                        ),
                        ','.join(
                            [
                                '{QUERY_POS}:{QUERY_END}'.format(**df_merge.loc[index])
                            ] + (list(
                                df.loc[cluster_support[index]].apply(lambda row: '{QUERY_POS}:{QUERY_END}'.format(**row), axis=1)
                            ) if cluster_support[index] else [])
                        ),
                        ','.join(
                            [
                                df_merge.loc[index]['QUERY_STRAND']
                            ] + (list(
                                df.loc[cluster_support[index], 'QUERY_STRAND']
                            ) if cluster_support[index] else [])
                        )
                        
                    ],
                    index=['INDEX', 'TIG_N', 'TIG_SUPPORT', 'TIG_COORD', 'TIG_STRAND']
                ) for index in df_merge.index
            ],
            axis=1
        ).T
        
        df_merge_support.set_index('INDEX', inplace=True)
        
        df_merge = pd.concat([df_merge, df_merge_support], axis=1).sort_values(['#CHROM', 'POS'])
        
        # Move SEQ to end
        if 'SEQ' in df_merge.columns:
            df_merge = analib.variant.order_variant_columns(df_merge, tail_cols=['SEQ', ])
        
        df_dropped = df.loc[sorted(dropped_key)]
        
        # Save BED
        print('Writing...')
        sys.stdout.flush()
        
        df_merge.to_csv(output.bed, sep='\t', index=False, compression='gzip')
        
        # Save dropped bed
        df_dropped.to_csv(output.bed_dropped, sep='\t', index=False, compression='gzip')

# pg_tiling_bed
#
# Make a tiling path for the most central aligned locations in each contig. This defines the regions
# for each contig where variants should be called.
rule pg_tiling_bed:
    input:
        bed='{asm_name}/align/aligned_tig_{hap}.bed.gz'
    output:
        bed='{asm_name}/align/central_tiling_path_{hap}.bed'
    wildcard_constraints:
        hap='h(1|2)'
    params:
        sge_opts='-l mfree=4G -l h_rt=4:00:00'
    run:
        
        # Read contig alignments and prioritize
        df_tig = pd.read_csv(input.bed, sep='\t')

        df_tig.sort_values(['#CHROM', 'POS'], inplace=True)

        df_tig.reset_index(drop=True, inplace=True)

        this_chrom = None

        last_interval = None  # [0] POS, [1] END, [2] ID

        interval_list = list()

        # Traverse sorted alignment records
        for index, row in df_tig.iterrows():
            
            # Handle chrom change
            if row['#CHROM'] != this_chrom:
                if last_interval is not None:
                    interval_list.append(pd.Series(
                        [this_chrom, last_interval[0], last_interval[1], last_interval[2]],
                        index=['#CHROM', 'POS', 'END', 'ID']
                    ))
                
                this_chrom = row['#CHROM']
                last_interval = None
            
            # Compare to last interval
            if last_interval is not None:
            
                if row['END'] <= last_interval[1]:
                    continue  # Alignment is shorter than the last interval, cannot factor into tiling path
                
                if last_interval[1] > row['POS']:
                    # Contig alignments overlap
                    
                    # Find endpoint
                    endpoint = (last_interval[1] + row['POS']) // 2
                    
                    interval_list.append(pd.Series(
                        [this_chrom, last_interval[0], endpoint, last_interval[2]],
                        index=['#CHROM', 'POS', 'END', 'ID']
                    ))
                    
                    # Set last interval
                    last_interval = (endpoint, row['END'], row['ID'])
                
                else:
                    # Contig alignments do not overlap
                    interval_list.append(pd.Series(
                        [this_chrom, last_interval[0], last_interval[1], last_interval[2]],
                        index=['#CHROM', 'POS', 'END', 'ID']
                    ))
                    
                    # Set last interval
                    last_interval = (row['POS'], row['END'], row['ID'])
            
            else:
                # Set last interval
                last_interval = (row['POS'], row['END'], row['ID'])
        
        # Add final interval
        if last_interval is not None:
            interval_list.append(pd.Series(
                [this_chrom, last_interval[0], last_interval[1], last_interval[2]],
                index=['#CHROM', 'POS', 'END', 'ID']
            ))
        
        # Merge and write BED
        df_tile = pd.concat(interval_list, axis=1).T
        
        df_tile.to_csv(output.bed, sep='\t', index=False)


### PrintGaps ###

# pg_sv
#
# Get SVs for one alignment.
rule pg_sv:
    input:
        bam=lambda wildcards: ASM_ALN_DICT[wildcards.aln_name]
    output:
        bed='temp/{aln_name}/pg/raw/sv_{hap}.bed'
    params:
        sge_opts='-l mfree=4G -l h_rt=24:00:00'
    shell:
        """samtools view -h {input.bam} | """
        """python3 {PRINT_GAPS} """
            """{REF_FA} /dev/stdin """
            """--condense 20 """
            """--minq 0 """
            """--outFile {output.bed}"""

# pg_indel
#
# Get SVs for one sample.
rule pg_indel:
    input:
        bam=lambda wildcards: ASM_ALN_DICT[wildcards.aln_name]
    output:
        bed='temp/{aln_name}/pg/raw/indel_{hap}.bed'
    params:
        sge_opts='-l mfree=4G -l h_rt=24:00:00'
    shell:
        """samtools view -h {input.bam} | """
        """python3 {PRINT_GAPS} """
            """{REF_FA} /dev/stdin """
            """--minLength 0 --maxLength 50 """
            """--removeAdjacentIndels """
            """--onTarget """
            """--condense 0 """
            """--minq 0 """
            """--outFile {output.bed}"""

# pg_snv
#
# Get SNVs for one sample.
rule pg_snv:
    input:
        bam=lambda wildcards: ASM_ALN_DICT[wildcards.aln_name]
    output:
        bed='temp/{aln_name}/pg/raw/snv_{hap}.bed'
    params:
        sge_opts='-l mfree=4G -l h_rt=24:00:00'
    shell:
        """samtools view -h {input.bam} | """
        """python3 {PRINT_GAPS} """
            """{REF_FA} /dev/stdin """
            """--minLength 0 --maxLength 0 """
            """--minq 0 """
            """--condense 0 """
            """--snv {output.bed} """
            """> /dev/null"""

### Genome coverage ###

# align_single_hap_win_merge
#
# Intersect haplotype regions.
rule align_single_hap_win_merge:
    input:
        bed1='{asm_name}/align/depth_1/regions_h1.bed',
        bed2='{asm_name}/align/depth_1/regions_h2.bed'
    output:
        bed='{asm_name}/align/depth_1/regions_h12.bed'
    params:
        sge_opts='-l mfree=1G -l h_rt=1:00:00'
    shell:
        """bedtools intersect -header -a {input.bed1} -b {input.bed2} """
        """> {output.bed}"""

# align_single_hap_win
#
# Get locations represented by a single contig.
rule align_single_hap_win:
    input:
        bed='{asm_name}/align/genomecov_{hap}.bed.gz'
    output:
        bed='{asm_name}/align/depth_1/regions_{hap}.bed'
    wildcard_constraints:
        hap='h(1|2)'
    params:
        sge_opts='-l mfree=2G -l h_rt=1:00:00'
    shell:
        """zcat {input.bed} | """
        """awk -vOFS="\\t" '($4 == 1) {{print $1, $2, $3}}' | """
        """bedtools merge -d {CONTIG_ALIGN_MERGE_DIST} """
        """> {output.bed}"""

# align_genomecov
#
# Get genome coverage.
rule align_genomecov:
    input:
        bed='{asm_name}/align/aligned_tig_{hap}.bed.gz'
    output:
        bed='{asm_name}/align/genomecov_{hap}.bed.gz'
    params:
        sge_opts='-l mfree=4G -l h_rt=4:00:00'
    shell:
        """{{ \n"""
        """    echo -e "#CHROM\tPOS\tEND\tDEPTH"; \n"""
        """    bedtools genomecov -bga -i {input.bed} -g {REF_FAI}; \n"""
        """}} | """
        """gzip > {output.bed}"""

# align_merge_h12_read_bed
#
# Alignment table for all reads.
rule align_merge_h12_read_bed:
    input:
        bed1='{asm_name}/align/aligned_tig_h1.bed.gz',
        bed2='{asm_name}/align/aligned_tig_h2.bed.gz'
    output:
        bed='{asm_name}/align/aligned_tig_h12.bed.gz'
    params:
        sge_opts='-l mfree=4G -l h_rt=1:00:00'
    run:
        
        # Read
        df1 = pd.read_csv(input.bed1, sep='\t', low_memory=False)
        df2 = pd.read_csv(input.bed2, sep='\t', low_memory=False)
        
        # Merge & sort
        df = pd.concat([df1, df2], axis=0)
        
        df.sort_values(['#CHROM', 'POS'], inplace=True)
        
        # Write
        df.to_csv(output.bed, sep='\t', index=False)

# align_get_read_bed
#
# Get alignment BED for one part (one aligned cell or split BAM) in one assembly.
rule align_get_read_bed:
    input:
        aln=lambda wildcards: ASM_ALN_DICT[wildcards.asm_name]
    output:
        bed='{asm_name}/align/aligned_tig_{hap}.bed.gz'
    params:
        sge_opts='-l mfree=6G -l h_rt=4:00:00'
    wildcard_constraints:
        hap='h(0|1|2)'
    run:
        
        # Get records
        clip_l = 0
        clip_r = 0

        record_list = list()

        with pysam.AlignmentFile(input.aln, 'rb') as in_file:
            for record in in_file:

                # Skipped unmapped reads
                if record.is_unmapped:
                    continue

                # Read tags
                tags = dict(record.get_tags())

                # Get clipping
                cigar_tuples = record.cigartuples

                l_index = 0 if cigar_tuples[0][0] != 5 else 1
                r_index = -1 if cigar_tuples[-1][0] != 5 else -2

                clip_l = cigar_tuples[l_index][1] if cigar_tuples[l_index][0] == 4 else 0
                clip_r = cigar_tuples[r_index][1] if cigar_tuples[r_index][0] == 4 else 0

                # Save record
                record_list.append(
                pd.Series(
                    [
                        record.reference_name,
                        record.reference_start,
                        record.reference_end,
                        record.query_name,

                        tags['RG'] if 'RG' in tags else 'NA',

                        record.mapping_quality,
                        clip_l,
                        clip_r,

                        str(record.is_reverse),
                        '0x{:04x}'.format(record.flag),
                        
                        wildcards.hap
                    ],
                    index=[
                        '#CHROM', 'POS', 'END', 'ID',
                        'RG',
                        'MAPQ', 'CLIP_L', 'CLIP_R',
                        'REV', 'FLAGS', 'HAP'
                    ]
                ))

        # Merge records
        df = pd.concat(record_list, axis=1).T

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')
