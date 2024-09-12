"""
Rules for writing VCF output.
"""

import pavlib

global VCF_PATTERN
global get_config

global shell


_VCF_SVTYPE = {
    'sv': ('ins', 'del', 'inv'),
    'indel': ('ins', 'del'),
    'snv': ('snv', )
}

_VCF_INPUT_PATTERN_BED = 'results/{asm_name}/bed_{merge}/{filter}/{varsvtype}.bed.gz'

_VCF_INPUT_PATTERN_FA = 'results/{asm_name}/bed_{merge}/{filter}/fa/{varsvtype}.fa.gz'

# Make VCF file.
rule vcf_write_vcf:
    input:
        bed=lambda wildcards: [
            _VCF_INPUT_PATTERN_BED.format(asm_name=wildcards.asm_name, merge=wildcards.merge, filter='pass', varsvtype=varsvtype)
                for varsvtype in ('snv_snv', 'svindel_ins', 'svindel_del', 'sv_inv')
        ],
        fa=lambda wildcards: [
            _VCF_INPUT_PATTERN_FA.format(asm_name=wildcards.asm_name, merge=wildcards.merge, filter='pass', varsvtype=varsvtype)
                for varsvtype in ('svindel_ins', 'svindel_del')
        ],
        ref_tsv='data/ref/contig_info.tsv.gz'
    output:
        vcf='vcf/{merge}/{asm_name}.vcf.gz'
    wildcard_constraints:
        merge='merged|hap'
    run:

        # Get a dictionary of input files.
        #
        # input_dict:
        #   * key: Tuple of...
        #     * [0]: varsvtype
        #     * [1]: "pass" or "fail"
        #   * Value: Tuple of...
        #     * [0]: Variant BED file name.
        #     * [1]: Variant FASTA file name (None if no variant sequences are not used in the VCF).
        input_dict = dict()

        for varsvtype in ('snv_snv', 'indel_ins', 'indel_del', 'sv_ins', 'sv_del', 'sv_inv'):

            # Pass
            input_dict[(varsvtype, 'pass')] = (
                _VCF_INPUT_PATTERN_BED.format(
                    asm_name=wildcards.asm_name, merge=wildcards.merge, filter='pass', varsvtype=varsvtype
                ),
                _VCF_INPUT_PATTERN_BED.format(
                    asm_name=wildcards.asm_name, merge=wildcards.merge, filter='pass', varsvtype=varsvtype
                ) if varsvtype in {'svindel_ins', 'svindel_del'} else None
            )

            # Fail
            input_dict[(varsvtype, 'fail')] = (
                _VCF_INPUT_PATTERN_BED.format(
                    asm_name=wildcards.asm_name, merge=wildcards.merge, filter='fail', varsvtype=varsvtype
                ),
                _VCF_INPUT_PATTERN_BED.format(
                    asm_name=wildcards.asm_name, merge=wildcards.merge, filter='fail', varsvtype=varsvtype
                ) if varsvtype in {'indel_ins', 'indel_del', 'sv_ins', 'sv_del'} else None
            )

        # Write VCF
        pavlib.vcf.write_merged_vcf(
            asm_name=wildcards.asm_name,
            input_dict=input_dict,
            reference_file=get_config(wildcards, 'reference')
        )

        # Write tabix index if possible
        try:
            shell("""tabix {output.vcf} && touch -r {output.vcf} {output.vcf}.tbi""")
        except:
            pass
