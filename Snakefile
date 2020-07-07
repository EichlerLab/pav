"""
Make variant calls from aligned contigs.
"""


import collections
import gzip
import intervaltree
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pysam
import subprocess
import sys
import re
from scipy import stats
import shutil
import tempfile

from Bio import SeqIO
import Bio.bgzf


#
# Global constants
#

PIPELINE_DIR = os.path.dirname(workflow.snakefile)


#
# Parameters
#

configfile: config.get('config_file', 'config.json')


### Parameters from config ###

# Reference FASTA & FAI
REF_FA = config.get(
    'reference',
    '/net/eichler/vol26/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa'
)

if not os.path.isfile(REF_FA):
    raise RuntimeError('Reference file does not exist or is not a file: {}'.format(REF_FA))

REF_FAI = config.get(
    'ref_fai',
    REF_FA + '.fai'
)

if not os.path.isfile(REF_FAI):
    raise RuntimeError('Reference FAI file does not exist or is not a file: {}'.format(REF_FAI))


# Environment source file for shell commands
ENV_FILE = config.get(
    'env_source', None
)

if ENV_FILE is None:
    ENV_FILE = os.path.join(PIPELINE_DIR, 'config/setenv.sh')

if not os.path.isfile(ENV_FILE):
    raise RuntimeError('Shell configuration source not found: {}'.format(ENV_FILE))


#
# Assembly library and dependency imports
#

sys.path.append(PIPELINE_DIR)
sys.path.append(os.path.join(PIPELINE_DIR, 'dep'))

import asmlib
import analib
import kanapy


#
# Read sample config
#

ASM_TABLE_FILE_NAME = config.get('assembly_table', None)

if ASM_TABLE_FILE_NAME is None and os.path.isfile('assemblies.tsv'):
    ASM_TABLE_FILE_NAME = 'assemblies.tsv'

if ASM_TABLE_FILE_NAME is None:
    # Empty table
    ASM_TABLE = pd.DataFrame([], columns=('NAME', 'HAP1', 'HAP2'))

else:
    # Read table

    if not os.path.isfile(ASM_TABLE_FILE_NAME):
        raise RuntimeError('Missing assembly table: {}'.format(ASM_TABLE_FILE_NAME))

    ASM_TABLE = pd.read_csv(ASM_TABLE_FILE_NAME, sep='\t', header=0)

    if 'NAME' not in ASM_TABLE.columns:
        raise RuntimeError('Missing assembly table column: NAME')

    if ('HAP1' not in ASM_TABLE.columns or 'HAP2' not in ASM_TABLE.columns) and 'HAP0' not in ASM_TABLE.columns:
        raise RuntimeError('Assembly table must contain both "HAP1" and "HAP2" columns, or the "HAP0" column: {}'.format(', '.join(ASM_TABLE.columns)))

    DUP_ASM_LIST = [assembly_name for assembly_name, count in collections.Counter(ASM_TABLE['NAME']).items() if count > 1]

    if DUP_ASM_LIST:
        raise RuntimeError('Found {} duplicate assembly names in table "{}": {}'.format(
            len(DUP_ASM_LIST), ASM_TABLE_FILE_NAME, ', '.join(DUP_ASM_LIST)
        ))

    del(DUP_ASM_LIST)

    ASM_TABLE.set_index('NAME', inplace=True, drop=False)


#
# Rules
#

shell.prefix('set -euo pipefail; source {}; '.format(ENV_FILE))


### Wildcard constraints ###

wildcard_constraints:
    asm_name='[A-Za-z_\-0-9]+'

### Default rule ###

localrules: pg_all

# pg_all
#
# Make all files for all samples.
rule pg_all:
    input:
        bed=expand(
            'results/{asm_name}/bed/{vartype}_{svtype}_h12.bed.gz',
            asm_name=ASM_TABLE['NAME'],
            vartype=('sv', 'indel'),
            svtype=('ins', 'del')
        ),
        bed_snv=expand(
            'results/{asm_name}/bed/snv_snv_h12.bed.gz',
            asm_name=ASM_TABLE['NAME']
        )


### Includes ###

include: 'rules/input_functions.snakefile'
include: 'rules/align.snakefile'
include: 'rules/call.snakefile'
include: 'rules/call_inv.snakefile'
include: 'rules/call_lg.snakefile'
include: 'rules/tracks.snakefile'
include: 'rules/figures.snakefile'
