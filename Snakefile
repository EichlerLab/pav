"""
Make variant calls from aligned contigs.
"""


import collections
import gc
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
REF_FA = 'data/ref/ref.fa.gz'

REF_FAI = REF_FA + '.fai'

# Environment source file for shell commands

env_file = config.get('env_source', 'setenv.sh')

ENV_FILE = os.path.join(PIPELINE_DIR, 'config', env_file)

if not os.path.isfile(ENV_FILE):
    ENV_FILE = None

#
# Assembly library and dependency imports
#

sys.path.append(PIPELINE_DIR)  # pavlib
sys.path.append(os.path.join(PIPELINE_DIR, 'dep', 'svpop'))  # svpoplib
sys.path.append(os.path.join(PIPELINE_DIR, 'dep'))  # kanapy

import pavlib
import svpoplib
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

if ENV_FILE:
    shell.prefix('set -euo pipefail; source {}; '.format(ENV_FILE))
else:
    shell.prefix('set -euo pipefail; ')


### Wildcard constraints ###

wildcard_constraints:
    asm_name='[A-Za-z_\-0-9\.]+'

### Default rule ###

localrules: pav_all

# pav_all
#
# Make all files for all samples.
rule pav_all:
    input:
        bed=expand(
            'pav_{asm_name}.vcf.gz',
            asm_name=ASM_TABLE['NAME']
        )

### Includes ###

include: 'rules/data.snakefile'
include: 'rules/input_functions.snakefile'
include: 'rules/align.snakefile'
include: 'rules/call.snakefile'
include: 'rules/call_inv.snakefile'
include: 'rules/call_lg.snakefile'
include: 'rules/tracks.snakefile'
include: 'rules/figures.snakefile'
include: 'rules/vcf.snakefile'