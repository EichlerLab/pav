#!/usr/bin/env python3

"""
Print the PAV version to standard output.
"""

import os
import sys

# Find pipeline directory
PIPELINE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

constants_file = os.path.join(PIPELINE_DIR, 'pavlib/constants.py')

if not os.path.isfile(constants_file):
    print('<unknown>')
    print('Unknown version: Cannot locate pavlib/constants.py relative to scripts/print_version.py in the PAV directory structure', file=sys.stderr)
    sys.exit(1)

with open(constants_file, 'rt') as in_file:
    exec(in_file.read())

print(get_version_string())
