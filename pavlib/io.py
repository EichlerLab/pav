"""
I/O utilities.
"""

import gzip
import os
import pysam

class PlainOrGzFile:
    """
    Read a plain or a gzipped file using context guard.

    Example:
        with PlainOrGzReader('path/to/file.gz'): ...
    """

    def __init__(self, file_name, mode='rt'):

        if file_name is None:
            raise RuntimeError('File name is missing')

        file_name = file_name.strip()

        if not file_name:
            raise RuntimeError('File name is empty')

        if mode is not None:
            mode = mode.strip()

            if not mode:
                mode = 'rt'
        else:
            mode = 'rt'

        self.file_name = file_name
        self.is_gz = file_name.strip().lower().endswith('.gz')

        self.mode = mode

        self.file_handle = None

    def __enter__(self):

        if self.is_gz:
            self.file_handle = gzip.open(self.file_name, self.mode)
        else:
            self.file_handle = open(self.file_name, self.mode)

        return self.file_handle

    def __exit__(self, exc_type, exc_value, traceback):

        if self.file_handle is not None:
            self.file_handle.__exit__(exc_type, exc_value, traceback)
            self.file_handle = None

class FastaReader:
    """
    Accepts a FASTA file name or an open FASTA file (pysam.FastaFile) and provides a context-guard for the file.

    Examples:
        with FastaReader('path/to/file.fa.gz'): ...
        with FastaReader(fasta_file): ...  # fasta_file is a pysam.FastaFile
    """

    def __init__(self, file_name):

        if file_name is None:
            raise RuntimeError('File name or open FASTA file is missing')

        if isinstance(file_name, str):
            file_name = file_name.strip()

            if not file_name:
                raise RuntimeError('File name is empty')

            if not os.path.isfile(file_name):
                raise RuntimeError(f'File name does not exist or is not a regular file: {file_name}')

            self.is_pysam = False

            self.file_name = file_name
            self.file_handle = None

        elif isinstance(file_name, pysam.FastaFile):
            self.is_pysam = True

            self.file_name = "<pysam.FastaFile Object>"
            self.file_handle = file_name

        else:
            raise RuntimeError(f'File name or open FASTA file is not a string or a pysam.FastaFile: {file_name} (type "{type(file_name)}")')

        self.file_handle = None

        self.is_open = False

    def __enter__(self):

        if self.is_open:
            raise RuntimeError(f'Enter called: File is already open by this context guard: {self.file_name}')

        if not self.is_pysam:
            self.file_handle = pysam.FastaFile(self.file_name)

        self.is_open = True

        return self.file_handle

    def __exit__(self, exc_type, exc_value, traceback):

        if not self.is_open:
            raise RuntimeError(f'Exit called: File is not open by this context guard: {self.file_name}')

        if not self.is_pysam:
            self.file_handle.__exit__(exc_type, exc_value, traceback)

        self.is_open = False
