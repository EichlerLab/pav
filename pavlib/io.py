"""
I/O utilities.
"""

import gzip

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
            self.file_handle.__exit__()
            self.file_handle = None
