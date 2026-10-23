"""
IO utilities for file handling.
"""

class FileHandler:
    """
    A simple context manager for opening and closing files safely.
    """

    def __init__(self, infile, mode='r', encoding='utf-8'):
        self.infile = infile
        self.mode = mode
        self.encoding = encoding
        self.file_handle = None

    def __enter__(self):
        self.file_handle = open(self.infile, self.mode, encoding=self.encoding)
        return self.file_handle

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file_handle:
            self.file_handle.close()
        return False

