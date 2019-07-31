from __future__ import (
    absolute_import,
    print_function,
    unicode_literals
)
import sys


def print_stderr(message):
    """
    Simply print str message to stderr
    """
    print(message, file=sys.stderr)


def force_exit(message):
    """
    Exit the program due to some error. Print out message to stderr
    """
    print_stderr(message)
    sys.exit(1)
