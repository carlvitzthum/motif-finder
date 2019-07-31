from __future__ import absolute_import
from . import (
    mf_main,
    __version__
)
import argparse
import sys


def main():
    """
    Execute the program from the command line
    """
    main_help = ("Main program for motif finding. Provide a query file and "
                 "optional motif file(s) with the -m argument")
    parser = argparse.ArgumentParser(prog='motiffinder', description=main_help)
    parser.add_argument('queries', help="File containing settings and queries")
    motifs_help = ("File containing proteins used to generate motifs."
                   "Can specify multiple")
    parser.add_argument('-m', '--motifs', action='append', help=motifs_help)
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + __version__)
    args = parser.parse_args()
    mf_main(args.queries, args.motifs)


if __name__ == '__main__':
    main()
