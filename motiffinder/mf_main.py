# Carl Vitzthum
# For honors thesis under Andrea Tilden, spring 2016
# Main code for motif generating and finding program
from __future__ import (
    absolute_import,
    print_function,
    unicode_literals
)
from .mf_generate import MotifGenerator
from .mf_query import Query
from .mf_utils import force_exit
from . import static
import json
try:
    import importlib.resources as pkg_resources
except ImportError:
    # for Python < 3.7
    import importlib_resources as pkg_resources


def read_in(filename, motif_files):
    """
    use: input txt file containing all queries formatted as follows:

    Any lines including '#' are comment lines (not read). By convention, # should be the first character

    FIRST TWO LINES: desired output filename and motifs to be searched for in ALL queries
    (This is ignoring commented lines)

    # COMMENTS (NAME, DATE, etc.)
    File: <output_filename.tsv> THIS SHOULD BE .TSV
    motifs: <Your motifs separated by commas>
        Example: FRHGLRLHDNPALLAAL,VPR

    motif input rules:
        For multiple possible amino acids: list in brackets.
            Example: A[LKY]R may ALR, AKR, or AYR
        For repeating amino acids of variable lengths: list in curly braces, from minimum to maximum number of repeats
        separated by a dash.
            Example: AY{1-3}R expands to AYR, AYYR, or AYYYR.
        For placeholders, use 'x'. These are wildcard positions and may be filled by any amino acid.
        For optional characters, follow by '?'.
            Example: AY?R may be AYR or AR

    GENERATED MOTIFS:
        It is an option to automatically generate motifs from a number of given protein sequences (as FASTAs)
        To do this, type 'generate' in the motifs: line.
        If generating motifs, provide any number of additional text files containing
        protein sequences to the program use the `-m` option.
        Each grouping of FASTAs will be used individually to generate motifs
        Then, you must input the input the motif finding thresholds in this order:
            [minimum pairwise match length, % gaps allowed, % different aa's allowed, % similar aa's allowed,
            proportion of reference sequences consensus motif must be found in, minimum length of consensus motifs]
            The first four thresholds are alignment thresholds, the fifth is the sampling threshold, and the sixth
            is the consensus threshold.
        These numbers are int, float, float, float, float, int
        Good thresholds to start with are 8, 0, 0.2, 0.8, 1.0, 6
        In file, input like this:
        Thresholds: <THRESHOLD VALUES>

        Then continue with the rest of the input file

    Next, list all your queries (as many as you want) like this
        Query: Example gene
        Sequence: ATGGCCACGCGAGGGGCGAATGTGATTTGGTTTCGCCATGGATTGCGCCTCCA
                  TGATAATCCCGCTCTATTGGCCGCCCTCGCCGATAAGGATCAGGGTATAGCCCTAATTCC
                  CGTTTTCATATTCGATGGAGAGAGTGCAGgtaagaactgagttttcgaaaattttatttc
                  gttgcgcacacatacatgtaggtatccaaattaaattgattacaggtaaccttgaactta
                  ccttgcactttcctaagttattattcctagthjbjhb zxtttgcgaaaattgtgcacgaaaaa
                  actttgaaaatgcataattttttttaaatttataaatcttataaaaccgttaaatgctct
                  tccttaacagGTACCAAGAATGTGGG
        Exons: 1-142, 368-383
        Options: ranges, flex2

    THERE SHOULD BE AT LEAST TWO QUERIES with SEQUENCES for alignment purposes
    More info on the query fields:

    Query: <any defining info, such as gene name/position>
    Sequence: <sequence in NUCLEOTIDES with exons in uppercase in introns in lowercase, unless exon ranges are to
        be used (see options). If so, nucleotide case does not matter
        Protein sequence may also be provided. If doing so, add protein to options (see below)>
    Exons: <exon coordinates. Ranges written as start-stop, or a single nucleotide
              separate coordinates with commas. First nucleotide is 1>
    Options: <any options, see below. Separate entries by commas. Order doesn't matter>
              reference - use this query as the reference sequence for all alignments (identity similarity scores)
                    if this option is not used for any query, the first query will be reference
              ranges - use if exon ranges, rather than upper/lowercase, will be used to determine exons in the sequence
              reverse - will reverse translate the given sequence (use if - strand)
              transcript - use if NT sequence is mRNA. Will exclude this query from the exon-based analysis
              protein - use if the input sequence is a PROTEIN sequence. Will obviously skip the translation process
              flex NUM - use if you want to expand the input motifs and look for new motifs.
                    NUM is in integer between 0 and 5. (If not flex is provided, 0 is automatically used, corresponding
                    the basic input motifs). 1 causes the least expansion of input motifs searched for the given
                    query, while 5 causes the most expansion. The higher the flex score, the more time it will take
                    to compute. BE WARNED - flex 5 can take a long time. Among other things, the flex algorithm works
                    by allowing amino acid substiutions (dictated by BLOSUM matrix), allowing a limited number of gaps,
                    and allowing a limited number of repeats. Higher flex scores should be used for more distant queries,
                    or if the input motifs are quite strict.
    """
    f_name = ''
    motifs = []
    check_idx = 0
    queries = []
    query_idx = -1
    generate = False
    thresholds = []
    frames = []
    with open(filename, 'r') as infile:
        for line in infile.readlines():
            if line.startswith('#'):
                continue
            l_list = line.split(':')
            if l_list[0] == '\n': # skip blank lines
                continue
            for entry in l_list:
                raw_str = entry.strip()
                t_str = raw_str.lower()  # case insenstive
                if t_str == 'file':
                    check_idx = 0  # parsing title
                elif t_str == 'motifs':
                    check_idx = 1
                elif t_str == 'query':
                    check_idx = 2
                    queries.append([])
                    query_idx += 1
                elif t_str == 'thresholds':
                    if generate:
                        check_idx = 3
                elif check_idx == 2:
                    queries[query_idx].append(raw_str)
                elif check_idx == 0:
                    f_name += raw_str
                elif check_idx == 3:
                    t_entries = raw_str.split(',')
                    if len(t_entries) != 6:
                        thresholds == [8, 0, 0.2, 0.8, 1.0, 6]
                        print('MAIN: using default thresholds %s' % thresholds)
                        continue
                    for i in range(len(t_entries)):
                        thresh = t_entries[i].replace(' ','')
                        if i == 0 or i == 5:
                            thresholds.append(int(thresh))
                        else:
                            thresholds.append(float(thresh))
                elif check_idx == 1:
                    if motif_files and t_str == 'generate':
                        print('MAIN: Will generate motifs...')
                        generate = True  # for parsing reference, frame sequences
                        for frame_file in motif_files:
                            frames.append(read_frames(frame_file))
                    else:
                        print('MAIN: Will use provided motifs...')
                        d_entries = raw_str.split(',')
                        for mot in d_entries:
                            mot = mot.replace(' ', '')
                            motifs.append(mot)
    return f_name, motifs, queries, frames, thresholds


def read_frames(filename):
    """
    Opens a second text file to read in protein FASTAs for motif generation.
    This file should consist of JUST protein FASTAs.
    # may be used for comments.
    Example FASTA:
    >CRY1_Drosophila_pseudoobscura
    MVPRGANVLWFRHGLRLHDNPALLAALEEKDQGIPLIPVFIFDGESAGTKSVGYNRMRFL
    LDSLQDLDEQLQSATEGRGRLFVFEGEPTLIFRRLHEQVRLHKICAELDCEPIWNERDES
    ARLLCRELGIEYVEKVSHTLWDPRLVIETNGGIPPLTYQMFLHTVQIIGVPPRPAIDAHI
    NDATFIQLAPELRQHLGCFDQVPNPEHFNIYSDNMGFLAKINWRGGETQALALLEERLKV
    ERNAFERGYYLPNQANPNIQEAPKSMSAHLRFGCLSVRRFYWSVHDLFENVQLAACVRGV
    QIEGGAHITGQLIWREYFYTMSVNNPNYDRMEGNEICLTIPWAKPDENLLQRWRLGQTGF
    PLIDGAMRQLLAEGWLHHTLRNTVATFLTRGGLWQSWEPGLKHFLKYLLDADWSVCAGNW
    MWVSSSAFERLLDSSLVTCPVALAKRLDPEGVYIRRYVPELKNLPKEYIHEPWRLSAEQQ
    VKFECLIGVHYPERIIDLSKAVKRNMMAMTALRNSLITPPPHCRPSNEEEVRQFFWLANY
    :param filename: filename of FASTAs test file (string)
    :return: protein sequences (list of strings)
    """
    frames = []
    with open(filename, 'r') as infile:
        for line in infile.readlines():
            if line.startswith('#'):
                continue
            if line[0] == '\n':  # skip blank lines
                continue
            raw_str = line.strip()
            if '>' in raw_str:
                frames.append('')
            else:
                frames[-1] += raw_str
    return frames


def mf_main(query_file, motif_files=None):
    """
    Main function that coordinates MotifGenerator and Query
    :param query_file: filepath for settings, queries, etc.
    :param motif_files: array of files containing protein seqs for motif generation
    :return: MF_query.Query
    """
    f_name, motifs, queries, frames, thresholds = read_in(query_file, motif_files)
    if frames:
        motifs = []
        for frame_set in frames:
            # run motif generation, letting each frame be reference once
            generator = MotifGenerator(frame_set, thresholds[:-2], thresholds[-2:])
            motifs += generator.run()
    if not motifs:
        force_exit(message='ERROR: NO MOTIFS FOUND! Try adjusting generator thresholds or input motifs')
    # read in blosum and flex multiplicity
    blosum = json.loads(pkg_resources.read_text(static, 'blosum.json'))
    flex_multiplicity = json.loads(pkg_resources.read_text(static, 'flex.json'))
    return Query(f_name, motifs, queries, blosum, flex_multiplicity)
