# Carl Vitzthum
# For honors thesis under Andrea Tilden, spring 2016
# Class code for generating consensus motifs
# Look to mf_main.py for use
from __future__ import (
    absolute_import,
    print_function,
    unicode_literals
)
from .mf_utils import force_exit
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Emboss.Applications import NeedleCommandline
import re
import csv
import operator
import itertools
import subprocess
import random
from tempfile import NamedTemporaryFile


class MotifGenerator(object):
    def __init__(self, frame_seqs, FD_thresholds, CD_thresholds):
        """
        Initialization function
        :param frame_seqs: frame protein sequences (list of strings)
        :param FD_thresholds: list of integers and float thresholds used in evalutate_frame(). See that method
        :param CD_thresholds: [float, int] read from MF_main
        :return:
        """
        self.frame_seqs = frame_seqs
        self.FD_thresholds = FD_thresholds
        self.CD_thresholds = CD_thresholds
        # make [int, int]... convert proportion of # frames to integer
        self.CD_thresholds[0] = int(round(self.CD_thresholds[0] * len(frame_seqs)))
        # adjust these to fine-tune needle alignment
        self.needle_gapopen = 10
        self.needle_gapextend = 0.5

    def run(self):
        """
        Runs the motif generator. A random number of frames are sampled to make consensus motifs; this number is equal
        to the sampling threshold * total number of frames. For each frame, a pairwise alignment is performed between it
        and all other frames and possible motifs are identified. These are compiled and used to make a final list of
        consensus motifs, which are filtered for redundancy
        :return: consensus motifs (list of strings)
        """
        total_motifs = []
        # each sequence gets used as the reference sequence one
        for i in range(len(self.frame_seqs)):
            ref_seq = self.frame_seqs[i]
            remain_frames = self.frame_seqs[:i] + self.frame_seqs[i+1:]
            # n of the remaining frames are used, where n = the sampling threshold (CD_thresholds[0])
            # TODO: multiprocess so that sample is not needed?
            # TODO: Ideally this code should run the same every time
            sampled_frames = random.sample(remain_frames, self.CD_thresholds[0]-1)
            ref_motifs = []
            for f_seq in sampled_frames:
                r_init, b_init, m_init = self.global_align(ref_seq,f_seq)
                r_mots, f_mots = self.evaluate_frame(r_init, b_init, m_init, self.FD_thresholds)
                ref_motifs.append([r_mots, f_mots])
            ref_motifs = sorted(ref_motifs, key=lambda x: len(x[0]), reverse=True)
            # add all found motifs from reference sequence at front of ref_motifs
            frame_mots = []
            for rd in ref_motifs:
                frame_mots += rd[0]
            ref_motifs.insert(0, self.evaluate_ref(frame_mots, ref_seq))
            total_motifs += self.make_consensus(ref_motifs, self.CD_thresholds)
        total_motifs = self.reduce_redundancy(total_motifs)
        print('GENERATE: finished!')
        return total_motifs

    def evaluate_frame(self, seq_a, seq_b, metrics, thresholds):
        """
        Using the output from global align method, uses give threshold scores to search the two sequences and
        their alignment metrics and return motifs that pass.
        :param seq_a: reference sequence (a string)
        :param seq_b: frame sequence (a string)
        :param metrics: needle alignment metrics (a string)
        :param thresholds: a list of integers/floats. First (int) is minimum length of found motifs. Second (float) is
                           the percentage of gaps allowed in motifs (as a % of motif length). Third (float) is % of
                           '.' scores (needle metric) allowed in motifs. Fourth (float) is % of ':' scores allowed
                           in motifs. '|' scores always allowed.
        :return: a_mots: indexes of the reference sequence for each passing motif found, in [start,stop] form
                 b_mots: corresponding amino acid motif found from the frame sequence
        """
        # evaluation thresholds:
        # should I use percentages of length?
        if len(seq_a) != len(metrics):
            force_exit('ERROR: Lengths of Needle metrics and asis input not equal')
        min_len = thresholds[0]
        th_gaps = thresholds[1] # at this point, this should be 0
        th_diffs = thresholds[2]
        th_sims = thresholds[3]

        # a_mots are reference
        a_mots = [] # listing indeces of a_mots
        b_mots = [] # listing sequences of b_mots
        skip_point = 0
        for i in range(len(metrics)):
            if metrics[i] != '_' and i >= skip_point:
                f_t, f_g, f_d, f_s = 0, 0, 0, 0
                for k in range(len(metrics[i:])):
                    if metrics[i+k] == '|':
                        f_t += 1
                    elif metrics[i+k] == ':':
                        f_s += 1
                        if f_s > int(round(th_sims*f_t)): # see if sims exceed threshold
                            if f_t >= min_len:
                                # get indexes of reference seq location excluding gaps, '-'
                                a_mots.append([i - seq_a[:i].count('-'), i + k - seq_a[:i + k].count('-')])
                                b_mots.append(seq_b[i:i + k])
                                skip_point = i + k
                                break
                            else:
                                break
                        f_t += 1
                    elif metrics[i+k] == '.':

                        f_d += 1
                        if f_d > int(round(th_diffs*f_t)): # see if diffs exceed threshold
                            if f_t >= min_len:
                                a_mots.append([i - seq_a[:i].count('-'), i + k - seq_a[:i + k].count('-')])
                                b_mots.append(seq_b[i:i + k])
                                break
                            else:
                                break
                        f_t += 1
                    elif metrics[i+k] == '_':
                        f_g += 1
                        if f_g > int(round(th_gaps*f_t)): # see if gaps exceed threshold
                            if f_t >= min_len:
                                a_mots.append([i - seq_a[:i].count('-'), i + k - seq_a[:i + k].count('-')])
                                b_mots.append(seq_b[i:i + k])
                                skip_point = i + k
                                break
                            else:
                                break
                        f_t += 1
        for i in range(len(b_mots)):
            b_mots[i] = list(b_mots[i])
        return a_mots, b_mots

    def evaluate_ref(self, mots, ref_seq):
        """
        Using all mots from reference-frame motif comparison, find the corresponding motif sequences from the ref
        :param mots: list of [start_idx, stop_idx] motifs from evaluate frame
        :param ref_seq: reference sequence
        :return: [list of motifs, motifs sequences]
        """
        seqs = []
        for mot in mots:
            seqs.append(list(ref_seq[mot[0]:mot[1]]))
        return [mots, seqs]

    def make_consensus(self, l_pairs, find_thresh):
        """
        Collapses pairwise alignment evaluations to create consensus motifs. These motifs should be used as
        input motifs for DF_query input file
        :param l_pairs: output from a list created from multiple evaluate_frames runs
        :param find_thresh: List [integer number of pairwise alignments the motif must be found in to be returned,
                                 minumum length of motif found to be returned]
        :return: list of found consensus motifs
        """
        #use with output of evaluate initial: each pair is [seq_a indeces, seq_b sequence]
        n_pair = []
        for i in range(len(l_pairs)):
            if len(n_pair) == 0:
                n_pair = l_pairs[i] + [[1]*len(l_pairs[i][0])]
                continue
            for j in range(len(l_pairs[i][0])):
                matched = False
                l_rng = range(l_pairs[i][0][j][0],l_pairs[i][0][j][1])
                for k in range(len(n_pair[0])):
                    n_rng = range(n_pair[0][k][0],n_pair[0][k][1])
                    found_rng = sorted(list(set(n_rng) & set(l_rng)), key = int)
                    if len(found_rng) > 0:
                        f_pair = [min(found_rng),max(found_rng)+1]
                        f_seq = []
                        ends_n = [min(found_rng)-min(n_rng), max(found_rng)-max(n_rng)]
                        ends_l = [min(found_rng)-min(l_rng), max(found_rng)-max(l_rng)]
                        seq_n = n_pair[1][k][ends_n[0]:ends_n[1]] if ends_n[1] < 0 else n_pair[1][k][ends_n[0]:]
                        seq_l = l_pairs[i][1][j][ends_l[0]:ends_l[1]] if ends_l[1] < 0 else l_pairs[i][1][j][ends_l[0]:]
                        for x in range(len(seq_n)):
                            if seq_n[x] == seq_l[x]:
                                f_seq.append(seq_n[x].upper())
                            else:
                                f_seq.append("".join(set(seq_n[x].upper() + seq_l[x].upper())))
                        if len(f_seq) >= find_thresh[1]:
                            n_pair[0][k] = f_pair
                            n_pair[1][k] = f_seq
                            n_pair[2][k] += 1
                            matched = True
                if not matched:
                    if len(l_pairs[i][1][j]) >= find_thresh[1]:
                        n_pair[0].append([l_pairs[i][0][j][0],l_pairs[i][0][j][1]])
                        n_pair[1].append(l_pairs[i][1][j])
                        n_pair[2].append(1)
        # add multi substiutions back in
        final_mots = []
        for c in range(len(n_pair[2])):
            if n_pair[2][c] >= find_thresh[0]:
                found = ''
                for i in range(len(n_pair[1][c])):
                    if len(n_pair[1][c][i]) > 1:
                        found += '[' + ''.join(sorted(n_pair[1][c][i])) + ']'
                    else:
                        found += n_pair[1][c][i]
                final_mots.append(found)
        return list(set(final_mots))

    def reduce_redundancy(self, l):
        """
        Takes a list of strings, l, and returns a list eliminating all repeats and combining wherever possible
        :param l: list of strings
        :return: a reduced list of strings
        """
        l = list(set(l))
        for i in range(len(l)):
            for j in range(len(l)):
                l[i], l[j] = self.evaluate_overlap(l[i],l[j])
        l = list(set(l))
        return l

    def evaluate_overlap(self, s1, s2):
        """
        NOT a perfect function. Attempts to consolidate strings s1 and s2, taking into account the inidivudal
        substiutions (e.g. [ADF]) within them. Has a full containment and an overlapping ends case
        :param s1: string 1
        :param s2: string 2
        :return: 2 strings as return1, return2
        """
        if s1 == s2:
            return s1, s2
        mod_s1, s1_subs = self.extract_substitutions(s1)
        mod_s2, s2_subs = self.extract_substitutions(s2)
        s3 = ''
        if mod_s1 in mod_s2:
            overlap_l = 0
            overlap_r = len(mod_s2)
            # align from left
            while mod_s1 in mod_s2[overlap_l+1:]:
                overlap_l += 1
            # align from right
            while mod_s1 in mod_s2[:overlap_r-1]:
                overlap_r -= 1
            sub_idx1 = 0
            sub_idx2 = 0
            for c in range(len(mod_s2)):
                if mod_s2[c] != '&':
                    s3 += mod_s2[c]
                elif c > overlap_l and c < overlap_r:
                    concat_subs = "".join(set(s1_subs[sub_idx1])|set(s2_subs[sub_idx2]))
                    sub_idx1 += 1
                    sub_idx2 += 1
                    s3 += '[' + concat_subs + ']'
                else:
                    s3 += '[' + s2_subs[sub_idx2] + ']'
                    sub_idx2 += 1
        else: # just try overlap from left
            best_overlap = 0
            overlap = 1
            while len(mod_s1)-overlap > 0:
                if mod_s1[len(mod_s1)-overlap:] == mod_s2[:overlap]:
                    best_overlap = overlap
                overlap += 1
            if best_overlap < 4:
                return s1, s2
            sub_idx1 = 0
            sub_idx2 = 0
            for c in range(len(mod_s1) + len(mod_s2) - best_overlap):
                if c < len(mod_s1)-best_overlap:
                    if mod_s1[c] != '&':
                        s3 += mod_s1[c]
                    else:
                        s3 += '[' + s1_subs[sub_idx1] + ']'
                        sub_idx1 += 1
                elif c >= len(mod_s1):
                    c = c-len(mod_s1) + best_overlap
                    if mod_s2[c] != '&':
                        s3 += mod_s2[c]
                    else:
                        s3 += '[' + s2_subs[sub_idx2] + ']'
                        sub_idx2 += 1
                else:
                    if mod_s1[c] != '&':
                        s3 += mod_s1[c]
                    else:
                        concat_subs = "".join(set(s1_subs[sub_idx1])|set(s2_subs[sub_idx2]))
                        sub_idx1 += 1
                        sub_idx2 += 1
                        s3 += '[' + concat_subs + ']'
        return s3, s3

    def extract_substitutions(self, s):
        """
        Takes a string, s,  with substitutions and returns a string with & placeholders. Also returns a list of
        substitution possibilities in the order they were removed
        :param s: string
        :return: modified string, list of substitution strings
        """
        s = list(s)
        subs = []
        t_sub = ''
        found = False
        del_list = [] #indeces to be deleted
        for i in range(len(s)):
            if s[i] == '[':
                found = True
                s[i] = '&'
            elif s[i] == ']':
                del_list.append(i)
                subs.append(t_sub)
                found = False
                t_sub = ''
            elif found:
                t_sub += s[i]
                del_list.append(i)
        for i in reversed(del_list):
            del s[i]
        s = ''.join(s)
        return s, subs

    def global_align(self, aseq, bseq):
        """
        Perform a global alignment using EMBOSS needle with input string sequences aseq and bseq
        Creates a file needle.txt for the output
        TODO: maybe combine with Query.global_align
        """
        with NamedTemporaryFile(mode='w+') as temp:
            needle_cline = NeedleCommandline(asequence="asis:"+aseq,
                                             bsequence="asis:"+bseq,
                                             datafile='EBLOSUM62',
                                             gapopen=self.needle_gapopen,
                                             gapextend=self.needle_gapextend,
                                             outfile=temp.name)
            child = subprocess.Popen(str(needle_cline), shell=True, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
            child.wait()
            ret = child.returncode
            if ret == 0:
                temp.seek(0)
                needle_res = self.read_needle_out(temp.readlines())
            else:
                print('ERROR: Non-zero return code from needle alignment (generate)')
                needle_res = ('', '', '')
        return needle_res

    def read_needle_out(self, needle_lines):
        """
        Process the lines from the needle file created by global_align method
        :param needle_lines: list of lines from needle output file
        :return: 3 strings: seq_a and seq_b from alignment, and str score with spaces replaced by '_'
        """
        seq_a = ''
        metrics = ''
        seq_b = ''
        for i in range(len(needle_lines)):
            # find all reference seq lines (aseq in Needle)
            if needle_lines[i][0:5] == 'asis ' and len(needle_lines[i+1]) > 1:
                seq_a += needle_lines[i][21:-8].replace(' ', '_')
                metrics += needle_lines[i+1][21:-1].replace(' ', '_')
                seq_b += needle_lines[i+2][21:-8].replace(' ', '_')
        return (seq_a, seq_b, metrics)
