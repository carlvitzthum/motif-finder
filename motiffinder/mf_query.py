# Carl Vitzthum
# For honors thesis under Andrea Tilden, spring 2016
# Class code for finding motifs within individual queries
# Look to DF_main.py for use
from __future__ import (
    absolute_import,
    print_function,
    unicode_literals
)
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Emboss.Applications import NeedleCommandline
import re
import csv
import operator
import itertools
import subprocess
from tempfile import NamedTemporaryFile


class Query(object):
    def __init__(self, f_name, motifs, queries, b, fm):
        """
        Initialization function. In addition to filling fields, runs the motif finding program through run()
        :param f_name: output filename
        :param motifs: motifs to be searched for, either from input or motif generator
        :param queries: each query sequence. Can be protein, annotated DNA, or unannotated DNA
        :param b: input blosum matrix (dictionary of dictionaries)
        :param fm: input flex multiplicity (dictionary of dictionaries)
        :return: N/A
        """
        self.f_name = f_name
        self.motifs = motifs
        self.queries = queries
        self.blosum = b
        self.flex_multiplicity = fm
        self.flex_dict = {} # created so flex will only have to run once
        self.whole_write_list = [['Query', 'Whole model', 'Identity (%)','Similarity (%)', 'Motifs only',
                                  'Concatenated motifs', 'Identity (%)','Similarity (%)', 'Concatenated non-motifs',
                                  'Identity (%)','Similarity (%)','Motifs found',
                                  'IE score','SE score']]
        # IE stands for identity enrichment, SE stands for similarity enrichment
        self.exon_write_list = [[]] # each list within the list corresponds to an exon number
        self.mot_write_dict = {}
        self.ref_query = 0
        self.ref_name = ''
        # adjust these to fine-tune needle alignment
        self.needle_gapopen = 10
        self.needle_gapextend = 0.5
        # TODO: move run() into it's own method
        for i in range(len(self.queries)):
            self.run(self.queries[i], i)
        self.whole_write_list = self.add_whole_alignments(self.whole_write_list, self.ref_query)
        self.mot_write_dict = self.add_mot_alignments(self.mot_write_dict, self.ref_name)
        self.exon_write_list = self.add_exon_alignments(self.exon_write_list, self.ref_name)
        self.whole_tsv_write(self.whole_write_list)
        self.mot_tsv_write(self.mot_write_dict)
        self.exon_tsv_write(self.exon_write_list)
        print('QUERY: finished!')

    def run(self, input_q, idx):
        """
        Sets up and searches for motifs with an input query. Builds the models that are eventually passed to TSV
        writing function
        :param input_q: input query
        :param idx: index of input query
        :return: N/A
        """
        info, sequence, exons, options, f_score = self.build_info(input_q)
        info = self.format_queries(f_score, info, options, idx)
        if 'protein' in options:
            sequence, exons = sequence, []
            options.append('transcript') if 'transcript' not in options else None
        else:
            sequence, exons = self.build_orf(sequence, exons, options)
        if f_score not in self.flex_dict:
            self.flex_dict[f_score] = self.run_flex(self.motifs, f_score)
        mot_results, mot_map = self.find_motifs(sequence, self.flex_dict[f_score])
        orf_model = self.whole_orf_models(info, sequence, mot_results)
        self.mot_models(info, mot_map)
        if 'transcript' not in options:
            exon_coords = self.exon_to_aa_ranges(exons, options)
            self.exon_models(info, orf_model, exon_coords, mot_results)

    def build_info(self, query):
        """Parsing code along used along with read_in to read the input file and source the correct info out of it
        Using a parsing index to see what is currently being read"""
        info = ''
        sequence = ''
        exons = [] # used for exon index ranges
        options = []
        flex_score = 0 # returns flex score for the query. Default is 0
        check_idx = 0  # 0 is query parsing, 1 is sequence parsing,
        #              2 is exon range parsing, 3 is option parsing
        for entry in query:
            raw_str = entry.strip()
            t_str = raw_str.lower()
            if t_str == 'sequence':
                check_idx = 1
            elif t_str == 'exons':
                check_idx = 2
            elif t_str == 'options':
                check_idx = 3
            elif check_idx == 0:
                info += raw_str
            elif check_idx == 1:
                sequence += raw_str.replace(' ', '')
            elif check_idx == 2:
                e_entries = raw_str.split(',')
                for e_range in e_entries:
                    e_split = e_range.strip().split('-')
                    if len(e_split) > 1:
                        if self.is_num(e_split[0]) == True and self.is_num(e_split[0]) == True:
                            exons.append([int(e_split[0]), int(e_split[1])])
                    else:
                        if self.is_num(e_split[0]) == True:
                            exons.append([int(e_split[0])])
            elif check_idx == 3:
                # handle options case-insenstive (with `t_str`)
                temp_o = t_str.strip().split(',')
                for opt in temp_o:
                    options.append(opt.strip())
                    if 'flex' in options[-1]:
                        flex_score = int(options[-1][-1])
        return info, sequence, exons, options, flex_score

    def format_queries(self, f_score, info, options, idx):
        """
        Returns a formatted title for each query. Adds '_fx' suffix if flex is used on this motif.
        Adds 'REF_' prefix if the protein is the reference sequence
        :param f_score: flex score (int)
        :param info: previous auery title
        :param options: list of options used in the query (strings)
        :param idx: index of the query
        :return: the formatted info string
        """
        if f_score > 0:
            info = info +'_fx'+str(f_score) # add 'fx + flex score" to query name if flex > 0
        if 'reference' in options: # see if this is the reference query (reference in options)
            self.ref_query = (idx+1)
            info = 'REF_' + info # append REF to query if reference in options
            self.ref_name = info
        return info

    def build_orf(self,in_seq,e_ranges, options):
        """
        Uses the given sequence, exon ranges (if entered) and options (to see if 'ranges' or 'reverse' are desired) to
        create a protein model and return component AA sequences corresponding to exons
        Runs 3 times to handle frameshift possibilities. Will use the frameshift that produces the longest ORF
        :param in_seq: nucleotide sequence (string)
        :param e_ranges: exon ranges, in form [start,stop] for each exon (list of lists of ints)
        :param options: user specified options for the query (list of strings)
        :return:
        """
        coords = []
        trans_seq = ''
        frame_shift_sort = []
        idv_exons = []
        if 'transcript' in options:
            trans_seq, coords = self.extract_exons(in_seq.upper())
        elif 'ranges' in options:
            for i in range(len(e_ranges)):
                temp = ''
                coords.append([len(trans_seq)])
                if len(e_ranges[i]) == 1:
                    temp = in_seq[e_ranges[i][0]-1].upper()
                else:
                    temp = in_seq[(e_ranges[i][0]-1):e_ranges[i][1]].upper()
                trans_seq += temp
                coords[i].append(len(trans_seq))
        else:
            trans_seq, coords = self.extract_exons(in_seq)
        for x in [0,1,2]: # frame shift by shifting all NTs x places to left. Sort for longest resulting ORF and use its
            transcript = self.trim_to_triplet(trans_seq[x:]) # cut off any hanging end nucleotides from 3' end
            translation = (Seq(transcript, IUPAC.unambiguous_dna).reverse_complement().translate()
                           if 'reverse' in options else Seq(transcript, IUPAC.unambiguous_dna).translate())
            orf = self.find_best_orf(translation)[0][2]
            frame_shift_sort.append([orf,len(orf),x])
        frame_shift_sort = sorted(frame_shift_sort, key = operator.itemgetter(1), reverse=True)
        used_orf = frame_shift_sort[0][0]
        used_shift = frame_shift_sort[0][2]
        for i in range(len(coords)): # add the resulting NT exons from the frameshift
            start = coords[i][0] + used_shift
            stop = coords[i][1] + used_shift  if (coords[i][1] + used_shift) <= len(trans_seq) else len(trans_seq)
            idv_exons.append(trans_seq[start:stop])
        return used_orf, idv_exons

    def extract_exons(self,in_seq):
        """Takes a given sequence string with lowercase introns and uppercase exons and returns the transcript sequence
        and a list of exon coordinates within the sequence
        """
        total_str = ''
        temp_str = ''
        cap_parse = False
        coords = []
        for i in range(len(in_seq)):
            if in_seq[i].isupper():
                if cap_parse is False:
                    cap_parse = True
                temp_str += in_seq[i]
                if i == (len(in_seq)-1): #last character
                    coords.append([len(total_str), len(total_str)+len(temp_str)])
                    total_str += temp_str
            else:
                if cap_parse is True:
                    coords.append([len(total_str), len(total_str)+len(temp_str)])
                    total_str += temp_str
                    temp_str = ''
                    cap_parse = False
        return total_str, coords

    def find_best_orf(self, prot_seq):
        """used in annotate to find the all orfs in the given protein sequence.
        order is defined by longest ORF (starting with M) and secondarily by those ending with *
        """
        seq_str = str(prot_seq)
        results = []
        used_ranges = []
        for i in range(len(seq_str)):
            if seq_str[i] == 'M':
                already_used = False
                for m_range in used_ranges: # test to see if orf identified is encapsualted within an already existing orf
                    if i in range(m_range[0],m_range[1]):
                        already_used = True
                if already_used:
                    continue
                met_idx = i
                stop_sort = 1
                try:
                    stop_idx = seq_str.index('*', i)
                except ValueError:
                    stop_idx = -1
                    stop_sort = 0
                orf_str = seq_str[met_idx:stop_idx]
                used_ranges.append([met_idx,stop_idx])
                results.append([len(orf_str),stop_sort,orf_str])
        # sort by length of ORF first and then if it contains a stop codon;
        # two equal length ORFs -> the one with a stop codon at end will be preferred
        results = sorted(results, key = operator.itemgetter(0,1), reverse = True)
        return results if len(results) > 0 else [[0,0,'']]

    def run_flex(self, mots, score):
        """Confusing function...
        In essence, removes all substiutions ([..]'s) and repetitions ({..}'s) from each input motif and expands
        it based on the input flex score. The substutions allowed based on flex score are from the BLOSUM62 matrixs,
        although the threshold to allow substitutions is arbitrary. The repetitions allowed are quite arbitrary,
        although they are based on the prevalence of any given amino acid repeating in proteins in nature
        Returns a list of "flexed" motifs that are now much more flexible when put through the search algorithm
        """
        ret_mots = []
        for pre_mot in mots:
            mot = self.process_motifs(pre_mot)
            if score == 0:
                ret_mots.append([pre_mot, [mot]])
                continue
            flex_mots = []
            gaps = [] # list of combinations of tuples containing gap possibilities
            subs = [] # for substitutions, e.g. [xyz]. These are bookmarked with '%'
            reps = [] # for repetitions, e.g. {1,3}. These are bookmarked with #
            # move subs and reps to their lists
            # there should only ever be reps in subs, e.g. [ab{1,3}c], never the other way around
            t_sub, t_rep = '', ''
            mot = list(mot)
            del_list = []  #indices to be deleted
            for i in range(len(mot)):
                if mot[i] == '[':
                    t_sub += mot[i]
                    mot[i] = '%'
                elif mot[i] == ']':
                    t_sub += mot[i]
                    del_list.append(i)
                    subs.append(t_sub)
                    t_sub = ''
                elif mot[i] == '{':
                    if mot[i-1] == ']':
                        t_rep += subs[-1]
                    else:
                        t_rep += mot[i - 1]
                    t_rep += mot[i]
                    mot[i] = '#'
                elif mot[i] == '}':
                    t_rep += mot[i]
                    del_list.append(i)
                    reps.append(t_rep)
                    t_rep = ''
                elif len(t_rep) > 0:
                    t_rep += mot[i]
                    del_list.append(i)
                elif len(t_sub) >0:
                    t_sub += mot[i]
                    del_list.append(i)
            for i in reversed(del_list):
                del mot[i]
            mot = ''.join(mot)
            # run BLOSUM subsitution on each substitution in the motif
            for i in range(len(subs)):
                l_sub = list(subs[i])
                for j in range(len(l_sub)):
                    l_sub[j] = self.blosum_substitute(l_sub[j], score, False)
                l_sub = ''.join(l_sub)
                subs[i] = ''.join(sorted(set(l_sub), key=l_sub.index))
            # run BLOSUM substitution on each case of repetition
            for i in range(len(reps)):
                l_rep = list(reps[i])
                for j in range(len(l_rep)):
                    if l_rep[j] in self.blosum:
                        l_rep[j] = self.blosum_substitute(l_rep[j], score, True) if l_rep[1] == '{' else self.blosum_substitute(l_rep[j], score, False)
                l_rep = ''.join(l_rep)
                reps[i] = ''.join(sorted(set(l_rep), key=l_rep.index))
            len_mot = len(mot)
            n_gaps = int(round(len_mot * (0.04 * (score - 1) + 0.02)))
            if n_gaps > 0:
                gaps = list(itertools.combinations(range(len_mot + n_gaps), n_gaps))
            gaps.append(('x')) if n_gaps == 0 else None # case of 0 gaps
            for opt in gaps: # iterate through all possible combinations of gap options
                t_mot = list(mot)
                for i in opt:
                    sub_count = 0
                    rep_count = 0
                    if i != 'x':
                        # .? is the regex method for a wildcard insertion that doesn't have to be present
                        t_mot.insert(i, '.?')
                    for j in range(len(t_mot)):
                        t_char = t_mot[j]
                        t_mot[j] = self.blosum_substitute(t_char, score) #BLOSUM substitution
                        if t_char in self.flex_multiplicity: # deal with multiplicity
                            t_mot[j] += ((t_char + '?') * self.flex_multiplicity[t_char][score])
                        if t_char == '%': # add back in the substitutions removed at the beginning
                            t_mot[j] = subs[sub_count]
                            sub_count += 1
                        elif t_char == "#": # add the multiplicity/repetition back in
                            t_mot[j] = reps[rep_count]
                            rep_count += 1
                flex_mots.append(''.join(t_mot))
            ret_mots.append([pre_mot, flex_mots])
        return ret_mots

    def process_motifs(self, mot):
        """Simple motif processing function. Replaces all cases of 'x' with '.'. Returns processed motif as string.
        """
        mot = mot.upper().replace('X', '.')
        mot = mot.replace('-', ',')
        return mot

    def find_motifs(self, seq, r_mots):
        """Using regular expressions modules, finds all motifs within input list (mots) in the given sequence string (seq)
        Returns a list of indices of the found motifs. Each is start,stop.
        Also returns mot_map, a list of lists with first element being the original motif name and the second being the
            actual corresponding motif found (after flex) from the sequence
        """
        results = []
        mot_map = []
        # Look for each motif in the sequence. Save the starting and ending indices in a list
        for r_mot in r_mots:
            mot = r_mot[0]
            p_mot = r_mot[1]
            for i in range(len(p_mot)):
                for match in re.finditer(p_mot[i], seq):
                    start = match.start()
                    end = match.end()
                    results.append([start, end, [mot]]) # in terms of coordinates in seq
        results = self.consolidate_motifs(results)
        # TODO: can results and mot_map be combined?
        for i in range(len(results)):
            mot_map.append([results[i][2], seq[results[i][0]:results[i][1]]])
            results[i] = [results[i][0],results[i][1]]
        results = sorted(results, key = operator.itemgetter(0))
        return results, mot_map

    def consolidate_motifs(self, mot_coords):
        """
        Consolidates redundant motifs by combining them if they overlap
        Uses indices of starting and stopping position of each motif
        Redundant motifs are deleted. Returns filtered list of coordinates

        TODO: Fix this janky code
        """
        # test lower/upper bound
        # test if sequence in sequence
        # consolidate into better candidate
        # add placeholder in same spot in list (so indices don't get messed up)
        # remove all placeholders
        for i in range(len(mot_coords)):
            for j in range(len(mot_coords)):
                if i == j:
                    continue
                d1, d2 = mot_coords[i], mot_coords[j]
                if len(d1) == 0 or len(d2) == 0:
                    continue
                if d1[1] < d2[0] or d2[1] < d1[0]:
                    continue
                if d1[0] <= d2[0] and d1[1]:
                    if d1[1] >= d2[1]: # d1 encapsulates d2 or is equal. remove d2
                        # add original motif tag for mot_map and eventual mot_models function
                        mot_coords[i][2] = list(set(mot_coords[i][2]) | set(mot_coords[j][2]))
                        mot_coords[j] = []
                    elif d1[1] < d2[1]: # if d2 goes further than d1
                        mot_coords[i][1] = mot_coords[j][1]
                        mot_coords[i][2] = list(set(mot_coords[i][2]) | set(mot_coords[j][2]))
                        mot_coords[j] = []
                elif d1[0] > d2[0]:
                    if d1[1] >= d2[1]: # d1 begins after and ends after d2
                        mot_coords[i][0] = mot_coords[j][0]
                        mot_coords[i][2] = list(set(mot_coords[i][2]) | set(mot_coords[j][2]))
                        mot_coords[j] = []
                    elif d1[1] < d2[1]: # d2 encapsulates d1
                        mot_coords[j][2] = list(set(mot_coords[i][2]) | set(mot_coords[j][2]))
                        mot_coords[i] = []
        # TODO: is that what I meant by this?
        # mot_coords = filter(None, mot_coords)
        mot_coords = [mot for mot in mot_coords if mot]
        return mot_coords

    def whole_orf_models(self, info, orf, mot_coords):
        """Takes in the whole ORF sequence and motif coordinates and returns two models of the whole ORF and a motif model:
        1. string sequence with non-motif aa's in lowercase, motif aa's in uppercase (orf_model)
        2. string with lengths of non-motif aa's between motifs in the total orf. motifs are explicity written
        3. string of all concatenated motifs
        """
        orf_model = '' # non-motif aa's are lowercase; motif aa's are uppercase
        mot_model = '' # non-motif aa regions are representing in aa length, motifs are uppercase
        concat_mots = '' # all motifs concatenated in order
        concat_nons = orf.upper() # used to model all non-motif aa's
        mot_indices = []
        for i in range(len(mot_coords)):
            mot_indices += range(mot_coords[i][0], mot_coords[i][1])
            mot_model += '_' + str(mot_coords[i][0] - mot_coords[i - 1][1]) + '_' if i > 0 else str(mot_coords[i][0]) + '_'
            mot_model += orf[mot_coords[i][0]:mot_coords[i][1]]
            mot_model += '_' + str(len(orf) - mot_coords[i][1]) if i == (len(mot_coords) - 1) else ''
            concat_mots += orf[mot_coords[i][0]:mot_coords[i][1]]
            if len(orf) > mot_coords[i][1]:
                concat_nons = concat_nons[:mot_coords[i][0]] + (mot_coords[i][1] - mot_coords[i][0]) * '_' + concat_nons[mot_coords[i][1]:]
            else:
                concat_nons = concat_nons[:mot_coords[i][0]] + (mot_coords[i][1] - mot_coords[i][0]) * '_'
        mot_num = str(len(mot_coords))
        for i in range(len(orf)):
            orf_model += orf[i].upper() if i in mot_indices else orf[i].lower()
        concat_nons = concat_nons.replace('_','')
        self.whole_write_list.append([info, orf_model,'','', mot_model, concat_mots, '', '', concat_nons, '', '',
                                      mot_num, '', '', '', ''])
        return orf_model

    def mot_models(self, query, mot_map):
        """Takes motif name - sequence mappings (from find_motifs function) and adds them to mot_dict
        if they are not there already (motif names are key). If motif naime is already in mot_dict, add
        the sequence to the list of found motifs
        """
        for entry in mot_map:
            for mot_entry in entry[0]:
                if mot_entry in self.mot_write_dict:
                    self.mot_write_dict[mot_entry].append([query, entry[1]])
                else:
                    self.mot_write_dict[mot_entry] = [[query, entry[1]]]

    def exon_to_aa_ranges(self, idv_exons, options):
        """Finds the resulting slices in the final ORF that are contributed by each original exon. Does this by iterateively
        adding exons, translating, finding the ORFs, and then indexing based on those findings.
        Returns a list of exon to aa coordinates in the final ORF, where each entry in the list is [start, stop] with
        starting and stop aa coordinates
        """
        concat_exons = ''
        aa_exon_coords = [[0]] # the coordinates of each exon in the final AA sequence (as an ORF)
        for i in range(len(idv_exons)):
            concat_exons += idv_exons[i]
            t_str = self.trim_to_triplet(concat_exons)
            t_trans = Seq(t_str, IUPAC.unambiguous_dna).reverse_complement().translate() if 'reverse' in options else Seq(t_str, IUPAC.unambiguous_dna).translate()
            sum_ORF = self.find_best_orf(t_trans)[0][2]
            if i > 0:
                aa_exon_coords.append([aa_exon_coords[i-1][1]])
            aa_exon_coords[i].append(len(sum_ORF))
        return aa_exon_coords

    def exon_models(self, info, orf_model, exon_coords, mot_coords):
        """Takes the complete orf_model (from whole_orf_models()), exon coordinates, and motif coordinates. Forms an
        exon based analysis of the protein, listing exon number, the sequence (with motifs in uppercase), and all component
        motifs in each exon. Returns a list of lists, each of which is all the above info for an exon.
        """
        exon_list = []
        for i in range(len(exon_coords)):
            exon_seq = orf_model[exon_coords[i][0]:exon_coords[i][1]]
            exon_list.append(['Exon ' + str(i+1), exon_seq, '', '']) #[exon number, exon seq, motifs, concat motifs]
            next_seq = orf_model[exon_coords[i+1][0]:exon_coords[i+1][1]] if i < (len(exon_coords)-1) else ''
            prev_seq = orf_model[exon_coords[i-1][0]:exon_coords[i-1][1]] if i > 0 else ''
            for mot in mot_coords:
                m_seq = orf_model[mot[0]:mot[1]]
                mot_str = ''
                if m_seq in exon_seq:
                    mot_str += m_seq
                # elif cases below deal with motifs that span two exons
                elif m_seq in (exon_seq+next_seq):
                    for j in range(len(m_seq)-1):
                        j+=1
                        t_this = exon_seq[-1*j:]
                        t_next = next_seq[:len(m_seq) - j]
                        if (t_this+t_next) == m_seq:
                            mot_str += t_this + '...'
                elif m_seq in (prev_seq+exon_seq):
                    for j in range(len(m_seq)-1):
                        j+=1
                        t_prev = prev_seq[-1*j:]
                        t_this = exon_seq[:len(m_seq) - j]
                        if (t_prev+t_this) == m_seq:
                            mot_str += '...' + t_this
                if len(mot_str) > 0:
                    exon_list[i][2] += (mot_str + ', ')
                    exon_list[i][3] += mot_str.replace('.', '')
        for j in range(len(exon_list)):
            if j >= len(self.exon_write_list):
                self.exon_write_list.append([[info, exon_list[j][0], exon_list[j][1],'','', exon_list[j][2][:-2], exon_list[j][3],'','']])
            else:
                self.exon_write_list[j] += [[info, exon_list[j][0], exon_list[j][1],'','', exon_list[j][2][:-2], exon_list[j][3],'','']]

    def add_whole_alignments(self, write_list, ref_query):
        """Add identity and similarity scores to output TSV write list using EMBOSS needle global alignment
        The first alignment is for total ORF models. The reference query will have 100% identity and simiarity
        The second alignment is for concatendated motifs.
        """
        for i in range(len(write_list)-1):
            i += 1
            if i != ref_query:
                if len(write_list[ref_query][1]) > 0 and len(write_list[i][1]) > 0:
                    write_list[i][2:4] = self.global_align(write_list[ref_query][1], write_list[i][1])
                else:
                    write_list[i][2:4] = ['', '']
                if len(write_list[ref_query][5]) > 0 and len(write_list[i][5]) > 0:
                    write_list[i][6:8] = self.global_align(write_list[ref_query][5], write_list[i][5])
                else:
                    write_list[i][6:8] = ['', '']
                if len(write_list[ref_query][8]) > 0 and len(write_list[i][8]) > 0:
                    write_list[i][9:11] = self.global_align(write_list[ref_query][8], write_list[i][8])
                else:
                    write_list[i][9:11] = ['', '']
            else:
                if len(write_list[ref_query][1]) > 0:
                    write_list[i][2:4] = ['100','100']
                else:
                    write_list[i][2:4] = ['-','-']
                if len(write_list[ref_query][5]) > 0:
                    write_list[i][6:8] = ['100', '100']
                else:
                    write_list[i][6:8] = ['-','-']
                if len(write_list[ref_query][8]) > 0:
                    write_list[i][9:11] = ['100','100']
                else:
                    write_list[i][9:11] = ['-', '-']
            # compute enrichment calculations for identity and similarity
            if self.is_float(write_list[i][2]) and self.is_float(write_list[i][6]):
                i_enrich = float(write_list[i][6])/float(write_list[i][2])
                s_enrich = float(write_list[i][7])/float(write_list[i][3])
                d = int(write_list[i][11])
                write_list[i][12:14] = [str(round(d*i_enrich)),str(round(d*s_enrich))]
            else:
                write_list[i][12:14] = ['0','0'];
        return write_list

    def add_mot_alignments(self, mot_dict, ref_name):
        """Add EMBOSS needle global identity and similarity alignment scores for each motif. Each motif is aligned
        against the motif from the reference sequence; if the reference sequence does not contain the given motif,
        'NO_REF' is written. Identity and similarity scores for the reference seq will always be 100.
        IF THERE ARE MULTIPLE REF SEQ ENTRIES for one motif (high possibility if flex is used):
        the alignments to EACH ref seq will be given, separated by slashes ('/') in the order the ref_seqs appear
        """
        for key in mot_dict:
            ref_seq = []
            for entry in mot_dict[key]:
                if entry[0] == ref_name:
                    ref_seq.append(entry[1])
            for entry in mot_dict[key]:
                if len(ref_seq) > 0:
                    if len(entry[1]) > 0:
                        if entry[0] == ref_name:
                            entry += ['100','100']
                        else:
                            id_str = ''
                            sim_str = ''
                            for r_seq in ref_seq:
                                [t_i, t_s] = self.global_align(entry[1],r_seq)
                                id_str += t_i + '/'
                                sim_str += t_s + '/'
                            entry += [id_str[:-1], sim_str[:-1]]
                    else:
                        entry += ['','']
                else:
                    entry += ['NO_REF', 'NO_REF']
        return mot_dict

    def add_exon_alignments(self, write_list, ref_name):
        """Given exon TSV write list, adds identity and similarity scores from EMBOSS needle global alignments
        First alignment is between reference whole exon model and whole exon model from each other query, for each exon
        Second alignment is between reference exon concatenated motifs and concat motifs from all other queries of same
        exon.
        If aligning to self, identity and similarity scores will be 100. If no alignment is possibe, - will be displayed
        If the reference query does not contain the given exon, NO_REF will be listed
        """
        for h in range(len(write_list)):
            ref_found = False
            ref_idx = 0
            for i in range(len(write_list[h])):
                if write_list[h][i][0] == ref_name:
                    ref_found = True
                    ref_idx = i
            for i in range(len(write_list[h])):
                if ref_found:
                    if write_list[h][i][0] != ref_name:
                        if len(write_list[h][i][2]) > 0 and len(write_list[h][ref_idx][2]) > 0:
                            write_list[h][i][3:5] = self.global_align(write_list[h][i][2],write_list[h][ref_idx][2])
                        else:
                           write_list[h][i][3:5] = ['', '']
                        if len(write_list[h][i][6]) > 0 and len(write_list[h][ref_idx][6]) > 0:
                            write_list[h][i][7:] = self.global_align(write_list[h][i][6],write_list[h][ref_idx][6])
                        else:
                           write_list[h][i][7:] = ['', '']
                    else:
                        if len(write_list[h][i][2]) > 0 and len(write_list[h][ref_idx][2]) > 0:
                            write_list[h][i][3] = '100'
                            write_list[h][i][4] = '100'
                        else:
                            write_list[h][i][3:5] = ['', '']
                        if len(write_list[h][i][6]) > 0 and len(write_list[h][ref_idx][6]) > 0:
                            write_list[h][i][7] = '100'
                            write_list[h][i][8] = '100'
                        else:
                            write_list[h][i][7:] = ['', '']
                else:
                    write_list[h][i][3:5] = ['NO_REF', 'NO_REF']
                    write_list[h][i][7:] = ['NO_REF', 'NO_REF']
        return write_list

    def whole_tsv_write(self, contents):
        """Writes output TSV for the whole ORF model for each query. Includes a model where selected motifs are in caps,
        with a list of specific motifs found and then a concatenated sequence for motifs. Also gives identity and similarity
        score for each whole model and concatenated motif model
        """
        with open('WHOLE_' + self.f_name, 'w') as tsv:
            for entry in contents:
                write_str = ''
                for element in entry:
                    write_str += element + '\t'
                tsv.write(write_str + '\n')

    def mot_tsv_write(self, mot_dict):
        """Writes output TSV for a motif based analysis of each query. Lists the query and found motif in that query,
        under the header of the generic motif listed in the input file. Also writes identity and similarity scores
        """
        with open('MOTS_' + self.f_name, 'w') as tsv:
            tsv.write('Query\tMotif\tIdentity (%)\tSimilarity (%)\n')
            for key in mot_dict:
                write_str = ''
                write_str += 'Input motif\t' + key + '\n'
                for mot in mot_dict[key]:
                    write_str += mot[0] + '\t' + mot[1] + '\t' + mot[2] + '\t' + mot[3] + '\n'
                tsv.write(write_str + '\n')

    def exon_tsv_write(self,contents):
        """Writes output TSV for the exon-based analysis information
        Aside from header, each line has the query name, exon numeber, exon protein seq, any motifs, all motifs
        concatenated, and the similarity and identity scores within each concatenated motifs for whole exon and concatt
        motif sequences. All are tab separated
        """
        with open('EXONS_' + self.f_name, 'w') as tsv:
            tsv.write('Query\tExon number\tExon model\tIdentity (%)\tSimilarity (%)\tMotifs only\tConcatenated motifs\tIdentity (%)\tSimilarity (%)\n')
            for exon in contents:
                for entry in exon:
                    if len(entry[2]) > 0: # don't write if it's an empty exon (may come up with 'ranges' option)
                        write_str = ''
                        for element in entry:
                            write_str += element + '\t'
                        tsv.write(write_str + '\n')
                tsv.write('\n')

    def blosum_substitute(self, char, score, out_sub=True):
        """Using globally defined BLOSUM matrix, return all substitutions of given character that are above threshold
        BLOSUM score (dependent on flex score)
        Set in_sub to False only when you are running BLOSUM subsitution of AA's already in brackets
        """
        ret = char
        if char in self.blosum:
            # TODO: check this
            # for x in ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z']:
            for x in self.blosum[char]:
                if x != char and self.blosum[char][x] >= self.flex_params(score):
                    ret += x
        if len(ret) > len(char) and out_sub:
            ret = '[' + ret + ']'
        return ret

    def flex_params(self, score):
        """Returns a flex threshold based on input type
        type is 'blosum' or 'gap'
        Blosum is for BLOSUM threshold, gap is for number of gaps allowed on either end of sequence
        """
        if score == 1:
            return 2
        elif score <= 3:
            return 1
        elif score < 5:
            return 0
        else:
            return -1

    def trim_to_triplet(self, in_str):
        """Trims the given str to make sure the number of characters is divisible by 3. If it is not, trims characters from
        end to make it so. Returns trimmed string.
        """
        if len(in_str)%3 == 1:
                in_str = in_str[:-1]
        elif len(in_str)%3 == 2:
            in_str = in_str[:-2]
        return in_str

    def is_num(self, s):
        """Test to see if given string is an int. Return True if so.
        """
        try:
            int(s)
            return True
        except ValueError:
            return False

    def is_float(self, s):
        """Test to see if given string is a float. Return True if so.
        """
        try:
            float(s)
            return True
        except ValueError:
            return False

    def global_align(self, aseq, bseq):
        """
        Perform a global alignment using EMBOSS needle with input sequences aseq and bseq
        Creates a file needle.txt for the output
        Returns the identity and similarity scores as strings in a list
        TODO: maybe combine with MotifGenerator.global_align
        :param aseq: first protein sequence (string)
        :param bseq: second protein sequence (string)
        :return: list of identity and similarity scores
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
                needle_res = self.read_score(temp.readlines())
            else:
                print('ERROR: Non-zero return code from needle alignment (query)')
                needle_res = ['ERR','ERR']
        # scores: [identity, similarity] as strings
        return needle_res

    def read_score(self, needle_lines):
        """
        Process the lines from the needle file created by global_align method
        Reads and returns identity and similarity scores as strings
        """
        scores = ['','']  # identity score, similarity score
        for line in needle_lines:
            # needle output in bytes
            if "Identity" in line:
                scores[0] = re.search('\((.*?)\)',line).group(1).replace('%', '')
            elif "Similarity" in line:
                scores[1] = re.search('\((.*?)\)',line).group(1).replace('%', '')
        return scores
