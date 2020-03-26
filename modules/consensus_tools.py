import numpy as np
import collections

def clean_seq(seq, strip_gap_flag):
    ''' Takes in a sequence and cleans any character that is not in the standard
        20 amino acid alphabetand replaces it with a gap.
        :Arguments:
            - seq = sequence to be cleaned
            - strip_gap_flag = flag to strip gaps (length_filter) or not (consensus)
        :Returns:
            - new_seq = a string object of the cleaned-up sequence
    '''
    alphabet = set('ACDEFGHIKLMNPQRSTVWY')
    gap = '-'
    
    if strip_gap_flag:
        seq = seq.upper().replace('-','').replace('.','')
    else:
        seq = seq.upper()
    new_seq = list()
    for aa in seq:
        if aa in alphabet:
            new_seq.append(aa)
        else:
            new_seq.append(gap)
    return(''.join(new_seq))

def read_fasta(in_seq_set, strip_gap_flag):
    ''' Reads a FASTA file and returns a list of sequences and a list of ids from
        the file.
        :Arguments:
            - in_seq_set = FASTA file
            - strip_gap_flag = flag to strip gaps (length_filter) or not (consensus)
        :Returns:
            - fasta_seqs = list of sequences in the input FASTA file
            - fasta_ids = list of ids in the input FASTA file
    '''
    fasta_seqs = list()
    fasta_ids = list()
    test_seq= ''
    for line in in_seq_set:
        if line[0] == '>':
            fasta_ids.append(line.rstrip())
            if len(test_seq) > 0:
                fasta_seqs.append(clean_seq(test_seq, strip_gap_flag))
            test_seq = ''
        else:
            test_seq += line.rstrip()
    fasta_seqs.append(clean_seq(test_seq, strip_gap_flag))
    return(fasta_seqs, fasta_ids)

def calc_lengths(in_seqs):
    ''' Calculates the lengths of all sequences in a list and returns a list of
        all lengths.
        :Arguments:
            - in_seqs = alignment sequences
        :Returns:
            - lengths = list of sequence lengths
    '''
    lengths = [len(i) for i in in_seqs]
    return(lengths)

def is_fasta_aligned(in_seqs):
    ''' Tests if FASTA file is aligned by determining if all sequences have the
        same length. Returns True or False.
        :Arguments:
            - in_seqs = alignment sequences
        :Returns:
            - aligned_boolean = boolean value is True if sequences are aligned
              (all sequences have the same length), False if sequences are not
              aligned
    '''
    seq_lengths = calc_lengths(in_seqs)
    aligned_boolean = seq_lengths.count(seq_lengths[0]) == len(seq_lengths)
    return(aligned_boolean)

def marginal_frequencies(in_seqs, res):
    ''' Determines the residue frequencies at each position in the alignment.
        Returns a L x q matrix where L is # of positions in alignment and
        q is the number of characters is amino acid alphabet
        :Arguments:
            - in_seqs = alignment sequences
            - res = amino acid alphabet
        :Returns:
            - matrix = numpy array of residue frequencies at all positions
    '''
    len_seqs = len(in_seqs[0])
    matrix = np.zeros((len_seqs,len(res)))
    for i in range(len_seqs):
        Fi = collections.Counter(aa[i] for aa in in_seqs)
        for z in Fi:
            res_index = res.index(z)
            matrix[i][res_index] = Fi[z] / sum(Fi.values())
    return(matrix)
    
def seq_entropy(marginals):
    ''' Calculates the sequence entropy (a measure of conservation of the position) of all positions in the alignment.
        Note: a base 2 log is used.
        Note: gaps are not considered in sequence entropy calculation.
        Returns a vector of sequence entropy for each position in alignment.
        :Arguments:
            - marginals = matrix of residue frequencies
        :Returns:
            - entropies = list of the sequence entropy at each position in the alignment
    '''
    entropies = list()
    for i, j in enumerate(marginals):
        ent_i = 0
        for k in j[:-1]:
            if k != 0:
                ent_i -= k * np.log2(k)
        entropies.append(ent_i)
    return(entropies)
    
def make_fasta(in_ids, in_seqs):
    ''' Converts list of sequences and list of sequence IDs to a format for export
        as a FASTA file.
        :Arguments:
            - in_ids = output sequence ids
            - in_seqs = output sequences
        :Returns:
            - out_fasta = list of ids and sequences for output as FASTA file
    '''
    out_fasta = list()
    for i, (j, k) in enumerate(zip(in_ids, in_seqs)):
        if i != len(in_seqs) - 1:
            out_fasta.append(''.join([j, '\n', k, '\n']))
        else:
            out_fasta.append(''.join([j, '\n', k]))
    return(out_fasta)
