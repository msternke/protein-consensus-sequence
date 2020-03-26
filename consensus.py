#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import collections
import sys
import csv
import modules

if __name__ == '__main__':

    #set path to current directory
    path = os.getcwd() + '/'
    
    #make directory for outputs if it does not exist
    if not os.path.exists('Outputs/'):
        os.makedirs('Outputs/')

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', required=True, metavar = 'Input FASTA file.', help='Filename for FASTA alignment')
    parser.add_argument('-o', default='out.txt', metavar = 'Output gap stripped FASTA file name', help='Output FASTA filename. If not given will use name of input FASTA file as template to name output files.')
    parser.add_argument('-c', default='0', metavar = 'Method for removing insertions', help='Desired method for removing insertions. 1 = Positions with gap frequencies < threshold (0.5 default, change with -t flag). 2 =  Positions with residue as most frequent character. 3 = Positions with residues in a specific sequence. If not given will ask for user input upon running script. See README for further explantion of methods.')
    parser.add_argument('-t', type=float, default = 0.5, metavar = 'Gap frequency threshold', help='Gap frequecy threshold to define a consensus positions. Only valid for Option 1 for removing insertions. Must be a value between 0 and 1 (default: 0.5)')
    parser.add_argument('-f', action='store_true', help='Include flag to prevent saving images of MSA data analysis.')

    args = parser.parse_args()
    
    #Show help if no arguments passed
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    #FIGURE SETTINGS
    mpl.rcParams['axes.titlesize'] = 18
    mpl.rcParams['axes.labelsize'] = 18
    mpl.rcParams['xtick.labelsize'] = 14
    mpl.rcParams['ytick.labelsize'] = 14
    mpl.rcParams['axes.facecolor'] = 'FFFFFF'
    mpl.rcParams['axes.edgecolor'] = '000000'
    mpl.rcParams['axes.linewidth'] = 1.0
    mpl.rcParams['axes.labelweight'] = 'regular'
    mpl.rcParams['xtick.major.pad'] = 3
    mpl.rcParams['ytick.major.pad'] = 3
    plt.rcParams['font.family'] = 'sans-serif'
    
    print()
    print(f'Reading file: {args.i}')

    #read FASTA and populate lists for sequences and IDs
    seqs = list()
    names = list()
    
    #flag to strip gaps when cleaning sequence
    strip_gaps = False
    
    try:
        with open(args.i, 'r') as n:
            seqs, ids = modules.read_fasta(n, strip_gaps)
    except FileNotFoundError:
        print(f'Could not find file: {args.i}')
        sys.exit(1)
        
    #If no sequences were added to list, file was not in FASTA format
    if len(seqs) < 2:
        print('Provided file is not in FASTA format.')
        sys.exit(1)
        
    #Ensures sequences are aligned (all sequences have same length)
    if not modules.is_fasta_aligned(seqs):
        print('Sequences in provided FASTA are not aligned')
        sys.exit(1)

    num_seqs = len(seqs)

    #standard 20 amino acid alphabet
    res_list = list('ACDEFGHIKLMNPQRSTVWY-')

    #Matrix for residue frequencies at each position
    marginals = modules.marginal_frequencies(seqs, res_list)

    #Calculating consensus sequences for method based on user input
    consensus_positions = list()
    consensus_sequence = list()

    #Get user input for how to determine positions to include from MSA (handling gaps)
    #Users must eneter '1', '2', or '3'
    #Exits script if they do not enter valid input in five tries
    consensus_choice = args.c
    nTries = 5
    while consensus_choice not in ['1', '2', '3'] and nTries > 0:
        consensus_choice = input('\nWhich method for determining consensus positions do you want?\n***Note: Option 1 recommended.***\n1: Positions with gap frequencies < threshold (0.5 default, change with -t flag)\n2: Positions with residue as most frequent character\n3: Positions with residues in a specific sequence (you will give sequence ID)\n')
        nTries -= 1
        if nTries == 0:
            print('Invalid responses. Must enter 1, 2, or 3. Exiting program...')
            sys.exit(1)
    print()

    #Includes all positions for which gap frequency < 0.5
    #In other words, the sum of frequencies of 20 residues is > 0.5
    #This option can include positions where a gap is the most frequent occurrence
    if consensus_choice == '1':
        threshold = args.t
        if threshold < 0 or threshold > 1:
            print('Gap frequency threshold must be between 0 and 1')
            sys.exit(1)
    
        for i, j in enumerate(marginals):
            if j[-1] < threshold:
                max_res = np.argmax(j[:-1])
                consensus_positions.append(i)
                consensus_sequence.append(res_list[max_res])
                
    #Includes all positions for which a residue is the most frequent ocurrence
    #In other words, eliminates positions for which the most frequent occurrence is a gap
    elif consensus_choice == '2':
        for i, j in enumerate(marginals):
            max_res = np.argmax(j)
            if max_res != len(j) - 1:
                consensus_positions.append(i)
                consensus_sequence.append(res_list[max_res])
                
    #Includes all positions occupied by residues in a user-defined reference sequence
    else:
        #Take user input for sequence ID of reference sequence
        ref_seq = input('What is the ID of your reference sequence? ')
        
        if len(ref_seq) == 0:
            print('\nNo sequence ID was given...')
            sys.exit(1)
        
        #Users may not give sequence ID with leading '>' character
        #IDs from FASTA all have leading '>' character, so must add to reference ID
        if ref_seq[0] != '>':
            ref_seq = ''.join(['>', ref_seq])
            
        try:
            ref_index = ids.index(ref_seq)
        except ValueError:
            print('\nCould not find sequence ID in set. Check your MSA for the correct ID...')
            sys.exit(1)
        
        for i, (j, k) in enumerate(zip(seqs[ref_index], marginals)):
            if j != '-':
                consensus_positions.append(i)
                max_res = np.argmax(k[:-1])
                consensus_sequence.append(res_list[max_res])
                
    consensus_sequence = ''.join(consensus_sequence)
    print(f'Consensus sequence: {consensus_sequence}')
    print()

    #Removes positions with high gap frequencies from all sequences in alignment
    #Calculates residue frequencies for gap stripped alignment
    seqs_gap_stripped = list()
    for i in seqs:
        seqs_gap_stripped.append(''.join([i[j] for j in consensus_positions]))
    marginals_gap_stripped = modules.marginal_frequencies(seqs_gap_stripped, res_list)

    #Calculates sequence entropies for all poisitions in the gap stripped alignment
    seq_entropies = modules.seq_entropy(marginals_gap_stripped)

    #Creates matrix of residue frequencies at all positions in gap stripped alignment
    #in a format for exporting as a CSV
    out_marginals = []
    for i, j in enumerate(marginals_gap_stripped):
        if i == 0:
            out_marginals.append([''] + res_list)
            out_marginals.append([f'Position {i+1}'] + list(j))
        else:
            out_marginals.append([f'Position {i+1}'] + list(j))
            
    #Calculating gap frequencies
    #Gap is last element in each row of marginal frequencies
    gap_frequencies = [i[-1] for i in marginals_gap_stripped]

    #Setting path for file outputs
    #If the user does not supply a name for the output, use input file name
    #If user used a file from the Inputs directory, must remove 'Inputs/' from path
    out_file_name = args.o
    if out_file_name == 'out.txt':
        out_file_prefix = os.path.splitext(args.i)[0]
        if 'Inputs' in out_file_prefix:
            out_file_prefix = out_file_prefix.split('/')[1]
    else:
        out_file_prefix = out_file_name.split('.')[0]
        
    #Save figures if -f flag is not given
    if not args.f:
        #Save figure of alignment sequence entropies
        fig, ax = plt.subplots()
        ax.stem(np.arange(1, len(seq_entropies) + 1), seq_entropies, use_line_collection = True)
        ax.set_xlabel('Sequence position')
        ax.set_ylabel('Sequence entropy (bits)')
        ax.set_ylim(0, 4.5)
        fig.savefig(f'{path}Outputs/{out_file_prefix}_sequenceEntropies.png', bbox_inches = 'tight', dpi = 300)
        plt.close()
        
        #Save figure of gap stripped alignment gap frequencies
        fig, ax = plt.subplots()
        ax.stem(np.arange(1, len(gap_frequencies) + 1), gap_frequencies, use_line_collection = True)
        ax.set_xlabel('Sequence position')
        ax.set_ylabel('Gap frequency')
        ax.set_ylim(0, 1)
        fig.savefig(f'{path}Outputs/{out_file_prefix}_gapFrequencies.png', bbox_inches = 'tight', dpi = 300)
        plt.close()
                
    #Create a FASTA of the gap stripped alignment
    out_fasta = modules.make_fasta(ids, seqs_gap_stripped)

    #Save files
    with open(f'{path}Outputs/{out_file_prefix}_residueFrequencies.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(out_marginals)
        
    #Write summary of analysis to file
    with open(f'{path}Outputs/{out_file_prefix}_consensus_output.txt', 'w', newline='') as f:
        f.write(f'Determined consensus sequence for: {args.i}\n\n')
        f.write('Parameters used:\n')
        if consensus_choice ==  1:
            f.write(f'Method for removing insertions: {consensus_choice}\n')
            f.write(f'Gap frequency threshold: {threshold}\n\n')
        else:
            f.write(f'Method for removing insertions: {consensus_choice}\n\n')
        f.write(f'>Consensus_sequence\n{consensus_sequence}\n\n')
        f.write(f'MSA sequence entropy per residue: {np.mean(seq_entropies):.3f}\n\n')
        f.write(f'Wrote gap stripped alignment to: {out_file_prefix}_gapStrip.txt\n\n')
        f.write(f'Wrote CSV of residue frequencies to: {out_file_prefix}_residueFrequencies.csv\n\n')
        if args.f:
            f.write(f'Wrote plot of sequence entropies to: {out_file_prefix}_sequenceEntropies.png\n\n')
            f.write(f'Wrote plot of gap frequencies to: {out_file_prefix}_gapFrequencies.png\n\n')

    with open(f'{path}Outputs/{out_file_prefix}_gapStrip.txt', 'w') as f:
        f.writelines(out_fasta)
        
    print(f'Wrote summary of output to: Outputs/{out_file_prefix}_consensus_output.txt')
