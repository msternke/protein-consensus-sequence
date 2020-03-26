#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
import modules

if __name__ == '__main__':

    #Set path to current directory
    path = os.getcwd() + '/'
    
    #make directory for outputs if it does not exist
    if not os.path.exists('Outputs/'):
        os.makedirs('Outputs/')

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', required=True, metavar = 'Input FASTA file', help='Filename for FASTA alignment.')
    parser.add_argument('-o', metavar = 'Output FASTA filename. If not given will use name of input FASTA file as template to name output files.', default='out.txt', help='Output FASTA filename')
    parser.add_argument('-t', metavar = 'Length filtering threshold', type=float, default=0.3, help='Sequence length filtering threshold value (default: 0.3 for removing sequences that deviate +/- 30% from median sequence length of set, change with -t flag). Must be a value between 0 and 1.')
    parser.add_argument('-f', action='store_true', help='Include flag to prevent saving image of sequence lengths histogram')
    parser.add_argument('-a', action='store_true', help='Include flag to keep gaps in sequences for output FASTA alignment (not recommended for further curation).')

    args = parser.parse_args()

    # show help if no arguments passed
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
    print()

    #read FASTA and populate lists for sequences and IDs
    seqs = list()
    names = list()
    
    #flag to strip gaps when cleaning sequence
    strip_gaps = True
    
    try:
        with open(args.i, 'r') as n:
            seqs, ids = modules.read_fasta(n, strip_gaps)
    except FileNotFoundError:
        print(f'Could not find file: {args.i}')
        sys.exit(1)
        
    #If no sequences were added to list, file was not in FASTA format
    if len(seqs) < 1:
        print('Provided file is not in FASTA format.')
        sys.exit(1)

    num_seqs = len(seqs)
    print(f'Number of sequences in initial sequence set: {num_seqs}')
    print()

    #Calculating lengths of all sequences in alignment
    #Note: gaps are removed before calculating sequence lengths
    lengths = modules.calc_lengths(seqs)
    
    threshold = args.t
    if threshold < 0 or threshold > 1:
        print('Sequence filtering threshold must be between 0 and 1')
        sys.exit(1)

    #Calculating statistics for length filtering
    med_length = int(np.median(lengths))
    lower_threshold = int(round(med_length - med_length * threshold))
    upper_threshold = int(round(med_length + med_length * threshold))

    print(f'Median sequence length: {med_length} residues')
    print(f'Sequence length lower boundary: {lower_threshold} residues')
    print(f'Sequence length upper boundary: {upper_threshold} residues')
    print()

    #Create FASTA of sequences that fall within threshold values
    #Gaps are remove gaps from sequences unless user indicates otherwise with -a flag
    out_fasta = list()
    if args.a:
        for i, (j, k, l) in enumerate(zip(lengths, ids, seqs)):
            if j >= lower_threshold and j <= upper_threshold:
                if i != num_seqs - 1:
                    out_fasta.append(f'{k}\n{l}\n')
                else:
                    out_fasta.append(f'{k}\n{l}')
    else:
        for i, (j, k, l) in enumerate(zip(lengths, ids, seqs)):
            if j >= lower_threshold and j <= upper_threshold:
                if i != num_seqs - 1:
                    out_fasta.append(f"{k}\n{l.replace('-','')}\n")
                else:
                    out_fasta.append(f"{k}\n{l.replace('-','')}")
                    

    num_seqs_out = len(out_fasta)
    print(f'Number of sequences in final sequence set: {num_seqs_out}')
    print()

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

    #Save figure if -f flag is not given
    if not args.f:
        fig, ax = plt.subplots()
        ax.hist(lengths, bins='doane', color='b', edgecolor='k', alpha=0.65)
        ax.set_xlabel('Sequence length')
        ax.set_ylabel('Count')
        ax.axvline(med_length - threshold*med_length, color='k', linestyle='dashed')
        ax.axvline(med_length + threshold*med_length, color='k', linestyle='dashed')
        fig.savefig(f'{path}Outputs/{out_file_prefix}_sequence_length_hist.png', bbox_inches = 'tight', dpi = 300)
        plt.close()
    
    #Save length filtered alignment
    with open(f'{path}Outputs/{out_file_prefix}_lengthFiltered.txt', 'w') as f:
        f.writelines(out_fasta)
        
    #Write summary of analysis to file
    with open(f'{path}Outputs/{out_file_prefix}_lengthFiltered_output.txt', 'w', newline='') as f:
        f.write(f'Filtered sequence set: {args.i}\n\n')
        f.write('Parameters used:\n')
        f.write(f'Sequence length filtering threshold: {threshold}\n\n')
        f.write(f'Number of sequences in initial alignment: {num_seqs}\n\n')
        f.write(f'Median sequence length: {med_length} residues\n')
        f.write(f'Sequence length lower boundary: {lower_threshold} residues\n')
        f.write(f'Sequence length upper boundary: {upper_threshold} residues\n\n')
        f.write(f'Number of sequences in final alignment: {num_seqs_out}\n\n')
        f.write(f'Wrote length filtered sequence set to: {out_file_prefix}_lengthFiltered.txt')
        
    print(f'Wrote filtered alignment to: Outputs/{out_file_prefix}_lengthFiltered.txt')
