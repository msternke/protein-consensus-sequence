# protein-consensus-sequence

## Table of contents
* [Overview](#Overview)
* [Requirements](#Requirements)
* [Installation](#Installation)
* [length_filter.py](#length_filter.py)
* [consensus.py](#consensus.py)
* [Workflow](#Workflow)
* [References](#References)


## Overview
Protein consensus sequence design has been shown to be a successful strategy for engineering highly stable proteins that retain their biological activities. A protein consensus sequences is composed of the most frequent residue at all positions in a multiple sequence alignment (MSA) of homologous protein sequences. All that is needed to design a protein consensus sequence is an MSA for the target protein family and basic coding scripts to determine residue frequencies at all positions in the MSA. Applying preprocessing steps to a sequence set can improve sequence alignment and the resulting consensus sequence. Here we have made available a script (length_filter.py) to assist in preprocessing a sequence set by filtering sequences by sequence length prior to sequence alignment, and a script (consensus.py) to determine residue frequencies at all positions in an MSA, filter residue insertions from the MSA, and determine a consensus sequence.

## Requirements
Both scripts require Python3.6 or newer. Scripts were written using Python3.6 on a MacOSX (Unix) system.

Both scripts require the following non-standard packages:    
    numpy  
    matplotlib

These packages can be installed using pip from the command line by running:
```
pip3 install numpy matplotlib
```
For information on installing or using pip if necessary see: [pip](https://pip.pypa.io/en/stable/installing/).

## Installation
To install protein-consensus-sequence, clone or download this repository to your local computer.

Both length_filter.py and consensus.py are run from the command line. Customizable parameters are given as command line options to allow the user to use the scripts as desired. The following is a tutorial on how to use both the length_filter.py and consensus.py scripts.

NOTE: Sequence sets and MSAs must be in the FASTA format.

NOTE: Both length_filter.py and consensus.py must be made executable before first use. This can be done from the command line by:
```
chmod +x length_filter.py
chmod +x consensus.py
```

## length_filter.py
The script length_filter.py will take in a sequence set (must be in FASTA format) and remove sequence truncations (very short sequences) and anamolously-long sequences (often from large insertions) from the set. This filtering is done by removing sequences that deviate from the median sequence length of the set by +/- a threshold length percentage. A default threshold of 30% deviation is used, but users can change this using the '-t' flag. This processesing step helps create better sequence alignments.

### Steps to run length_filter.py
1. Open command line interface and enter protein-consensus-sequence directory by:
```
cd <filepath_to_protein-consensus-sequence>/protein-consensus-sequence
```
2. If not does so already, change file permissions to executable (see Installation section)
3. Place desired sequence set in Inputs folder
4. Run length_filter.py script (example is for default parameters) from command line by:
```
./length_filter.py -i Inputs/<sequence_set_file_name>
```
5. A curated sequence set file, a histogram of sequence lengths, and a summary file will be saved in the Outputs folder

To use as an example, a sequence set of phosphoglycerate kinase (PGK) obtained from the Interpro database is in the Inputs folder. To run this example, run the following from the command line:
```
./length_filter.py -i Inputs/pgk_IPR015824.txt
```

### Command line options
Customizable options can be given as flags on the command line. To see the various available flag options run
```
./length_filter.py -h
```
The following table gives an overview of the flag options as well.

| Flag     |   Description   |  
|:-:|:-:| 
| -i | Filename for FASTA alignment.|
| -o | Output FASTA filename. If not given will use name of input FASTA file as template to name output files. |
| -t | Sequence length filtering threshold value (default: 0.3 for removing sequences that deviate +/- 30% from median sequence length of set). Must be a value between 0 and 1.|
| -f | Include flag to prevent saving image of sequence lengths histogram |
| -a | Include flag to keep gaps in sequences for output FASTA alignment (not recommended for further curation). |

For example, to run the PGK sequence set with a different length filtering threshold value, run the following from the command line:
```
./length_filter.py -i Inputs/pgk_IPR015824.txt -t 0.5
```

## consensus.py
The script consensus.py will take in an MSA (must be in FASTA format), calculate the residue frequencies at all positions in the MSA, filter insertion positions (see note below for explanation insertion filtering methods), and determine a consensus sequence for the MSA. 

NOTE: consensus.py can be run directly on an un-curated MSA obtained from a database such as Pfam (i.e. skipping steps 2-4 in Workflow section below). However, it is recommended that these curation steps be applied.

### Steps to run consensus.py

1. Open command line interface and enter protein-consensus-sequence directory by:
```
cd <filepath_to_protein-consensus-sequence>/protein-consensus-sequence
```
2. If not does so already, change file permissions to executable (see Installation section)
3. Place desired MSA in Inputs folder
4. Run consensus.py script (example is for default parameters) from command line by:
```
./consensus.py -i Inputs/<MSA_file_name>
```
5. A insertion filtered (gap stripped) alignment of all sequences in the MSA, a CSV file of residue frequencies at all positions in the insertion filtered alignment, a plot of the gap frequencies for all positions, a plot of sequence entropties (a measure of position conservation), and a summary file will be saved in the Outputs folder.

To use as an example, a previously curated MSA of homeodomains (HD) obtained from the Pfam database is in the Inputs folder. To run this example (using default parameters), run the following from the command line:
```
./consensus.py -i Inputs/homeodomain_PF00046_curated.txt
```

### Command line options

Customizable options can be given as flags on the command line. To see the various available flag options run
```
./consensus.py -h
```
The following table gives an overview of the flag options as well.

| Flag |   Description   |  
|:-:|:-:| 
| -i | Filename for FASTA alignment.|
|-c | Desired method for removing insertions. See below for explanation of methods |
| -o | Output FASTA filename. If not given will use name of input FASTA file as template to name output files. |
| -t | Gap frequecy threshold to define a consensus positions. Only valid for Option 1 for removing insertions. Must be a value between 0 and 1 (default: 0.5).|
| -f | Include flag to prevent saving images of MSA data analysis. |

### Insertion filtering options
Filtering an MSA for insertions is done differently by different groups. Filtering insertions can be viewed as "Which positions do I want to include in my consensus sequence?" The consensus.py script allows for users to choose which method they prefer. The available options are:   

1. (Recommended) Remove positions for which the gap frequency is > 50%. Or put conversely, keep positions for which the gap frequency is < 50%. Note that for this strategy positions included in the consensus sequence can have a gap as the most frequence character in the alignment, however more than half of the sequences in the alignment will have a residue at all included positions.
2. Remove positions for which a gap is the most frequent character. Or put conversely, keep positions for which a residue is the most frequent character.
3. Remove positions that contain gaps in a user-specified reference sequence. For this option, the user will be asked to provide the ID of the references sequence. Note that the reference sequence and ID must be present in the MSA.

## Workflow
The scripts available in protein-consensus-sequence do not encompass all steps in the process of designing a protein consensus sequence. Rather they are intended to be used in combination with a program to filter for sequence redundancy and a multiple sequence alignment program. The workflow for the entire process of determining a consensus sequence is:

1. Obtain a sequence set of target protein family from a database such as:   
    [Pfam](http://pfam.xfam.org/)  
    [Interpro](https://www.ebi.ac.uk/interpro/)  
    [SMART](http://smart.embl-heidelberg.de/)  
2. Filter sequence set by sequence lengths with length_filter.py script.
3. Filter sequence set for redundant sequences a program such as:  
    [CDHIT](http://weizhongli-lab.org/cdhit_suite/cgi-bin/index.cgi?cmd=cd-hit)  
    [UCLUST](https://drive5.com/usearch/manual/uclust_algo.html)
4. Align curated MSA using:  
    [MAFFT](https://mafft.cbrc.jp/alignment/server/)  
    [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/)
5. Determine consensus sequence using consensus.py script.


## References
For more information on protein consensus sequence design see:

1. [Consensus sequence design as a general strategy to create hyperstable, biologically active proteins. Sternke et al. PNAS, 2019.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6561275/)
2. [Consensus protein design. Porebski and Buckle. Protein Eng Des Sel, 2016.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4917058/)
