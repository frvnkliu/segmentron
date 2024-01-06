#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 11:43:51 2023

@author: velaga01

This code preprocess and input DNA for segmentation by SERGIO. 
IWhen done it will take the DNA to be segmented and return:
    * A list of groups of repeats >20bp, i.e., a list of lists corresponding to 
    a group of repeats. E.g. [[AAAAAAAAAAAAAAAAAAAA, AAAAAcAAAAAAAAAAAAAA],
                            GGGGGGGGGGGGGGGGGGG, GGGGGtGGGGGGGGGGGGG]]. 
    The list is actually made with coordinates, i.e., [[(1, 202), (3462, 3503)], 
                                           [(4623, 4664), (9594, 9793)]
                                           
    The list also contains the strand of each repeat, i.e., 
    [[(1, 202, +), (3462, 3503, +)], [(4623, 4664, -), (9594, 9793, +)]
    
    *A list of coordinates for mini-repeats 8bp-19bp, these are also grouped,
    and have the strand information. 
    
"""
# Importing necessary modules from the Biopython package: 
# NcbiblastnCommandline for setting up and running BLAST searches, 
# and NCBIXML for parsing BLAST results.
import subprocess
from Bio.Blast import NCBIXML


#----------------Output--------------------------------

def return_repeats(sequence):
    with open("temp.fasta", "w") as f:
        f.write(">seq\n" + sequence)
    # ADD QUICK EXPLENATION FOR EACH PARAMETER
    # Run the BLAST search. The results are stored in the "blast_results.xml" file.
    blast_command = [
        'blastn',                 # Command for nucleotide-nucleotide BLAST
        '-query', 'temp.fasta',   # Input query file in FASTA format
        '-subject', 'temp.fasta', # Subject file to BLAST against
        '-out', 'blast_results.xml', # Output file in XML format
        '-outfmt', '5',           # Output format (5 for XML)
        '-gapopen', '5',          # Cost to open a gap
        '-gapextend', '5',        # Cost to extend a gap
        '-perc_identity', '80',   # Percentage identity threshold
        '-strand', 'both',        # Search both strands
        '-word_size', '20'        # Word size for initial matches
    ]
    # Execute the command using subprocess
    result = subprocess.run(blast_command, capture_output=True, text=True)

    # Check if the command was successful
    if result.returncode == 0:
        print("BLAST search completed successfully.")
    else:
        print("Error in BLAST search:")
        print(result.stderr)
        
    #Open output file
    result_handle = open("blast_results.xml")
    #Read output file
    blast_record = NCBIXML.read(result_handle)
    #List of matches found
    matches = []
    #Iterate over each alignment found
    for alignment in blast_record.alignments:
        for i, hsp in enumerate(alignment.hsps):
            #i>0 removes the sequence aligning to itself.
                if i > 0:
                    #For each alignment made there is a query and a subject, 
                    #the query is the sequence used to find homology, and the 
                    #subject is the found homology. They are both considered 
                    #repeats, so if you print the results you will see that the 
                    #subject/query pairs will be repeated with inverse orders 
                    #in the next alignment. So, we can retreive the coordinate
                    #of either all queries or all subject and get the same list 
                    #of matches. I chose to use the querries:
                    query_start = hsp.query_start
                    query_end = hsp.query_end
                    
                    if (query_start, query_end) not in matches:
                        matches.append((query_start, query_end))
                    
    result_handle.close()
    #Consolidate overlapping matches
    consolidated_matches = []
    matches.sort()
    for start, end in matches:
        if consolidated_matches and consolidated_matches[-1][1] >= start - 1:
            consolidated_matches[-1] = (consolidated_matches[-1][0], max(end, consolidated_matches[-1][1]))
        else:
            consolidated_matches.append((start, end))
    return(consolidated_matches)