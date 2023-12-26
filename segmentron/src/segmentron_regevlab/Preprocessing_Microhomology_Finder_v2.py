#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 11:22:41 2023

@author: velaga01

This is a first attempt at writting a function to find relevant microhomologies.

It takes two inputs, the sequence, and the already found repeats of 20bp or more.

It finds all the repeats in that sequence (top strand only) of 19bp (100% identity), 
then it finds all the repeats of 18bp and deletes all groups of repeats of this size
that are contained within the repeats of 19bp. This is done so that we are not
counting substrings of the 19bp as indepent repeats. This process is reapeated
with smaller window sizes, with repeates of the new samller size always being
trimmed from the list of bigger microhomologies (as explained above). 

When all the microrepeats of different sizes are found, all groups that are fully
contained within the repeats of 20bp or more are deleted.

Finally, the list is then further trimmed: 
if there is no gaps between the coordinates of different instances of a repeat,
the whole repeat group is deleted. 


"""


#This function takes a sequence and a window size and returns all the
#subsequences of that window size that appear more than once (only top strand).
#It also saves the number of appearances.
def find_repeated_subsequences(input_string, window_size):
    # Dictionary to store substrings, their counts, and positions
    substring_info = {}  
    
    # Create sliding window of the specified size
    for i in range(len(input_string) - window_size + 1):
        substring = input_string[i:i + window_size]
        
        # Check if the substring already exists in the dictionary
        if substring in substring_info:
            # Increment the count
            substring_info[substring]["count"] += 1
            # Add the position of the occurrence
            substring_info[substring]["positions"].append((i, i + window_size - 1))
        else:
            # Initialize the entry for the new substring
            substring_info[substring] = {
                "count": 1,
                "positions": [(i, i + window_size - 1)]
            }
    
    # Filter out substrings that appear only once
    repeated_substrings = {substring: info for substring, info in substring_info.items() if info["count"] > 1}

    # Filter out homologies that are spaced too far apart to matter
    positions = None
    valid_substrings = []
    for substring, info in repeated_substrings.items():
        positions = info["positions"]
        minimum_gap = 0
        for i in range(info["count"] - 1):
            if ((minimum_gap == 0) or (positions[i+1][0] - positions[i][1] < minimum_gap)):
                minimum_gap = positions[i+1][0] - positions[i][1]
        if minimum_gap < 3500:
            #If all the homologies are within a 1000 nucleotide area, the minimum gap is less than 3500
            if (positions[-1][1] - positions[0][0]) > 1000:
                #Include the set of homologies if it takes up an area greater than 1000 nucleotides
                valid_substrings.append((substring, info))
    valid_repeated_substrings = {substring: info for (substring, info) in valid_substrings}
    
    
    return {window_size: valid_repeated_substrings}


#This function checks if groups of repeats are fully contained within the 
#repeats of bigger sizes. This prevents counting substrings of bigger repeats
#as independent repeats.
def eliminate_contained_coordinates(new_small_homologies, bigger_homologies):
    # Retrieve all positions from bigger_homologies for all window sizes
    bigger_positions = [pos for window_size, substrings in bigger_homologies.items() for substring_info in substrings.values() for pos in substring_info["positions"]]
    # Iterate over elements in new_small_homologies
    for window_size, small_window_substrings in list(new_small_homologies.items()):
        # List to store substrings to be removed
        to_remove = []
        #Iterate over substrings
        for small_substring, small_substring_info in small_window_substrings.items():
            # Initialize positions
            positions = small_substring_info["positions"]
            #Initialize count
            count = small_substring_info["count"]
            # List to store positions that should be retained
            retained_positions = []
            # Assume the current position is not fully contained
            retain_position = True
            # Check if any coordinate is fully contained in bigger_positions
            for start, end in positions:
                for big_start, big_end in bigger_positions:
                    if start >= big_start and end <= big_end:
                        retain_position = False
                        count-=1
                        break
                        
                # If not fully contained, retain the position
                if retain_position:
                    retained_positions.append((start, end))
            # Update the positions for the substring
            #small_substring_info["positions"] = retained_positions

            # If there are no positions remaining, mark for removal
            if count<= 0:
                to_remove.append(small_substring)
        
        # Remove substrings with no positions remaining
        for small_substring in to_remove:
            del new_small_homologies[window_size][small_substring]
    
    # Update the counts to be the number of coordinates (positions) for each substring
    for window_size, small_window_substrings in new_small_homologies.items():
        for small_substring, small_substring_info in small_window_substrings.items():
            small_substring_info["count"] = len(small_substring_info["positions"])
    
    return new_small_homologies


def eliminate_microhomologies_in_big_repeats(bigger_homologies, big_repeats):
    # Iterate over elements in bigger_homologies
    for window_size, small_window_substrings in bigger_homologies.items():
        positions = None
        bigger_positions = None
        # Assume that all substrings are fully contained in larger repeats
        all_included = True
        # Check each set of substrings
        for small_substring, small_substring_info in small_window_substrings.items():
            # Initialize positions
            positions = small_substring_info["positions"]
            # Retrieve all positions from bigger_homologies for all window sizes
            bigger_positions = big_repeats
            # Check if any coordinate is fully contained in bigger_positions
            for start, end in positions:
                for big_start, big_end in bigger_positions:
                    # If any substring is not contained within a larger repeat, the entire set is preserved
                    if not (start >= big_start and end <= big_end):
                        all_included = False
                        break
            if (all_included):
                #Delete any set of substrings that is entirely contained within larger repeats
                del bigger_homologies[window_size][small_substring]
            else:
                all_included = True
    
    # Update the counts to be the number of coordinates (positions) for each substring
    for window_size, small_window_substrings in bigger_homologies.items():
        for small_substring, small_substring_info in small_window_substrings.items():
            small_substring_info["count"] = len(small_substring_info["positions"])
    
    return bigger_homologies

def is_overlap(range1, range2):
    start1, end1 = range1
    start2, end2 = range2

    # Check if the ranges overlap or are directly adjacent
    if end1 >= start2 - 1 :
        return True
    else:
        return False

def remove_repeats_that_form_a_single_sequence(repeat_dict):
    for window_size, repeats_info in repeat_dict.items():
        repeats_to_remove = []
        
        for repeat, info in repeats_info.items():
            positions = info["positions"]
            num_positions = len(positions)
            
            # Check if all positions within a repeat overlap
            all_overlap = all(is_overlap(positions[i], positions[j]) for i in range(num_positions) for j in range(i + 1, num_positions))
            
            if all_overlap:
                repeats_to_remove.append(repeat)
        
        # Remove the repeats that need to be removed
        for repeat in repeats_to_remove:
            del repeat_dict[window_size][repeat]

    return repeat_dict  # Add this line to return the modified dictionary

def preprocess_microhomologies(input):
    input_string = input.upper()
    biggest_window_size = 19
    bigger_homologies = find_repeated_subsequences(input_string, biggest_window_size)
    #print(bigger_homologies)
    #print('----------')
    new_small_homologies = None
    for new_small_window_size in range(18,7, -1):
        new_small_homologies = find_repeated_subsequences(input_string, new_small_window_size)
        # Eliminate coordinates contained in bigger_homologies and remove elements with 0 count
        new_small_homologies = eliminate_contained_coordinates(new_small_homologies, bigger_homologies)
        #Merging new dictionary into old dictionary
        bigger_homologies[new_small_window_size]=new_small_homologies[new_small_window_size]
    #print('\n', bigger_homologies)

    bigger_homologies = eliminate_microhomologies_in_big_repeats(bigger_homologies, [(24011, 24458)])
    #print('\n', bigger_homologies)

    bigger_homologies = remove_repeats_that_form_a_single_sequence(bigger_homologies)
    #print('\n', bigger_homologies)

    return bigger_homologies                   


def write_bed_file(repeat_dict, output_file):
    with open(output_file, 'w') as bed_file:
        unique_identifier = 1
        for window_size, substrings_info in repeat_dict.items():
            for repeat_name, info in substrings_info.items():
                positions = info["positions"]
                total_positions = len(positions)
                
                for i, (start, end) in enumerate(positions, start=1):
                    # Generate the unique identifier for the repeat
                    
                    # Generate the name field using the specified format
                    name_field = f"{unique_identifier}_{i}/{total_positions}"
                    
                    # Write the BED entry to the file
                    bed_file.write(f"chromosome\t{start}\t{end+1}\t{name_field}\n")
                unique_identifier += 1