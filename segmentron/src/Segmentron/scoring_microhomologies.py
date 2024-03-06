from tqdm import tqdm

#Scoring function that scores a segment relative to the position of microhomologies inside the segment
#Segments are penalized if the minimum distance between a microhomology and an endpoint is less than a certain constant distance
#Parameters only stores the microhomologies within this constant distance for each index
#Returns a positive value instead of a negative value so smaller values are preferred
def microhomology_score(parameters, starting_index, ending_index):
    score = 0
    contained = False
    segment_starting_index = max(0, starting_index - parameters["overlap"])
    homologies = parameters["homologies"]
    #List of checked homologies to prevent overpenalization from double counting
    checked_homologies = []
    #For each type of homology, check that each segment only contains it once
    #Penalize if any segment does contain more than one microhomology
    #First check homologies that are close to the front of the segment
    for homology in parameters["front_homologies"][starting_index]:
        contained = False
        for homology_starting_index, homology_ending_index in homologies[homology]:
            if (segment_starting_index <= homology_starting_index and homology_ending_index <= ending_index):
                if not contained:
                    contained = True
                else:
                    score = score + 225**2
        checked_homologies.append(homology)
    #Check homologies that are close to the end of the segment next
    for homology in parameters["back_homologies"][starting_index]:
        contained = False
        #Add a condition to skip homologies that are already included in front_homologies
        if homology in checked_homologies:
            continue
        for homology_starting_index, homology_ending_index in homologies[homology]:
            #Check to see if the homology is contained entirely within the segment
            if (segment_starting_index <= homology_starting_index and homology_ending_index <= ending_index):
                if not contained:
                    contained = True
                else:
                    score = score + 225**2
    return score

#Preprocessing function to store all microhomologies
#Helper function for relevant_microhomologies but can be used on its own
def microhomologies(parameters):
    search_sequence = parameters["sequence"]
    #Calculate the number of microhomologies that are in proximity to any given index
    microhomologies = preprocess_microhomologies(search_sequence, parameters["min_microhomology_length"], parameters["max_microhomology_length"], parameters["forbidden_regions"])
    #Temporary storage location for dictionary that will be stored in parameters
    homologies = {}
    #Go through all homologies and store their locations organized by sequence
    #Homologies are sorted by length and then sequence
    for homology_length, homology_list in microhomologies.items():
        for homology_sequence, homology_info in homology_list.items():
            homologies[homology_sequence] = homology_info["positions"]
    parameters["homologies"] = homologies

#Preprocessing function to store all microhomologies and their locations
def relevant_microhomologies(parameters):
    search_sequence = parameters["sequence"]
    sequence_length = len(search_sequence)
    safe_distance = parameters["microhomology_distance"]
    min_length = parameters["min_microhomology_length"]
    max_length = parameters["max_microhomology_length"]
    #Call microhomologies as a helper function
    microhomologies(parameters)
    homologies = parameters["homologies"]
    #Arrays are created with a bit of unneeded space to make indices easier to understand and read
    parameters["front_homologies"] = [None] * (sequence_length + 1)
    parameters["back_homologies"] = [None] * (sequence_length + 1)
    #Temporary storage locations for dictionaries that will be put into front_homologies and back_homologies
    current_front_homologies = []
    current_back_homologies = []
    #Calculate front and back homologies for index 0
    for homology, homology_positions in homologies.items():
        for start, end in homology_positions:
            if (0 <= start) and (start < safe_distance):
                current_front_homologies.append(homology)
                break
    parameters["front_homologies"][0] = current_front_homologies
    parameters["back_homologies"][0] = current_back_homologies
    #For each index, check through the entire dictionary of microhomologies
    tagged_for_removal = False
    homologies_to_remove = []
    subsequence = ""
    #Use a sliding window approach to finding all relevant homologies for any given index
    progress_bar = tqdm(desc = "indices preprocessed: ", total = (sequence_length) + 1)
    for index in range(1, sequence_length + 1):
        #Remove homologies relevant to index - 1 but not index
        #Handle front homologies
        for homology in current_front_homologies:
            tagged_for_removal = False
            for start, end in homologies[homology]:
                #If a homology might slip out of relevancy, check that it doesn't repeat twice in the relevant region
                if (start == index - 1):
                    tagged_for_removal = True
                if (tagged_for_removal) and (index <= start) and (start < index + safe_distance):
                    tagged_for_removal = False
                    break
            if tagged_for_removal:
                homologies_to_remove.append(homology)
        current_front_homologies = [homology for homology in current_front_homologies if homology not in homologies_to_remove]
        #Handle back homologies
        for homology in current_back_homologies:
            tagged_for_removal = False
            for start, end in homologies[homology]:
                #If a homology might slip out of relevancy, check that it doesn't repeat twice in the relevant region
                if (end == index - safe_distance):
                    tagged_for_removal = True
                if (tagged_for_removal) and (index - safe_distance < end) and (end <= index):
                    tagged_for_removal = False
                    break
            if tagged_for_removal:
                homologies_to_remove.append(homology)
        current_back_homologies = [homology for homology in current_back_homologies if homology not in homologies_to_remove]
        #Add homologies relevant to index
        #Handle start homologies
        for length in range(min_length, max_length):
            if (index + safe_distance - 1 + length > sequence_length):
                break
            subsequence = search_sequence[index + safe_distance - 1 : index + safe_distance - 1 + length]
            if (homologies.get(subsequence, 0)):
                for start, end in homologies[subsequence]:
                    if (start == index + safe_distance - 1 and subsequence not in current_front_homologies):
                        current_front_homologies.append(subsequence)
                        break
        #Handle back homologies
        for length in range(min_length, max_length):
            if (index - length < 0):
                break
            subsequence = search_sequence[index - length : index]
            if (homologies.get(subsequence, 0)):
                for start, end in homologies[subsequence]:
                    if (end == index and subsequence not in current_back_homologies):
                        current_front_homologies.append(subsequence)
                        break
        parameters["front_homologies"][index] = current_front_homologies
        parameters["back_homologies"][index] = current_back_homologies
        progress_bar.update(1)
    return 0

def preprocess_microhomologies(input, min_length, max_length, forbidden_regions):
    input_string = input.upper()
    bigger_homologies = find_repeated_subsequences(input_string, max_length)
    #print(bigger_homologies)
    #print('----------')
    new_small_homologies = None
    for new_small_window_size in range(max_length - 1, min_length - 1, -1):
        new_small_homologies = find_repeated_subsequences(input_string, new_small_window_size)
        # Eliminate coordinates contained in bigger_homologies and remove elements with 0 count
        new_small_homologies = eliminate_contained_coordinates(new_small_homologies, bigger_homologies)
        #Merging new dictionary into old dictionary
        bigger_homologies[new_small_window_size]=new_small_homologies[new_small_window_size]
    #print('\n', bigger_homologies)

    bigger_homologies = eliminate_microhomologies_in_big_repeats(bigger_homologies, forbidden_regions)
    #print('\n', bigger_homologies)

    bigger_homologies = remove_repeats_that_form_a_single_sequence(bigger_homologies)
    #print('\n', bigger_homologies)

    return bigger_homologies  

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
        for small_substring, small_substring_info in small_window_substrings.copy().items():
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
                del (bigger_homologies[window_size])[small_substring]
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