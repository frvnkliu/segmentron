import math
from tqdm import tqdm
import Preprocessing_Microhomology_Finder_v2

#Scoring function that scores a segment's length based off of Antonio's code
#Returns a positive value instead of a negative value so smaller values are preferred
def length_score(parameters, starting_index, ending_index):
    segment_length = ending_index - starting_index
    return math.exp(segment_length / parameters["max_length"])

#Scoring function that scores a segment relative to the forbidden regions that is based off of Antonio's code
#Returns a positive value instead of a negative value so smaller values are preferred
#Needs to be able to take in a parameter for overlap length eventually but for now, the overlap length is fixed at 337
def forbidden_region_score(parameters, starting_index, ending_index):
    score = 0
    forbidden_regions = parameters["forbidden_regions"]
    #Create a large forbidden region encompassing all the relevant forbidden regions
    big_forbidden_region_starting_index = len(parameters["sequence"])
    big_forbidden_region_ending_index = 0
    region_starting_index = 0
    region_ending_index = 0
    #Adjust starting and ending indices as you go through the forbidden regions
    for forbidden_region in forbidden_regions:
        region_starting_index, region_ending_index = forbidden_region
        if region_starting_index < ending_index and region_ending_index > starting_index:
            if region_starting_index < big_forbidden_region_starting_index:
                big_forbidden_region_starting_index = region_starting_index
            if region_ending_index > big_forbidden_region_ending_index:
                big_forbidden_region_ending_index = region_ending_index
    big_forbidden_region = (big_forbidden_region_starting_index, big_forbidden_region_ending_index + parameters["overlap"])
    big_forbidden_region_length = big_forbidden_region[1] - big_forbidden_region[0]
    if (big_forbidden_region_length < 0):
        return 0
    #Score is increased when the big forbidden region overlaps with the edge points of the segment
    elif (big_forbidden_region[0] <= ending_index and big_forbidden_region[1] >= ending_index) or (big_forbidden_region[0] <= starting_index and big_forbidden_region[1] >= starting_index):
        score = 450**2
    #Score is minimized when the forbidden regions are centered inside the segment
    else:
        score = (450**2)/(abs((big_forbidden_region[0]- starting_index) * (ending_index-big_forbidden_region[1])) + 1)
    return score

#Scoring function that scores a segment relative to the base pair composition of the overlap area at the tail end of the segment
#Score is minimized when the proportion of GC:AT is 50:50
#Returns a positive value instead of a negative value so smaller values are preferred
def overlap_composition_score(parameters, starting_index, ending_index):
    overlap_length = parameters["overlap"]
    gc_proportion = parameters["gc_overlap_composition"][ending_index - overlap_length] / overlap_length
    return math.exp(4 * abs(gc_proportion - 0.5)) - 1

#Scoring function that scores a segment relative to the types of forbidden regions contained inside the segment
#Score is minimized when all forbidden regions contained inside the segment are pairwise dissimilar
#Input dictionary must contain a list of lists of forbidden regions that are similar
#Returns a positive value instead of a negative value so smaller values are preferred
def forbidden_region_class_score(parameters, starting_index, ending_index):
    forbidden_region_classes = parameters["forbidden_region_classes"]
    if len(forbidden_region_classes) == 1:
        return 0
    region_starting_index = 0
    region_ending_index = 0
    score = 0
    contained = False
    #For each class, check that each segment does not intersect with more than one forbidden region of that class
    #Penalize heavily if any segment does intersect more than one forbidden region of any class
    for forbidden_region_class in forbidden_region_classes:
        contained = False
        for forbidden_region in forbidden_region_class:
            region_starting_index, region_ending_index = forbidden_region
            #Check to see if the forbidden region intersects at all with the segment
            if not(region_ending_index <= starting_index or region_starting_index >= ending_index):
                #Allow at most one forbidden region per class to intersect with the segment before penalization
                if not contained:
                    contained = True
                else:
                    score = score + 450**2
    return score

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

#Basic accumulation function that adds together scores
def addition_function(a, b):
    return a + b

#Preprocessing function to check every interval of length "overlap" in the sequence and store the number of G's or C's included
def GC_proportions(parameters):
    search_sequence = parameters["sequence"]
    sequence_length = len(search_sequence)
    overlap = parameters["overlap"]
    #Calculate GC proportions in every interval of length "overlap"
    overlap_gc_counts = [0] * (sequence_length - overlap + 1)
    #Calculate GC proportions in the first interval of length "overlap"
    count = search_sequence[0 : overlap].count("G") + search_sequence[0 : overlap].count("C")
    overlap_gc_counts[0] = count
    #Use a sliding window approach to storing the GC count in each interval of length "overlap"
    for i in range(1, sequence_length - overlap + 1):
        if search_sequence[i] == "G" or search_sequence[i] == "C":
            count = count - 1
        if search_sequence[i + overlap - 1] == "G" or search_sequence[i + overlap - 1] == "C":
            count = count + 1
        overlap_gc_counts[i] = count
    parameters["gc_overlap_composition"] = overlap_gc_counts
    return None

#Preprocessing function to store all microhomologies
#Helper function for relevant_microhomologies but can be used on its own
def microhomologies(parameters):
    search_sequence = parameters["sequence"]
    #Calculate the number of microhomologies that are in proximity to any given index
    microhomologies = Preprocessing_Microhomology_Finder_v2.preprocess_microhomologies(search_sequence)
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
        for length in range(8, 20):
            if (index + safe_distance - 1 + length > sequence_length):
                break
            subsequence = search_sequence[index + safe_distance - 1 : index + safe_distance - 1 + length]
            if (homologies.get(subsequence, 0)):
                for start, end in homologies[subsequence]:
                    if (start == index + safe_distance - 1 and subsequence not in current_front_homologies):
                        current_front_homologies.append(subsequence)
                        break
        #Handle back homologies
        for length in range(8, 20):
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