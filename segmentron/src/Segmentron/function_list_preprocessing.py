import math
from tqdm import tqdm
from . import Preprocessing_Microhomology_Finder_v2 as Preprocessing_Microhomology_Finder
from . import Preprocessing_Repeat_Finder_v2 as Preprocessing_Repeat_Finder

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
    microhomologies = Preprocessing_Microhomology_Finder.preprocess_microhomologies(search_sequence, parameters["min_microhomology_length"], parameters["max_microhomology_length"], parameters["forbidden_regions"])
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

#Preprocessing function to find forbidden regions
def forbidden_regions(parameters):
    parameters["forbidden_regions"] = Preprocessing_Repeat_Finder.return_repeats(parameters["sequence"])
    return 0