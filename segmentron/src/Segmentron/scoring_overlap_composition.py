import math

#Scoring function that scores a segment relative to the base pair composition of the overlap area at the tail end of the segment
#Score is minimized when the proportion of GC:AT is 50:50
#Returns a positive value instead of a negative value so smaller values are preferred
def overlap_composition_score(parameters, starting_index, ending_index):
    overlap_length = parameters["overlap"]
    gc_proportion = parameters["gc_overlap_composition"][ending_index - overlap_length] / overlap_length
    return math.exp(4 * abs(gc_proportion - 0.5)) - 1

#Preprocessing function to check every interval of length "overlap" in the sequence and store the number of G's or C's included
def GC_proportions(parameters):
    search_sequence = parameters["sequence"].upper()
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