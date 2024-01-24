import math
from tqdm import tqdm

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