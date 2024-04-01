import math

#Scoring function that scores a segment's length based off of Antonio's code
#Returns a positive value instead of a negative value so smaller values are preferred
#Smaller scores are assigned to longer segments
def length_score(parameters, starting_index, ending_index):
    #Default value for max_length is 3000
    max_length = parameters.get("max_length", 3000)
    segment_length = ending_index - starting_index
    return math.exp(segment_length / max_length)