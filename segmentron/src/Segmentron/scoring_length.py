import math

#Scoring function that scores a segment's length based off of Antonio's code
#Returns a positive value instead of a negative value so smaller values are preferred
def length_score(parameters, starting_index, ending_index):
    segment_length = ending_index - starting_index
    return math.exp(segment_length / parameters["max_length"])