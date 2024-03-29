import snapgene_reader as sgreader
import math
from typing import List, Callable
from tqdm import tqdm
#Multiprocessing is not compatible with web assembly
#from multiprocessing.managers import SharedMemoryManager
import time
#from tqdm.contrib.concurrent import process_map
from . import scoring_forbidden_regions
from . import scoring_length
from . import scoring_microhomologies
from . import scoring_overlap_composition
from . import scoring_accumulator
import json

class segmentron:
    #Define stored variable types
    #Preprocessing functions that will be used to collect important information before dynamic programming begins
    preprocessing_functions: List[Callable]
    #Segment scoring functions must return negative values since the class tries to maximize score
    segment_scoring_functions: List[Callable]
    #Accumulating function dictates how scores of multiple segments are combined
    accumulating_function: Callable
    #Stored parameter list containing all necessary information
    parameters: dict
    #Set to store the segmentation for future operations following every segmentation
    segmentation: List[int]
    #Set to store the score of the stored segmentation
    score = float
    #Stored arrays for usage in dynamic programming
    optimal_scores: List[float]
    optimal_cuts: List[int]
    
    
    def __init__(self, preprocessing_functions, segment_scoring_functions, accumulating_function, parameters):
        """Initializes a segmentron class with a given list of scoring functions, a function to combine these scores, and a set of parameters
        The parameters should be passed as a dictionary and all scoring functions should take a dictionary of parameters as arguments
        There are default values for parameters for max_length, min_length, and overlap length
        Parameters that are set during function calls are initially initialized as empty or 0 depending on type"""
        #Set stored variables based on initialization call
        #List of preprocessing functions to be used
        self.preprocessing_functions = preprocessing_functions
        #List of scoring functions to be used
        self.segment_scoring_functions = segment_scoring_functions
        #Dictionary of parameters to be passed to scoring functions
        self.parameters = dict(parameters)
        #Function to be used to compile multiple scores from multiple functions into a single value
        self.accumulating_function = accumulating_function
        #Ensure certain variables are included in parameters and if not included, provide default values
        #Max and min lengths are given constraint lengths for any segment with default values of 3000 and 1000 respectively
        if (self.parameters.get("max_length") is None):
            self.parameters["max_length"] = 3000
        if (self.parameters.get("min_length") is None):
            self.parameters["min_length"] = 1000
        #Overlap is the length of the desired overlap region
        if (self.parameters.get("overlap") is None):
            self.parameters["overlap"] = 100
        #Microhomology_distance is the minimum distance from the ends of the segment where a microhomology can be safely included
        if (self.parameters.get("microhomology_distance") is None):
            self.parameters["microhomology_distance"] = 20
        #Min_microhomology_length is the minimum length for a small microhomology to be recorded
        if (self.parameters.get("min_microhomology_length") is None):
            self.parameters["min_microhomology_length"] = 8
        #Max_microhomology_length is the maximum length for a small microhomology to be recorded
        if (self.parameters.get("max_microhomology_length") is None):
            self.parameters["max_microhomology_length"] = 19
        #Forbidden_region_classes is the number of types of forbidden regions to be considered
        if (self.parameters.get("forbidden_region_classes") is None):
            self.parameters["forbidden_region_class_count"] = 1
        if(self.parameters.get("block_size") is None):
            self.parameters["block_size"] = 10000
        #Set up storage areas for later use in dynamic programming
        #Set after reading from a file
        self.parameters["sequence"] = ""
        self.parameters["forbidden_regions"] = []
        #Parameters determining if forbidden_regions should be read from the file or generated
        if (self.parameters.get("forbidden_regions_from_file") is None):
            self.parameters["forbidden_regions_from_file"] = False
        if (self.parameters.get("forbidden_region_generation") is None):
            self.parameters["forbidden_region_generation"] = True
        #Parameters determining how forbidden_regions should be read from a file
        if (self.parameters.get("color") is None):
            self.parameters["color"] = "#ff0000"
        #Set to store the segmentation for future operations following every segmentation
        self.segmentation = []

    def encodeSegmentron(self, blockIndex):
        #encodes the current segmentation
        data = {
            "parameters": self.parameters,
            "segmentation": {
                "segmentation": self.segmentation,
                "optimal_scores": self.optimal_scores,
                "optimal_cuts": self.optimal_cuts,
                "blockIndex": blockIndex
            }
        }
        file_path = "segmentron.json"

        with open(file_path, 'w') as json_file:
            json.dump(self.to_dict(), json_file)
    
    #Score a segment using all stored segment scoring functions
    #Parameters input dictionary must contain all necessary parameters for all stored segment scoring functions
    def score_segment(self, parameters, starting_index, ending_index):
        """Returns a total score calculated using all stored scoring functions on the same dictionary of parameters
        Total score is obtained from using the accumulating_function to compile all the scores into a single value"""
        total_score = 0
        #Calculate score under each scoring function and combine them using the accumulating function
        for segment_scoring_function in self.segment_scoring_functions:
            total_score = self.accumulating_function(total_score, segment_scoring_function(parameters, starting_index, ending_index))
        return total_score
    
    #Segment a sequence by reading from a given filepath
    def segment_from_file(self, filepath, coarseness = 1, multiprocessing_cores = 0):
        """Given a file path to a valid snapgene file, the segmentron will read the sequence and forbidden regions from the file
        Since the given file may have multiple features, only features marked with the color red (#ff0000) will be counted as forbidden regions
        This function returns the optimal score and the corresponding optimal segmentation by calling a helper function
        This helper function will also store the sequence and forbidden regions that were read from the file"""
        self.parameters["filepath"] = filepath
        dictionary = sgreader.snapgene_file_to_dict(filepath)
        sequence = dictionary["seq"]
        self.parameters["sequence"] = sequence
        #Preprocess the sequence and store relevant information
        self.preprocess()
        forbidden_regions = self.parameters["forbidden_regions"]
        return self.segment(sequence.upper(), forbidden_regions, coarseness, multiprocessing_cores, False)

    #Segment a sequence after being given the sequence and a set of forbidden regions
    def segment(self, sequence, forbidden_regions, coarseness = 1, multiprocessing_cores = 0, preprocessing = True):
        """Given a sequence and a list of forbidden regions, the segmentron will find an optimal segmentation through dynamic programming
        This helper function is not supposed to be called unless it is through the segment_from_file function
        This helper function returns the optimal score and the corresponding optimal segmentation"""
        sequence_length = len(sequence)
        #Set up timer
        start_time = time.perf_counter()
        #Pull out sequence and other parameters from parameters list ahead of time
        #Store sequence and forbidden_regions for possible future operations
        self.parameters["sequence"] = sequence.upper()
        self.parameters["forbidden_regions"] = forbidden_regions
        #Preprocess the sequence and store relevant information if necessary
        if (preprocessing):
            self.preprocess()
        #Set up multiprocessing
        smm = None
        #Set up dynamic programming arrays for optimal scores
        """if (multiprocessing_cores > 0):
            smm = SharedMemoryManager()
            smm.start()
            #Store scores for the optimal segmentations for any subsequence starting from index 0
            self.optimal_scores = smm.ShareableList([None for i in range(sequence_length + 1)])
            #Store optimal cut for any index starting from index 0
            #For each index, only the last cut is stored to save space
            #The full segmentation can be obtained recursively
            self.optimal_cuts = smm.ShareableList([None for i in range(sequence_length + 1)])
        else:"""
        #Store scores for the optimal segmentations for any subsequence starting from index 0
        self.optimal_scores = [None for i in range(sequence_length + 1)]
        #Store optimal segmentation for any subsequence starting from index 0
        self.optimal_cuts = [None for i in range(sequence_length + 1)]
        #Set up progress bar using tqdm
        progress_bar = tqdm(desc = "indices completed: ", total = (sequence_length // coarseness) + 1)
        #Set up base case of an empty sequence
        self.optimal_scores[0] = 0
        self.optimal_cuts[0] = 0
        #Iterate over the rest of the sequence
        #If multiprocessing is used, use built in TQDM functions to create a multiprocessing pool to perform the dynamic programming
        """if (multiprocessing_cores > 0):
            #If coarseness > 1, certain indices will be skipped which will reduce runtime by a factor of about coarseness^2
            args = [current_index for current_index in range(0, sequence_length + 1, coarseness)]
            args.append(sequence_length)
            process_map(self.segment_subsequence, args, chunksize = math.floor(self.parameters["min_length"] / (multiprocessing_cores)), max_workers = multiprocessing_cores)
        #If only a single process is used, call subfunctions to perform the dynamic programming
        else:"""
        #If coarseness > 1, certain indices will be skipped which will reduce runtime by a factor of about coarseness^2
        for current_index in range(0, sequence_length + 1, coarseness):
            #For each index, check previous cases
            self.segment_subsequence(current_index)
            #Update progress bar
            progress_bar.update(1)
            if (current_index % 5000 == 0):
                print(f"A total of {current_index} indices have been processed")
        #Ensure that the entire sequence is sequenced in case it was skipped due to the coarseness setting
        self.segment_subsequence(sequence_length)
        progress_bar.update(1)
        #Shut down progress_bar and multiprocessing resources
        """if (multiprocessing_cores > 0):
            #Before shutting down the SharedMemoryManager, store the stored data into a normal list
            self.optimal_cuts = list(self.optimal_cuts)
            self.optimal_scores = list(self.optimal_scores)
            smm.shutdown
        else:"""
        progress_bar.close()
        #Store segmentation for possible future operations
        self.segmentation = []
        current_index = sequence_length
        current_cut = self.optimal_cuts[sequence_length]
        #Construct full segmentation from stored optimal cuts
        while not ((current_cut is None) or (current_cut == current_index)):
            self.segmentation.append(current_index)
            current_index = current_cut
            current_cut = self.optimal_cuts[current_index]
        #If there is no valid cut at any point, invalidate the entire segmentation
        if (current_index is None):
            self.segmentation = None
        else:
            self.segmentation.append(0)
            self.segmentation.reverse()
        print(f"This function has taken {time.perf_counter() - start_time} seconds to run")
        return self.optimal_scores[sequence_length], self.segmentation
    
    #Goes through all preprocessing functions given and collects all information necessary
    def preprocess(self):
        for preprocessing_function in self.preprocessing_functions:
            preprocessing_function(self.parameters)
        return None

    #Checks previous entries in stored arrays to find an optimal segmentation for a subsequence from 0 to a given ending_index
    def segment_subsequence(self, ending_index):
        """Given a specific index, this function checks all valid segments ending at that index
        This function stores the segment with the lowest score that is not None
        If no segment has a score that is not None, then no segment is stored and the stored score is set as None"""
        #Set up temporary storage areas
        optimal_score = None
        optimal_cut = None
        current_score = 0
        #Check through all possible starting indices that result in a segment of valid length
        for search_index in range(max(ending_index - self.parameters["max_length"], 0), ending_index - self.parameters["min_length"] + 1, 1):
            #Temporarily store the score of the segmentation corresponding to the current search_index for comparison
            current_score = self.score_segmentation(search_index, ending_index)
            #Stores the smallest score that is not None along with the corresponding index
            if not (current_score is None) and ((optimal_score is None) or (current_score < optimal_score)):
                optimal_score = current_score
                optimal_cut = search_index
        #If there is no valid cut, there is no valid segmentation at this index
        if not (optimal_cut is None):
            #Store the lowest possible score (if it exists) and the corresponding sequence of cuts
            self.optimal_scores[ending_index] = optimal_score
            self.optimal_cuts[ending_index] = optimal_cut

    #Checks specific segmentation, compares it to a stored segmentation, and stores the better segmentation
    def score_segmentation(self, starting_index, ending_index):
        """Given a specific segment from the starting_index to the ending_index, this function scores the corresponding segmentation
        Score is calculated by accumulating the score of the segment with the optimal score of the subsequence from 0 to the starting_index
        This function returns the accumulated score if the segmentation is valid and None if the segmentation is invalid
        This function is a helper function meant to be called by the _subsequence function during dynamic programming"""
        prior_score = self.optimal_scores[starting_index]
        #If the prior score of the region from 0 to starting_index is None, the segment is invalid or not possible
        if not (prior_score is None):
            #Return the score of the segment accumulated with the score of the prior subsequence
            return self.accumulating_function(self.score_segment(self.parameters, starting_index, ending_index), prior_score)
        #Return None if the segment is invalid or not possible
        return None
    
    #Takes a segmentation and returns the corresponding segments including overlap regions
    def segmentation_to_segments(self):
        """Using the stored segmentation, which is a list of cuts, this function returns the corresponding list of segments
        These segments include the regions that are meant to be overlapping"""
        segments = []
        segment_count = len(self.segmentation) - 1
        segments.append((self.segmentation[0], self.segmentation[1]))
        for i in range(1, segment_count):
            segments.append((self.segmentation[i] - self.parameters["overlap"], self.segmentation[i + 1]))
        return segments
    
    #Returns the subsequences corresponding to each segment
    def segment_subsequences(self):
        """Returns the subsequences of each segment corresponding to the stored segmentation
        This function uses the segmentation_to_segments function as a helper function"""
        segments = self.segmentation_to_segments()
        subsequences = []
        for segment in segments:
            subsequences.append(self.parameters["sequence"][segment[0] : segment[1]])
        return subsequences
    
    #Prints the output segmentation along with all the segments
    def print_results(self):
        """Function that will print out a list of all the cuts, the score of the segmentation, and a list of all the segments
        This function also prints some additional information such as the number of segments and the length of each segment not including overlap
        This function will return the printed information as well in order of printing
        This function returns the subsequences of each segment as well but will not print them due to length
        Returns the following in order, segmentation, final_score, segments, segment_count, segment_lengths, subsequences"""
        #Print and return the segmentation
        print(f"The segmentation is {self.segmentation}")
        final_score = self.optimal_scores[len(self.optimal_scores) - 1]
        #Print and return the final score
        print(f"The score of the segmentation is {final_score}")
        segments = self.segmentation_to_segments()
        #Print and return the segments corresponding to the segmentation
        print(f"The segments corresponding to this segmentation are {segments}")
        segment_count = len(segments)
        #Print and return the number of segments
        print(f"There are {segment_count} segments")
        segment_lengths = [0 for i in range(0, segment_count)]
        #Calculate, print, and return the length of each segment
        for i in range(0, segment_count):
            segment_lengths[i] = self.segmentation[i + 1] - self.segmentation[i]
        print(f"The lengths of all the segments not including overlap are {segment_lengths}")
        #Print and return the subsequences of each segment
        subsequences = self.segment_subsequences()
        print(f"The subsequences of all the segments are too long to be printed")
        return self.segmentation, final_score, segments, segment_count, segment_lengths, subsequences

    #Writes everything to a txt file
    def write_subsequences_to_txt(self, filepath):
        """Writes the segment label, starting index, ending index, and subsequence to a given txt file
        Every segment corresponding to the stored segmentation is included in this txt file
        This function uses segment_subsequences and segmentation_to_segments as helper functions"""
        segments = self.segmentation_to_segments()
        subsequences = self.segment_subsequences()
        #Set up header
        lines = [""] * (len(segments) + 1)
        lines[0] = [f"Segment Label\tStart\tEnd\tSubsequence\n"]
        for i in range(0, len(segments)):
            #Each line has the segment name, starting index, ending index, and subsequence of the segment separated by tabs
            lines[i] = f"Segment {i + 1}\t{segments[i][0]}\t{segments[i][1]}\t{subsequences[i]}\n"
        with open(filepath, "w") as f:
            for line in lines:
                f.write(line)
        return filepath
    
    #Writes the segmentation to a bed file that can be imported into SnapGene
    def write_segments_to_bed(self, filepath):
        """Writes the segment label, starting index, and ending index to a given bed file
        Every segment excluding the overlap corresponding to the stored segmentation is included in this bed file
        This function uses the stored segmentation when writing to the bed file"""
        with open(filepath, "w") as f:
            #Set up header
            f.write(f"#Chromosome Number\tStart\tEnd\tLabel\n")
            for i in range(0, len(self.segmentation) - 1):
                #Write segments to the bed file
                f.write(f"0\t{self.segmentation[i]}\t{self.segmentation[i + 1]}\tSegment {i + 1}\n")
        return filepath
    
    #Writes the segmentation and forbidden regions to a bed file that can be imported into SnapGene
    def write_segments_and_forbidden_regions_to_bed(self, filepath):
        """Writes the segment label, starting index, and ending index to a given bed file
        Every segment excluding the overlap corresponding to the stored segmentation is included in this bed file
        This function uses the stored segmentation when writing to the bed file"""
        with open(filepath, "w") as f:
            #Set up header
            f.write(f"#Chromosome Number\tStart\tEnd\tLabel\tScore\tStrand\tThick Start\tThick End\tRGB\n")
            for i in range(0, len(self.segmentation) - 1):
                #Write segments to the bed file
                f.write(f"0\t{self.segmentation[i]}\t{self.segmentation[i + 1]}\tSegment {i + 1}\t0\t+\t0\t0\t166,172,179\n")
            forbidden_regions = self.parameters["forbidden_regions"]
            for i in range(0, len(forbidden_regions) - 1):
                #Write forbidden regions to the bed file
                f.write(f"0\t{forbidden_regions[i][0]}\t{forbidden_regions[i][1]}\tForbidden Region {i + 1}\t0\t+\t0\t0\t255,0,0\n")
        return filepath
    
#Test code
if __name__ == "__main__":
    segment_scoring_functions = [scoring_length.length_score, scoring_forbidden_regions.forbidden_region_score, scoring_overlap_composition.overlap_composition_score, scoring_forbidden_regions.forbidden_region_class_score, scoring_microhomologies.microhomology_score]
    preprocessing_functions = [scoring_forbidden_regions.forbidden_regions, scoring_overlap_composition.GC_proportions, scoring_microhomologies.relevant_microhomologies]
    parameters = {
                    "max_length" : 3000,
                    "min_length" : 1000,
                    "overlap" : 100,
                    "microhomology_distance" : 20,
                    "min_microhomology_length" : 8,
                    "max_microhomology_length" : 19,
                    "forbidden_region_class_count" : 1,
                    "forbidden_regions_from_file" : False,
                    "forbidden_region_generation" : False,
                    "color" : "#ff0000"
                }
    segmenter = segmentron(preprocessing_functions, segment_scoring_functions, scoring_accumulator.addition_function, parameters)
    filepath = "./Hba_Sergio_Test.dna" 
    total_score, segmentation = segmenter.segment_from_file(filepath, multiprocessing_cores = 0, coarseness = 1)
    segmenter.print_results()
    segmenter.write_subsequences_to_txt("segmentation_multiprocessed.txt")
    segmenter.write_segments_to_bed("segmentation_multiprocessed.bed")
    segmenter.write_segments_and_forbidden_regions_to_bed("segmentation_and_forbidden_regions_multiprocessed.bed")
    parameters["forbidden_regions_from_file"] = True
    parameters["forbidden_region_generation"] = True
    segmenter = segmentron(preprocessing_functions, segment_scoring_functions, scoring_accumulator.addition_function, parameters)
    filepath = "./Hba_Sergio_Test.dna" 
    total_score, segmentation = segmenter.segment_from_file(filepath, multiprocessing_cores = 0, coarseness = 1)
    segmenter.print_results()
    segmenter.write_subsequences_to_txt("segmentation.txt")
    segmenter.write_segments_to_bed("segmentation.bed")
    segmenter.write_segments_and_forbidden_regions_to_bed("segmentation_and_forbidden_regions.bed")