import snapgene_reader as sgreader
import matplotlib.pyplot as plt
import numpy as np
import function_list_preprocessing as preprocessing
import function_list_scoring as scoring
import segmentron_v9 as segmentron

class scorer:
    """The scorer class will create bar plots representing the scores of a given segmentation of a strand of DNA under a given set of scoring functions
    Bar plots can be created for a single scoring function or a stacked bar plot for all stored scoring functions"""
    #Initialize scorer class with a given array of scoring functions, an array of labels corresponding to those functions, and a maximum segment length
    #The two arrays scoring_functions and function_labels must have the same length

    def __init__(self, preprocessing_functions, scoring_functions, function_labels, parameters):
        """Initializes a scorer class with a given list of functions and a corresponding list of labels for those functions
        The scorer class is also initialized with a parameter for the maximum length which may be expanded to accomodate additional functions
        Stored sequence, segmentation, and forbidden regions are initialized as empty"""
        #List of preprocessing functions that will be used
        self.preprocessing_functions = preprocessing_functions
        #List of scoring functions that will be used 
        self.scoring_functions = scoring_functions
        #List of scoring function labels for the graph's legend
        self.function_labels = function_labels
        self.parameters = dict(parameters)
        #These parameters will be set whenever the class reads from a file
        #Sequence that is segmented
        self.parameters["sequence"] = ""
        #Array of ordered pairs containing the start and end indices of all forbidden_regions
        self.parameters["forbidden_regions"] = []
        #Starting indices of all segments with the last element being the maximum index of the sequence
        self.segmentation = []

    #Goes through all preprocessing functions given and collects all infofrmation necessary
    def preprocess(self):
        for preprocessing_function in self.preprocessing_functions:
            preprocessing_function(self.parameters)
        return None

    #Read from a given segmentron object and score the stored segmentation within segmenter
    #Take the parameters from segmenter and score the stored segmentation
    #Links the parameters of file_scorer to the givern segmenter allowing for faster scoring after every segmentation
    def score_from_segmentron(self, segmenter):
        self.segmentation = segmenter.segmentation
        self.parameters = segmenter.parameters
        self.graph_score_breakdown()

    #Calculates the scores for each segment given a certain scoring function
    #Returns an array of scores under the given scoring function for each segment
    #Assumes that all scoring_functions take input as the starting index, ending index, sequence, and the forbidden regions
    def find_score(self, scoring_function):
        """Subfunction used to calculate an array of scores under a given scoring function
        Each entry in the array corresponds to a segment in the stored segmentation"""
        scores = np.zeros(len(self.segmentation) - 1)
        for i in range(len(scores)):
            starting_index = self.segmentation[i]
            ending_index = self.segmentation[i + 1]
            scores[i] = scoring_function(self.parameters, starting_index, ending_index)
        return scores
    
    #Calculates the total score for each segments based on all scoring functions given
    #Returns an array of the total scores for each segment
    #Assumes that all scoring_functions take input as the starting index, ending index, sequence, and the forbidden regions
    def find_total_score(self):
        """Function used to calculate an array of scores under all stored scoring functions
        Each entry in the array corresponds to a segment in the stored segmentation"""
        scores = np.zeros(len(self.segmentation) - 1)
        for scoring_function in self.scoring_functions:
            scores = scores + self.find_score(scoring_function)
        return scores
    
    #Creates a bar graph of the total score for each segment
    def graph_total_score(self):
        """Function used to create a bar graph of all the total scores under all the stored scoring functions
        Each bar corresponds to a segment in the stored segmentation"""
        scores = self.find_total_score()
        plt.bar(np.arange(1, len(scores) + 1, 1), scores, 0.5)
        plt.title("Total Score")
        plt.xlabel("Segment Number")
        plt.ylabel("Total Score")
        plt.show()

    #Creates a stacked bar graph of the breakdown of the total score between all scoring functions
    def graph_score_breakdown(self):
        """Function used to create a stacked bar graph of all the total scores under all the stored scoring functions
        Each bar corresponds to a segment in the stored segmentation
        Each layer corresponds to the subscore under a specific scoring function"""
        bottom = np.zeros(len(self.segmentation) - 1)
        segment_numbering = np.arange(1, len(bottom) + 1, 1)
        #Create a dictionary with labels for all input data
        scores = {}
        for i in range(len(self.scoring_functions)):
            scores.update({self.function_labels[i] : self.find_score(self.scoring_functions[i])}) 
        for score_type, score_list in scores.items():
            plot = plt.bar(segment_numbering, score_list, 0.5, bottom, label = score_type)
            bottom = np.add(bottom, score_list)
        #Bar graph formatting
        plt.xticks(segment_numbering)
        plt.title("Total Score")
        plt.xlabel("Segment Number")
        plt.ylabel("Score")
        plt.legend(loc = "upper left")
        plt.show()