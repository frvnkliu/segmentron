//import { segmentronCode } from './segemntron.js';

let pyodideLoaded = false;
let pyodide;
var segmented;
var file;
//var content;

async function loadPackages() {
    pyodide = await loadPyodide();

    //Imports with Micropip
    await pyodide.loadPackage("micropip");
    const micropip = pyodide.pyimport("micropip");
    await micropip.install("biopython");
    await micropip.install("snapgene_reader");
    await micropip.install("typing");
    await micropip.install("Bio");
    await micropip.install("https://test-files.pythonhosted.org/packages/9b/c2/20c374253e8754b48ab2ba7d14f97e408e0aa3ff4aa3a8580f51219762d0/Segmentron-8.0.0-py3-none-any.whl");
}

function readFileAsync(file) {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();

        reader.onload = (event) => {
            resolve(event.target.result);
        };

        reader.onerror = (error) => {
            reject(error);
        };
        reader.readAsArrayBuffer(file);
    });
}

function getParameters() {
    // Get the values from the input elements within the parameter section
    var blastValue = document.getElementById("blastCheckbox").checked;
    var gcCountValue = document.getElementById("gcCountCheckbox").checked;
    var function1Value = document.getElementById("function1Checkbox").checked;
    var function2Value = document.getElementById("function2Checkbox").checked;
    var maxLenValue = document.getElementById("maxLen").value;
    var microLenValue = document.getElementById("MicroLen").value;

    // Create a map with ID as key and corresponding values
    var parameterValues = new Map([
        ["blastCheckbox", blastValue],
        ["gcCountCheckbox", gcCountValue],
        ["function1Checkbox", function1Value],
        ["function2Checkbox", function2Value],
        ["maxLen", maxLenValue],
        ["MicroLen", microLenValue]
    ]);

    // Return the map
    return parameterValues;
}

async function segmentFile(){
    const content = await readFileAsync(file);
    const uint8ArrayContent = new Uint8Array(content);
    pyodide.FS.writeFile("/sequence.dna", uint8ArrayContent);
    pyodide.runPythonAsync(segmentronCode2);
}


function downloadSegmentsTxt(){
    if(!segmented){
        alert("Please wait for segmentation to finish");
        return;
    }
    const fileName = file.name.split(".")[0];
    const content = pyodide.FS.readFile("/segmentation.txt");
    const blob = new Blob([content],  {type: "text/plain"});
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${fileName}_segmentation.txt`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);

}
function downloadSegmentsBed() {
    if(!segmented){
        alert("Please wait for segmentation to finish");
        return;
    }
    // Split into multiple download buttons
    const fileName = file.name.split(".")[0];
    const content = pyodide.FS.readFile("/segmentation.bed");
    const blob = new Blob([content],  {type: "application/octet-stream"});
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${fileName}_segmentation.bed`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
}

function downloadSegments(){
    const exportType = document.getElementById('exportFormat').value;
    switch(exportType){
        case "bed":
            downloadSegmentsBed();
            break;
        case "txt":
            downloadSegmentsTxt();
            break;
    }
}

async function segment(){
    if (!pyodideLoaded) {
        // Pyodide is not yet loaded, display a message or handle it appropriately
        alert("Please wait for Pyodide and packages to be loaded.");
        return;
    }
    // Retrieve the file input element
    var fileInput = document.getElementById("upload");

    // Get the selected file
    file = fileInput.files[0];
    //var fastaSegment = document.getElementById("segment").value;

    if (file) {
        // File is selected, you can perform further actions
        console.log("Selected file:", file);
        document.getElementById("importSection").classList.add("hidden");
        document.getElementById("status").innerHTML = "Calculating Optimal Segmentation... (Please Wait, May Take a While)";
        document.getElementById("segmentSection").classList.remove("hidden");
        await segmentFile();
        segmented = true;
        document.getElementById("segmentSection").classList.add("hidden");
        document.getElementById("downloadSection").classList.remove("hidden");
        document.getElementById("status").innerHTML = "Segmentation Finished!";
        alert("Segments are ready!");
    } else {
        // No file is selected
        if (!fastaSegment.trim()) {
            // FASTA segment is empty, send an error message
            alert("Please enter a FASTA segment.");
            return; // Stop further execution
        }
    }



document.addEventListener("DOMContentLoaded", async function () {
    // Wait for the DOM content to be fully loaded
    // Add an event listener to the "Find Segments" button
    const segmentButton = document.getElementById("findSegments");
    segmentButton.addEventListener("click", segment);

    const downloadButton = document.getElementById("downloadButton");
    downloadButton.addEventListener("click", downloadSegments);
    await loadPackages();
    pyodideLoaded = true;
    var fileInput = document.getElementById("upload");
    document.getElementById("status").innerHTML = "Ready!";
    segmentButton.classList.remove("disabledButton");
 });

const segmentronCode2 = 
`
from Segmentron import function_list
from Segmentron import segmentron_v8 as segmentron
`

const segmentronCode1 = 
`
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def run_blast(query_sequence, database='nr', num_results=10):
    result_handle = NCBIWWW.qblast("blastp", database, query_sequence, alignments=num_results)
    blast_record = NCBIXML.read(result_handle)
    return blast_record

def print_blast_results(blast_record):
    for alignment in blast_record.alignments:
        print(f"*** Alignment: {alignment.title} ***")
        for hsp in alignment.hsps:
            print(f"E-value: {hsp.expect}")
            print(f"Identity: {hsp.identities}/{hsp.align_length} ({hsp.identities / hsp.align_length * 100:.2f}%)")
            print(f"Query: {hsp.query}")
            print(f"Subject: {hsp.sbjct}")
            print()

if __name__ == "__main__":
    # Replace the sequence below with your own protein sequence
    query_sequence = "MKTIIALSYIFCLVTGGLTLIGNILISLVVIVRFAKMKKLR"
    
    # Run BLAST
    blast_result = run_blast(query_sequence)
    
    # Print results
    print_blast_results(blast_result)
`

const segmentronCode =
`
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
    #print('\\n', bigger_homologies)

    bigger_homologies = eliminate_microhomologies_in_big_repeats(bigger_homologies, [(24011, 24458)])
    #print('\\n', bigger_homologies)

    bigger_homologies = remove_repeats_that_form_a_single_sequence(bigger_homologies)
    #print('\\n', bigger_homologies)

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
                    bed_file.write(f"chromosome\\t{start}\\t{end+1}\\t{name_field}\\n")
                unique_identifier += 1

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
    microhomologies = preprocess_microhomologies(search_sequence)
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

import snapgene_reader as sgreader
import math
from typing import List, Callable
from tqdm import tqdm
import time

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
    #Stored values for usage in dynamic programming
    optimal_cut: int

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
        self.parameters = parameters
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
        #Set up storage areas for later use in dynamic programming
        #Set after reading from a file
        self.parameters["sequence"] = ""
        self.parameters["forbidden_regions"] = []
        #Set to store the segmentation for future operations following every segmentation
        self.segmentation = []

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
    def segment_from_file(self, filepath, forbidden_region_classes = 1, coarseness = 1, multiprocessing_cores = 0):
        """Given a file path to a valid snapgene file, the segmentron will read the sequence and forbidden regions from the file
        Since the given file may have multiple features, only features marked with the color red (#ff0000) will be counted as forbidden regions
        This function returns the optimal score and the corresponding optimal segmentation by calling a helper function
        This helper function will also store the sequence and forbidden regions that were read from the file"""
        print(filepath)
        dictionary = sgreader.snapgene_file_to_dict(filepath)
        sequence = dictionary["seq"]
        features = dictionary["features"]
        color = ""
        forbidden_regions = []
        #Check if forbidden regions need to be grouped into classes
        if forbidden_region_classes == 1:
            for feature in features:
                #Assume that all forbidden regions will be the color red (#ff0000)
                color = feature["color"]
                if (color == "#ff0000"):
                    #Add the forbidden region to the list of all forbidden_regions
                    forbidden_regions.append((feature["start"], feature["end"]))
            forbidden_regions.sort()
            self.parameters["forbidden_region_classes"] = [forbidden_regions]
        else:
            region = (0, 0)
            region_name = ""
            class_tag = ""
            forbidden_region_classes_list = [[] for i in range(forbidden_region_classes)]
            for feature in features:
                #Assume that all forbidden regions will be the color red (#ff0000)
                color = feature["color"]
                if (color == "#ff0000"):
                    region = (feature["start"], feature["end"])
                    region_name = feature["name"]
                    #Add the forbidden region to the list of all forbidden_regions
                    forbidden_regions.append(region)
                    #Check the region_name to see which class of forbidden region it belongs in
                    #Names of features in the input file should be prefixed as "class_tag-number"
                    #class_tag should be an integer in the range [1, forbidden_region_classes]
                    class_tag = region_name[0 : region_name.index("-")]
                    forbidden_region_classes_list[int(class_tag) - 1].append(region)
            #Sort the forbidden_regions by starting index
            forbidden_regions.sort()
            for forbidden_region_class in forbidden_region_classes_list:
                forbidden_region_class.sort()
            parameters["forbidden_region_classes"] = forbidden_region_classes_list
        return self.segment(sequence, forbidden_regions, coarseness, multiprocessing_cores)

    #Segment a sequence after being given the sequence and a set of forbidden regions
    def segment(self, sequence, forbidden_regions, coarseness = 1, multiprocessing_cores = 0):
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
        #Preprocess the sequence and store relevant information
        self.preprocess()
        #Set up dynamic programming arrays for optimal scores
        if (multiprocessing_cores == 0):
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
        #If only a single process is used, call subfunctions to perform the dynamic programming
        #If coarseness > 1, certain indices will be skipped which will reduce runtime by a factor of about coarseness^2
        for current_index in range(0, sequence_length + 1, coarseness):
            #For each index, check previous cases
            self.segment_subsequence(current_index)
            #Update progress bar
            progress_bar.update(1)
        #Ensure that the entire sequence is sequenced in case it was skipped due to the coarseness setting
        self.segment_subsequence(sequence_length)
        progress_bar.update(1)
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
    
    #Goes through all preprocessing functions given and collects all infofrmation necessary
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
        lines[0] = [f"Segment Label\\tStart\\tEnd\\tSubsequence\\n"]
        for i in range(0, len(segments)):
            #Each line has the segment name, starting index, ending index, and subsequence of the segment separated by tabs
            lines[i] = f"Segment {i + 1}\\t{segments[i][0]}\\t{segments[i][1]}\\t{subsequences[i]}\\n"
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
            f.write(f"#Chromosome Number\\tStart\\tEnd\\tLabel\\n")
            for i in range(0, len(self.segmentation) - 1):
                #Write segments to the bed file
                f.write(f"0\\t{self.segmentation[i]}\\t{self.segmentation[i + 1]}\\tSegment {i + 1}\\n")
        return filepath

#Test code
if __name__ == "__main__":
    segment_scoring_functions = [length_score, forbidden_region_score, overlap_composition_score, forbidden_region_class_score, microhomology_score]
    preprocessing_functions = [GC_proportions, relevant_microhomologies]
    parameters = {
                    "max_length" : 3000,
                    "min_length" : 1000,
                    "overlap" : 100,
                    "microhomology_distance" : 20
                }
    segmenter = segmentron(preprocessing_functions, segment_scoring_functions, addition_function, parameters)
    filepath = "/sequence.dna"
    total_score, segmentation = segmenter.segment_from_file(filepath, forbidden_region_classes = 1, multiprocessing_cores = 0, coarseness = 1)
    segmenter.print_results()
    segmenter.write_subsequences_to_txt("/segmentation.txt")
    print("Written to segmentation.txt")
    segmenter.write_segments_to_bed("/segmentation.bed")
    print("Written to segmentation.bed")
 `;