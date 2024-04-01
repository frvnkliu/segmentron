# Importing necessary modules from the Biopython package: 
# NcbiblastnCommandline for setting up and running BLAST searches, 
# and NCBIXML for parsing BLAST results.
import subprocess
from Bio.Blast import NCBIXML
import snapgene_reader as sgreader

# Helper function to call BLAST and generate microhomologies
def return_repeats(sequence):
    with open("temp.fasta", "w") as f:
        f.write(">seq\n" + sequence)
    # ADD QUICK EXPLENATION FOR EACH PARAMETER
    # Run the BLAST search. The results are stored in the "blast_results.xml" file.
    blast_command = [
        "blastn",                 # Command for nucleotide-nucleotide BLAST
        "-query", "temp.fasta",   # Input query file in FASTA format
        "-subject", "temp.fasta", # Subject file to BLAST against
        "-out", "blast_results.xml", # Output file in XML format
        "-outfmt", "5",           # Output format (5 for XML)
        "-gapopen", "5",          # Cost to open a gap
        "-gapextend", "5",        # Cost to extend a gap
        "-perc_identity", "80",   # Percentage identity threshold
        "-strand", "both",        # Search both strands
        "-word_size", "20"        # Word size for initial matches
    ]
    # Execute the command using subprocess
    result = subprocess.run(blast_command, capture_output = True, text = True)

    # Check if the command was successful
    if result.returncode == 0:
        print("BLAST search completed successfully.")
    else:
        print("Error in BLAST search:")
        print(result.stderr)
        
    #Open output file
    result_handle = open("blast_results.xml")
    #Read output file
    blast_record = NCBIXML.read(result_handle)
    #List of matches found
    matches = []
    #Iterate over each alignment found
    for alignment in blast_record.alignments:
        for i, hsp in enumerate(alignment.hsps):
            #i>0 removes the sequence aligning to itself.
                if i > 0:
                    #For each alignment made there is a query and a subject, 
                    #the query is the sequence used to find homology, and the 
                    #subject is the found homology. They are both considered 
                    #repeats, so if you print the results you will see that the 
                    #subject/query pairs will be repeated with inverse orders 
                    #in the next alignment. So, we can retreive the coordinate
                    #of either all queries or all subject and get the same list 
                    #of matches. I chose to use the querries:
                    query_start = hsp.query_start
                    query_end = hsp.query_end
                    
                    if (query_start, query_end) not in matches:
                        matches.append((query_start, query_end))
                    
    result_handle.close()
    #Consolidate overlapping matches
    consolidated_matches = []
    matches.sort()
    for start, end in matches:
        if consolidated_matches and consolidated_matches[-1][1] >= start - 1:
            consolidated_matches[-1] = (consolidated_matches[-1][0], max(end, consolidated_matches[-1][1]))
        else:
            consolidated_matches.append((start, end))
    return(consolidated_matches)

#Scoring function that scores a segment relative to the types of forbidden regions contained inside the segment
#Score is minimized when all forbidden regions contained inside the segment are pairwise dissimilar
#Input dictionary must contain a list of lists of forbidden regions that are similar
#Returns a positive value instead of a negative value so smaller values are preferred
def forbidden_region_class_score(parameters, starting_index, ending_index):
    forbidden_region_classes = parameters.get("forbidden_region_classes")
    if not forbidden_region_classes or len(forbidden_region_classes) == 1:
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

#Helper function to read forbidden regions from a file
def forbidden_region_read(parameters):
    filepath = parameters["filepath"]
    dictionary = sgreader.snapgene_file_to_dict(filepath)
    features = dictionary["features"]
    forbidden_region_color = parameters.get("color", "#ff0000")
    color = ""
    forbidden_regions = []
    for feature in features:
        #Assume that all forbidden regions will be the color red (#ff0000)
        color = feature["color"]
        if (color == forbidden_region_color):
            #Add the forbidden region to the list of all forbidden_regions
            forbidden_regions.append((feature["start"], feature["end"]))
    #Sort the forbidden_regions by starting index
    forbidden_regions.sort()
    parameters["forbidden_regions"] = forbidden_regions
    parameters["forbidden_region_classes"] = [forbidden_regions]
    return None

#Helper function to read forbidden region classes from a file
def forbidden_region_classes_read(parameters):
    filepath = parameters["filepath"]
    #Default value for forbidden_region_class_count is 1
    forbidden_region_class_count = parameters.get("forbidden_region_class_count", 1)
    dictionary = sgreader.snapgene_file_to_dict(filepath)
    features = dictionary["features"]
    forbidden_region_color = parameters["color"]
    region = (0, 0)
    region_name = ""
    class_tag = ""
    forbidden_region_classes_list = [[] for i in range(forbidden_region_class_count)]
    for feature in features:
        #Assume that all forbidden regions will be the color red (#ff0000)
        color = feature["color"]
        if (color == forbidden_region_color):
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
    parameters["forbidden_regions"] = forbidden_regions
    for forbidden_region_class in forbidden_region_classes_list:
        forbidden_region_class.sort()
    parameters["forbidden_region_classes"] = forbidden_region_classes_list
    return None

#Preprocessing function to read forbidden regions from a file and generate forbidden regions if necessary
def forbidden_regions(parameters):
    #Default value for forbidden_region_class_count is 1
    forbidden_region_class_count = parameters.get("forbidden_region_class_count", 1)
    forbidden_regions = []
    #Default value for forbidden_region_generation is True
    forbidden_region_generation = parameters.get("forbidden_region_generation", True)
    #Default value for forbidden_regions_from_file is False
    forbidden_regions_from_file = parameters.get("forbidden_regions_from_file", False)
    if not (forbidden_region_generation or forbidden_regions_from_file):
        parameters["forbidden_regions"] = forbidden_regions
        parameters["forbidden_region_classes"] = [forbidden_regions]
        return 0
    #Check if forbidden regions need to be read from the file
    if (parameters["forbidden_regions_from_file"]):
        #Check if forbidden regions need to be grouped into classes
        if forbidden_region_class_count == 1:
            #Call helper function to read forbidden regions from the file
            forbidden_region_read(parameters)
        #Group forbidden regions into classes if necessary
        else:
            forbidden_region_classes_read(parameters)
    #Check if forbidden regions need to be generated
    if (parameters["forbidden_region_generation"]):
        parameters["forbidden_regions"].extend(return_repeats(parameters["sequence"]))
        parameters["forbidden_regions"].sort()
    return 0