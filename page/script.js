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
    await micropip.install("https://test-files.pythonhosted.org/packages/4a/f3/45dfbd7f98836cf9563a1f8f17a74db1c4b8ec5a3257d6ebd6813b73bd09/Segmentron-9.1.0-py3-none-any.whl");
    console.log("Packages Loaded");
/*const segmentronCode = 
`
import time
start_time = time.time()

# Import Segmentron and its modules
import Segmentron
import Segmentron.function_list_scoring as scoring
import Segmentron.function_list_preprocessing as preprocessing

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Finished import in {elapsed_time} seconds")

segment_scoring_functions = [scoring.length_score, scoring.forbidden_region_score, scoring.overlap_composition_score, scoring.forbidden_region_class_score, scoring.microhomology_score]
preprocessing_functions = [preprocessing.GC_proportions]
parameters = {
                "max_length" : 3000,
                "min_length" : 1000,
                "overlap" : 100,
                "microhomology_distance" : 20,
                "min_microhomology_length" : 8,
                "max_microhomology_length" : 19
            }
segmenter = Segmentron.segmentron(preprocessing_functions, segment_scoring_functions, scoring.addition_function, parameters)
filepath = "./sequence.dna"
print("Finished Scoring")
total_score, segmentation = segmenter.segment_from_file(filepath, forbidden_region_classes = 1, multiprocessing_cores = 0, coarseness = 1)
segmenter.print_results()
segmenter.write_subsequences_to_txt("segmentation_multiprocessed.txt")
segmenter.write_segments_to_bed("segmentation_multiprocessed.bed")
segmenter.write_segments_and_forbidden_regions_to_bed("segmentation_and_forbidden_regions_multiprocessed.bed")
`;
pyodide.runPythonAsync(segmentronCode);*/
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
    const parameterSelectors = Array.from(document.getElementById("parameterSelectors").children);
    const parameterValues = {};
    parameterSelectors.forEach(parameterDiv => {
        const parameterInput = parameterDiv.querySelector("input");
        parameterValues[parameterInput["name"]] = parameterInput.type === "checkbox" ? parameterInput.checked : parameterInput.value;
    });

    // Return the map
    return parameterValues;
}

async function segmentFile(){
    const content = await readFileAsync(file);
    const uint8ArrayContent = new Uint8Array(content);
    const param = getParameters()
    console.log(param);
    pyodide.FS.writeFile("/sequence.dna", uint8ArrayContent);
    const segmentronCode = 
`
import time
start_time = time.time()

# Import Segmentron and its modules
import Segmentron
import Segmentron.function_list_scoring as scoring
import Segmentron.function_list_preprocessing as preprocessing

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Finished import in {elapsed_time} seconds")

segment_scoring_functions = [scoring.length_score, scoring.forbidden_region_score, scoring.overlap_composition_score, scoring.forbidden_region_class_score]
preprocessing_functions = [preprocessing.GC_proportions]
parameters = {
                "max_length" : ${param["maxLen"]},
                "min_length" : ${param["minLen"]},
                "overlap" : ${param["overlap"]},
                "microhomology_distance" : ${param["microDist"]},
                "min_microhomology_length" : ${param["minMicroLen"]},
                "max_microhomology_length" : ${param["maxMicroLen"]}
            }
segmenter = Segmentron.segmentron(preprocessing_functions, segment_scoring_functions, scoring.addition_function, parameters)
filepath = "/sequence.dna"
print("Start Scoring")
total_score, segmentation = segmenter.segment_from_file(filepath, forbidden_region_classes = 1, multiprocessing_cores = 0, coarseness = 1)
segmenter.print_results()
segmenter.write_subsequences_to_txt("/segmentation.txt")
segmenter.write_segments_to_bed("/segmentation.bed")
segmenter.write_segments_and_forbidden_regions_to_bed("segmentation_and_forbidden_regions_multiprocessed.bed")
`;
    await pyodide.runPythonAsync(segmentronCode);
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
        /*if (!fastaSegment.trim()) {
            // FASTA segment is empty, send an error message
            alert("Please enter a FASTA segment.");
            return; // Stop further execution
        }*/
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