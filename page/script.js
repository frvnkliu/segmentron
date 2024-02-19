import { asyncRun } from "./py-worker.js";

let pyodideLoaded = false;

let pyodide;
var segmented;
var file;
async function loadPackages() {
    pyodide = await loadPyodide();
    console.log("Packages Loaded");
/*const script = 
    `
import statistics
import json
def obj(arg1, arg2):
    return {
        "arg1": arg1,
        "arg2": arg2
    }
print("Hrllo2")
obj(123, 245)
print(obj(123, 245))
json.dumps(obj(123, 245))
`;
const resultString = await pyodide.runPythonAsync(script);
const resObj = JSON.parse(resultString);
console.log(resObj);*/
};


/*
    Async file read given input: File file
*/
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


/*
    Retrieves parameters from the parameters section as an object
*/
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

/*
    Starts segmentation from a file
*/
async function segmentFile(){
    const param = getParameters();
    console.log(param);
    const segmentronCode = 
`import time
import json
start_time = time.time()

# Import Segmentron and its modules
import Segmentron
import Segmentron.function_list_scoring as scoring
import Segmentron.function_list_preprocessing as preprocessing

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Finished import in {elapsed_time} seconds")

start_time = time.time()
segment_scoring_functions = [scoring.length_score, scoring.forbidden_region_score, ${param["GCCount"]?"scoring.overlap_composition_score, ":""} scoring.forbidden_region_class_score, scoring.microhomology_score]
preprocessing_functions = [${param["GCCount"]?"preprocessing.GC_proportions, ":""}preprocessing.relevant_microhomologies]
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
segmenter.write_segments_and_forbidden_regions_to_bed("/segmentation_and_forbidden_regions_multiprocessed.bed")
end_time = time.time()
elapsed_time = end_time - start_time

print(f"Finished Segmenting in {elapsed_time} seconds")
json.dumps(segmenter.encodingJson())`;
async function main() {
    try {
        const { results, error } = await asyncRun(segmentronCode, file);
        if (results) {
            pyodide.FS.writeFile("/segmentation.txt", results["segmentation.txt"]);
            pyodide.FS.writeFile("/segmentation.bed", results["segmentation.bed"]);
            pyodide.FS.writeFile("/segmentation_and_forbidden_regions_multiprocessed.bed", results["segmentation_and_forbidden_regions_multiprocessed.bed"]);        
            segmented = true;
            document.getElementById("segmentSection").classList.add("hidden");
            document.getElementById("downloadSection").classList.remove("hidden");
            document.getElementById("status").innerHTML = "Segmentation Finished!";
            document.getElementById("inputFileName").innerHTML = `Input File: ${file.name.split(".")[0]}`;
            alert("Segments are ready!");
        } else if (error) {
        console.log("pyodideWorker error: ", error);
        }
    } catch (e) {
        console.log(
        `Error in pyodideWorker at ${e.filename}, Line: ${e.lineno}, ${e.message}`,
        );
    }
}
main();
}

/*
    Downloads .txt output file of segmentation
*/
function downloadSegmentsTxt(){
    if(!segmented){
        alert("Please wait for segmentation to finish");
        return;
    }
    const fileName = file ? file.name.split(".")[0]: "";
    const content = pyodide.FS.readFile("/segmentation.txt");
    const blob = new Blob([content],  {type: "text/plain"});
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${fileName}_segmentation.txt`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
}

/*
    Downloads .bed output file of segmentation
*/
function downloadSegmentsBed() {
    if(!segmented){
        alert("Please wait for segmentation to finish");
        return;
    }
    // Split into multiple download buttons
    const fileName = file ? file.name.split(".")[0]: "";
    const content = pyodide.FS.readFile("/segmentation.bed");
    const blob = new Blob([content],  {type: "application/octet-stream"});
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${fileName}_segmentation.bed`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
}

/*
    Downloads .bed output file of segmentation + forbidden regions
*/
function downloadSegmentsBedFR() {
    if(!segmented){
        alert("Please wait for segmentation to finish");
        return;
    }
    // Split into multiple download buttons
    const fileName = file ? file.name.split(".")[0]: "";
    const content = pyodide.FS.readFile("/segmentation_and_forbidden_regions_multiprocessed.bed");
    const blob = new Blob([content],  {type: "application/octet-stream"});
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${fileName}_segmentation_and_forbidden_regions_multiprocessed.bed`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
}

/*
    Downloads segments of code
*/
function downloadSegments(){
    const exportType = document.getElementById('exportFormat').value;
    switch(exportType){
        case "bed":
            downloadSegmentsBed();
            break;
        case "bedFR":
            downloadSegmentsBedFR();
            break;
        case "txt":
            downloadSegmentsTxt();
            break;
    }
}

/*
    Starts segmentation
*/
async function segment(){
    if (!pyodideLoaded) {
        // Pyodide is not yet loaded, display a message or handle it appropriately
        alert("Please wait for Pyodide and packages to be loaded.");
        return;
    }
    // Retrieve the file input element
    var fileInput = document.getElementById("upload");

    // Get the selected file
    
    //var fastaSegment = document.getElementById("segment").value;

    if (fileInput.files[0]) {
        file = fileInput.files[0];
        // File is selected, you can perform further actions
        console.log("Selected file:", file);
        document.getElementById("importSection").classList.add("hidden");
        document.getElementById("status").innerHTML = "Calculating Optimal Segmentation... (Please Wait, May Take a While)";
        document.getElementById("segmentSection").classList.remove("hidden");
        segmentFile();
    } else {
        // No file is selected
        /*if (!fastaSegment.trim()) {
            // FASTA segment is empty, send an error message
            alert("Please enter a FASTA segment.");
            return; // Stop further execution
        }*/
        alert("Please upload a file");

    }
}


document.addEventListener("DOMContentLoaded", async function () {
    // Wait for the DOM content to be fully loaded
    // Add an event listener to the "Find Segments" button
    const segmentButton = document.getElementById("findSegments");
    segmentButton.addEventListener("click", segment);

    const downloadButton = document.getElementById("downloadButton");
    downloadButton.addEventListener("click", downloadSegments);

    const restartButton = document.getElementById("restartButton");
    restartButton.addEventListener("click", ()=>{
        document.getElementById("importSection").classList.remove("hidden");
        document.getElementById("downloadSection").classList.add("hidden");
    });
    await loadPackages();
    pyodideLoaded = true;
    var fileInput = document.getElementById("upload");
    document.getElementById("status").innerHTML = "Ready!";
    segmentButton.classList.remove("disabledButton");
});
