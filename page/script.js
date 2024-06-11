var segmented;
var file;
var startTime;
var timerInterval;

const resultFiles = {};

const worker = new Worker('./webworker.js');
var currCt, totalLen;
worker.onmessage = function(event){
    console.log("Message received from worker");
    console.log(event);
    const results = event.data;
    if (results["type"] == "segmentation") {
        clearInterval(timerInterval);
        resultFiles["txt"] = results["segmentation.txt"];
        resultFiles["bed"] = results["segmentation.bed"];
        resultFiles["bedFR"] = results["segmentation_and_forbidden_regions_multiprocessed.bed"];

        pyodide.FS.writeFile("/segmentation.txt", results["segmentation.txt"]);
        pyodide.FS.writeFile("/segmentation.bed", results["segmentation.bed"]);
        pyodide.FS.writeFile("/segmentation_and_forbidden_regions_multiprocessed.bed", results["segmentation_and_forbidden_regions_multiprocessed.bed"]);        
        segmented = true;
        document.getElementById("findSegments").disabled = false;
        document.getElementById("downloadSection").classList.remove("hidden");
        document.getElementById("downloadSection").scrollIntoView({ behavior: 'smooth' });
        document.getElementById("inputFileName").innerHTML = `Input File: ${file.name.split(".")[0]}`;
        alert("Segments are ready!");
    }else if(results["type"] == "output"){
        //progressMsg.scrollTop = progressMsg.scrollHeight;
        if(results["msg"].indexOf("Progress") == 0){
            [,currCt, totalLen] = results["msg"].split("/").map(x => parseInt(x));
            setLoadingProgress(currCt, totalLen);
        }else{
            console.log("Message From Worker: ", results["msg"]);
            const progressMsg = document.getElementById("progressMsg")
            progressMsg.value += `\n${results["msg"]}`;
            if(results["msg"].indexOf("This function")==0){
                setLoadingProgress(totalLen, totalLen);
                document.getElementById('loadingBar').style.backgroundColor = '#007bff';
            }
        }
    }
}
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
    const parameterInputs = document.querySelectorAll("#parameterSection input");
    const parameterValues = {};
    parameterInputs.forEach(parameterInput => {
        parameterValues[parameterInput["name"]] = parameterInput.type === "checkbox" ? parameterInput.checked : Math.ceil(parameterInput.value);
    });
    // Return the map
    return parameterValues;
}

/*
    Starts segmentation from a file
*/

function setLoadingProgress(currCt, totalCt) {
    var loadingBar = document.getElementById('loadingBar');
    loadingBar.style.width = 100*currCt/totalCt+ '%';
    document.getElementById('percentText').innerHTML = `${Math.floor(100*currCt/totalCt)}%`;
    document.getElementById('progressText').innerHTML = `${currCt}/${totalCt} nt`;
}

async function segmentFile(){
    const param = getParameters();
    console.log(param);

//Add Function to print
    const segmentronCode = 
`import time
import json
start_time = time.time()

# Import Segmentron and its modules
import Segmentron
import Segmentron.scoring_forbidden_regions as scoring_forbidden_regions
import Segmentron.scoring_length as scoring_length
import Segmentron.scoring_microhomologies as scoring_microhomologies
import Segmentron.scoring_overlap_composition as scoring_overlap_composition
import Segmentron.scoring_accumulator as scoring_accumulator

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Finished import in {elapsed_time:.2f} seconds")

def verbose(i, total_length):
    if(i%1000==0):
        print(f'Progress/{i}/{total_length}')

start_time = time.time()
segment_scoring_functions = [scoring_length.length_score, scoring_forbidden_regions.forbidden_region_score,${param["GCCount"]?"scoring_overlap_composition.overlap_composition_score,":""} scoring_forbidden_regions.forbidden_region_class_score, scoring_microhomologies.microhomology_score]
preprocessing_functions = [scoring_forbidden_regions.forbidden_regions, ${param["GCCount"]?"scoring_overlap_composition.GC_proportions":""}, scoring_microhomologies.relevant_microhomologies]
parameters = {
                "max_length" : ${param["maxLen"]},
                "min_length" : ${param["minLen"]},
                "overlap" : ${param["overlap"]},
                "microhomology_distance" : ${param["microDist"]},
                "min_microhomology_length" : ${param["minMicroLen"]},
                "max_microhomology_length" : ${param["maxMicroLen"]},
                "forbidden_regions_from_xml": ${param["blast"]?"\"blast_results.xml\"":"None"},
                "forbidden_regions_from_input": True, 
                "forbidden_region_class_count" : 1,
                "forbidden_region_generation" : ${param["blast"]?"True":"False"},
                "color" : "#ff0000",
                "verbose" : verbose
            }

segmenter = Segmentron.segmentron(preprocessing_functions, segment_scoring_functions, scoring_accumulator.addition_function, parameters)
filepath = "/sequence.${file.name.split('.')[1]}"
print("Start Scoring")
total_score, segmentation = segmenter.segment_from_file(filepath, multiprocessing_cores = 0, coarseness = 1)
segmenter.print_results()
segmenter.write_subsequences_to_txt("/segmentation.txt")
segmenter.write_segments_to_bed("/segmentation.bed")
segmenter.write_segments_and_forbidden_regions_to_bed("/segmentation_and_forbidden_regions_multiprocessed.bed")
end_time = time.time()
elapsed_time = end_time - start_time

print(f"Finished Segmenting in {elapsed_time} seconds")
`;



//post message to web worker and activate timer
async function startWebWorkerSegment() {
    try {
        startTime = performance.now();
        // Update the HTML element every second
        timerInterval = setInterval(() => {
            const elapsedTime = (performance.now() - startTime) / 1000; // Convert to seconds
            document.getElementById("elapsedTime").innerText = `Elapsed Time: ${elapsedTime.toFixed(1)}s`;
        }, 1000);

        worker.postMessage({
            python: segmentronCode,
            file: file,
            blast: param["blast"]
        });
    } catch (e) {
        console.log(
        `Error in pyodideWorker at ${e.filename}, Line: ${e.lineno}, ${e.message}`,
        );
    }
}
startWebWorkerSegment();
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
    const content = resultFiles["txt"];//pyodide.FS.readFile("/segmentation.txt");
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
    const content = resultFiles["bed"]; //pyodide.FS.readFile("/segmentation.bed");
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
    const content = resultFiles["bedFR"];//pyodide.FS.readFile("/segmentation_and_forbidden_regions_multiprocessed.bed");
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
function segment(event){
    event.target.disabled  = true;
    // Retrieve the file input element
    var fileInput = document.getElementById("upload");

    // Get the selected file
    
    //var fastaSegment = document.getElementById("segment").value;

    if (fileInput.files[0]) {
        file = fileInput.files[0];
        // File is selected, you can perform further actions
        console.log("Selected file:", file);
        //document.getElementById("importSection").classList.add("hidden");
        document.getElementById("segmentSection").classList.remove("hidden");
        document.getElementById("segmentSection").scrollIntoView({ behavior: 'smooth' });
        segmentFile();
    } else {
        // No file is selected
        /*if (!fastaSegment.trim()) {
            // FASTA segment is empty, send an error message
            alert("Please enter a FASTA segment.");
            return; // Stop further execution
        }*/
        event.target.disabled  = false;
        alert("Please upload a file");

    }
}

const defaultParams = {
    blast: true,
    GCCount: true,
    minLen: 1000,
    maxLen: 3000,
    overlap: 100,
    microDist: 20,
    minMicroLen: 8,
    maxMicroLen: 19
};

document.getElementById("restoreParam").addEventListener("click", ()=>{
    const parameterInputs = document.querySelectorAll("#parameterSection input");
    parameterInputs.forEach(parameterInput => {
        parameterInput[parameterInput.type === "checkbox"?"checked":"value"] = defaultParams[parameterInput["name"]];
    });
});


document.addEventListener("DOMContentLoaded", async function () {
    // Wait for the DOM content to be fully loaded
    // Add an event listener to the "Find Segments" button
    const segmentButton = document.getElementById("findSegments");
    segmentButton.addEventListener("click", segment);

    const downloadButton = document.getElementById("downloadButton");
    downloadButton.addEventListener("click", downloadSegments);

    const restartButton = document.getElementById("restartButton");
    restartButton.addEventListener("click", ()=>{
        document.getElementById("importSection").scrollIntoView({ behavior: 'smooth' });
    });

    // Get references to the scroll buttons and the sections container
    const sectionsContainer = document.getElementById('sections');

    // Add event listeners to scroll buttons
    /*
    scrollLeftBtn.addEventListener('click', () => {
        const scrollAmount = sectionsContainer.offsetWidth * 0.6; // Scroll by 60% of the section width
        sectionsContainer.scrollBy({
            left: -scrollAmount,
            behavior: 'smooth'
        });
    });

    scrollRightBtn.addEventListener('click', () => {
        const scrollAmount = sectionsContainer.offsetWidth * 0.6; // Scroll by 60% of the section width
        sectionsContainer.scrollBy({
            left: scrollAmount,
            behavior: 'smooth'
        });
    });
    */

    var fileInput = document.getElementById("upload");
    segmentButton.classList.remove("disabledButton");
});