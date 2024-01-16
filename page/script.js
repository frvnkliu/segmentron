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
    await micropip.install("https://github.com/frvnkliu/segmentron/raw/9b6b8afe49e5c2be4ebb3d026770ba5a93535a8e/segmentron/dist/Segmentron-9.0.1-py3-none-any.whl");
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


const segmentronCode=
`
from Segmentron import function_list
from Segmentron import segmentron_v8 as segmentron
`;