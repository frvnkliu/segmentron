// webworker.js

// Setup your project to serve `py-worker.js`. You should also serve
// `pyodide.js`, and all its associated `.asm.js`, `.json`,
// and `.wasm` files as well:

async function loadPyodideAndPackages() {
  const pyodideScript = await import("https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.js");
    // Initialize Pyodide
  //self.pyodide = await pyodideScript.loadPyodide();
  self.pyodide = await loadPyodide();
  //Imports with Micropip
  await pyodide.loadPackage("micropip");
  const micropip = pyodide.pyimport("micropip");
  await micropip.install("biopython");
  await micropip.install("snapgene_reader");
  await micropip.install("typing");
  await micropip.install("Bio");
  await micropip.install("https://test-files.pythonhosted.org/packages/20/cd/fd0670165a01ce60350e114c15c242799348522e987475ec14614e7fd092/Segmentron-9.2.0-py3-none-any.whl");
  await self.pyodide.loadPackage(["numpy", "pytz"]);
}
let pyodideReadyPromise = loadPyodideAndPackages();

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

self.onmessage = async (event) => {
  // make sure loading is done
  await pyodideReadyPromise;
  // Don't bother yet with this line, suppose our API is built in such a way:
  const { id, python, file } = event.data;
  // The worker copies the context in its own "memory" (an object mapping name to values)
  /*for (const key of Object.keys(context)) {
    self[key] = context[key];
  }*/
  const content = await readFileAsync(file);
  const uint8ArrayContent = new Uint8Array(content);
  pyodide.FS.writeFile("/sequence.dna", uint8ArrayContent);
  // Now is the easy part, the one that is similar to working in the main thread:
  try {
    await self.pyodide.loadPackagesFromImports(python);
    await self.pyodide.runPythonAsync(python);
    const output1 = pyodide.FS.readFile("/segmentation.txt");
    const output2 = pyodide.FS.readFile("/segmentation.bed");
    const output3 = pyodide.FS.readFile("/segmentation_and_forbidden_regions_multiprocessed.bed");
    const results = 
    {
      "segmentation.txt": output1,
      "segmentation.bed": output2,
      "segmentation_and_forbidden_regions_multiprocessed.bed": output3
    };
    self.postMessage({ results, id });
  } catch (error) {
    self.postMessage({ error: error.message, id });
  }
};