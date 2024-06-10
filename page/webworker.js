// webworker.js

// Setup your project to serve `py-worker.js`. You should also serve
// `pyodide.js`, and all its associated `.asm.js`, `.json`,
// and `.wasm` files as well:

var emulator;
var data = "";

async function startV86(){
  await import("./libs/libv86.js");
  //v86
  emulator = new V86({
    wasm_path: `./libs/v86.wasm`,
    memory_size: 512 * 1024 * 1024,
    vga_memory_size: 8 * 1024 * 1024,
    initial_state: { url: `./libs/blastImage.zst` },
    filesystem: { baseurl: `./libs/debian-9p-rootfs-flat/` },
    autostart: true,
    disable_keyboard: true,
    disable_mouse: true
  });
  data = "";
}
let v86ReadyPromise = startV86();

//pyodide

async function loadPyodideAndPackages() {
  //await import("https://cdn.jsdelivr.net/pyodide/v0.21.0/full/");
  importScripts("https://cdn.jsdelivr.net/pyodide/v0.26.1/full/pyodide.js");
  self.pyodide = await loadPyodide();
  //Imports with Micropip
  await pyodide.loadPackage("micropip");
  const micropip = pyodide.pyimport("micropip");
  await micropip.install("biopython");
  await micropip.install("snapgene_reader");
  await micropip.install("typing");
  await micropip.install("Bio");
  await micropip.install("https://test-files.pythonhosted.org/packages/49/58/22e1594692e4824c80a70a763a3fc0b03ce9f784a0221e364e0491809adb/Segmentron-13.4.2-py3-none-any.whl");
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



async function runPython(python){
    await self.pyodide.loadPackagesFromImports(python);
    await self.pyodide.runPythonAsync(python);
    const output1 = pyodide.FS.readFile("/segmentation.txt");
    const output2 = pyodide.FS.readFile("/segmentation.bed");
    const output3 = pyodide.FS.readFile("/segmentation_and_forbidden_regions_multiprocessed.bed");
    const results = 
    {
      "type": "segmentation",
      "segmentation.txt": output1,
      "segmentation.bed": output2,
      "segmentation_and_forbidden_regions_multiprocessed.bed": output3
    };
    self.postMessage(results);
}

self.onmessage = async (event) => {
  try {
    // make sure loading is done
    await pyodideReadyPromise;

    pyodide.setStdout({ batched: (msg) => self.postMessage({type: "output", msg})  });
    // Don't bother yet with this line, suppose our API is built in such a way:
    const {python, file, blast} = event.data;
    // The worker copies the context in its own "memory" (an object mapping name to values)
    const content = await readFileAsync(file);
    const parts = file.name.split('.');
    if (parts.length === 1 || (parts[0] === "" && parts.length === 2)) {
      throw new Error("Missing File Extension");
    }

    const extension = parts[1];
    const uint8ArrayContent = new Uint8Array(content);
    pyodide.FS.writeFile(`/sequence.${extension}`, uint8ArrayContent);
    // Now is the easy part, the one that is similar to working in the main thread:
    if(blast){
      let blastStartTime = Date.now(); // Get the current time in milliseconds
      //.dna to .fa
      await pyodide.runPythonAsync(
`
from Bio import SeqIO
import snapgene_reader as sgreader

filepath = "/sequence.${extension}"
extension = "${extension}"

sequence = ""
if(extension == "dna"):
  dictionary = sgreader.snapgene_file_to_dict(filepath)
  sequence = dictionary["seq"]
elif(extension == "txt"):
  with open(filepath, "r") as file:
    sequence = file.read()
elif(extension == "fa" or extension == "fasta"):
  sequence = str(SeqIO.read(filepath, "fasta").seq)
else:
  raise ValueError("Unsupported File Type")

with open("/temp.fa", "w") as f:
  f.write(">seq\\n" + sequence)`
);
      self.postMessage({
        "type" : "output",
        "msg": "Starting Blast Query\n"
      });
      //create blast_results.xml
      const fa_file =pyodide.FS.readFile("/temp.fa");
      const fa_blob = new Blob([fa_file],  {type: "application/octet-stream"});
      const fa_content = await readFileAsync(fa_blob);
      const fa_uint8ArrayContent = new Uint8Array(fa_content);

      await(v86ReadyPromise);
      await emulator.run();
      while(true){
        try{
          await emulator.create_file("root/temp.fa", fa_uint8ArrayContent);
          break;
        }catch{
          console.log("9front File system not ready yet")
          await new Promise(resolve => setTimeout(resolve, 1000));
        }
      }

      var blastState = 2;
      emulator.add_listener("serial0-output-byte", async function(byte){
          var char = String.fromCharCode(byte);
          if(char === "\r")
          {
              return;
          }

          data += char;

          if(data.endsWith(":~#"))
          {
            if(blastState == 2){
              //sync
              blastState = 1;
              emulator.serial0_send("sync\n");
            }if(blastState ==1){
              blastState = 0;
              emulator.serial0_send("rm temp.fa\n");
              let blastEndTime = Date.now(); 
              self.postMessage({
                "type" : "output",
                "msg": `Finished Blast Query (${((blastEndTime-blastStartTime)/1000).toFixed(2)} seconds)\n`
              });
              //read file from v86

              var blastOutput;
              while(true){
                try{
                  blastOutput = await emulator.read_file("root/output.xml");
                  break;
                }catch{
                  console.log("File changes not propagated, trying again")
                  await new Promise(resolve => setTimeout(resolve, 5000));
                }
              }

              const blast_blob = new Blob([blastOutput], { type: 'application/octet-stream' });
              const blastContent = await readFileAsync(blast_blob);
              const blast_uint8ArrayContent = new Uint8Array(blastContent);
              //write contents into pyodide
              pyodide.FS.writeFile("blast_results.xml", blast_uint8ArrayContent);
              runPython(python);
            }
          }
      });
      emulator.serial0_send("blastn -query temp.fa -subject temp.fa -out output.xml -outfmt 5 -gapopen 5 -gapextend 5 -perc_identity 80 -strand both -word_size 20\n");
    }else{
      runPython(python);
    }
  } catch (error) {
    self.postMessage({error: error.message});
  }
};
