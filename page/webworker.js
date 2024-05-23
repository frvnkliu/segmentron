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

      /*wasm_path: `./libs/v86.wasm`,
      memory_size: 512 * 1024 * 1024,
      vga_memory_size: 8 * 1024 * 1024,
      screen_container: document.getElementById("screen_container"),
      initial_state: { url: `./libs/blastImage.bin` },
      filesystem: { baseurl: `${V86_ROOT}/images/debian-9p-rootfs-flat/` },
      autostart: true*/
  data = "";
}
let v86ReadyPromise = startV86();

//pyodide

async function loadPyodideAndPackages() {
  await import("https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.js");
  self.pyodide = await loadPyodide();
  //Imports with Micropip
  await pyodide.loadPackage("micropip");
  const micropip = pyodide.pyimport("micropip");
  await micropip.install("biopython");
  await micropip.install("snapgene_reader");
  await micropip.install("typing");
  await micropip.install("Bio");
  await micropip.install("https://test-files.pythonhosted.org/packages/d7/8b/6ef135c12b2050b1bc2e7642e487f45c6b48544b19246d5d16d0f1ef7ba7/Segmentron-12.0.0-py3-none-any.whl");
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
    /*for (const key of Object.keys(context)) {
      self[key] = context[key];
    }*/
    const content = await readFileAsync(file);
    const uint8ArrayContent = new Uint8Array(content);
    pyodide.FS.writeFile("/sequence.dna", uint8ArrayContent);
    console.log("Python");
    // Now is the easy part, the one that is similar to working in the main thread:
    if(blast){
      //.dna to .fa
      await pyodide.runPythonAsync(
`import snapgene_reader as sgreader
sequence = sgreader.snapgene_file_to_dict("/sequence.dna")["seq"]
with open("/temp.fa", "w") as f:
  f.write(">seq\\n" + sequence)`
      );
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

      //3await emulator.create_file("root/temp.fa", fa_uint8ArrayContent);

      var blastState = true;
      emulator.add_listener("serial0-output-byte", async function(byte){
          var char = String.fromCharCode(byte);
          if(char === "\r")
          {
              return;
          }

          data += char;

          if(data.endsWith(":~#"))
          {
            /*self.postMessage({
              "type" : "output",
              "msg": "OutputTest:\n" + data + "\n"
            });*/
            if(blastState){
              blastState = false;
              emulator.serial0_send("rm temp.fa\n");
              self.postMessage({
                "type" : "output",
                "msg": "Finished Blast Query\n"
              });
              //read file from v86

              var blastOutput;
              while(true){
                try{
                  blastOutput = await emulator.read_file("root/output.xml");
                  break;
                }catch{
                  console.log("File changes not propagated, trying again")
                  await new Promise(resolve => setTimeout(resolve, 10000));
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