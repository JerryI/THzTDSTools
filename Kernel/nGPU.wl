BeginPackage["JerryI`TDSTools`nGPU`", {"OpenCLLink`"}]

run;
readyQ;
load;
allocate;
unload;
get;

Begin["`Private`"]

root = DirectoryName[$InputFileName];
clFile = FileNameJoin[{root, "nGPU.cl"}];

readyQ = OpenCLQ;
load = OpenCLMemoryLoad;
allocate = OpenCLMemoryAllocate;
unload = OpenCLMemoryUnload;
get = OpenCLMemoryGet;

If[OpenCLQ[], Module[{},
  PrintTemporal[Style["Compiling to OpenCL", Italic] ];

  run = OpenCLFunctionLoad[File[clFile], "clRun", {
    {"Float", _, "InputOutput"}, 
    {"Float", _, "Input"}, 
    _Integer,
    _Integer,
    _Integer,
    _Integer,
    _Integer,
    "Float",
    "Float",
    "Float",
    "Float"
  }, 256];




  (* compiledQ = True; *)
] ];

End[]
EndPackage[]
