## CircularSTAR3D: a stack-based RNA 3D structural alignment tool for circular matching

### Terms

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)  

Where appropriate, please cite the following circularSTAR3D paper:
Chen et al. "CircularSTAR3D: a stack-based RNA 3D structural alignment tool for circular matching.

### Requirement
CircularSTAR3D is implemented by using java 1.8 and can be executed in 64-bit 
Linux. Two java packages, "commons-cli-1.2.jar" and "EJML-core-0.26.jar", 
are used in the program to support argument parsing and efficient 
matrix computation. The two jar files have been already included into the package.

 "Apache ant" is required to compile the source code.

Apache ant: http://ant.apache.org/bindownload.cgi 
For Linux system, ant can be installed directly by using command line.
Debian/Ubuntu: "apt-get install ant"
Fedora/CentOS: "yum install ant"

CircularSTAR3D extracts base-pair information from DSSR's annotations. The DSSR annotation files that were used to generate the results in our manuscript are provided along with CircularSTAR3D package. There are two ways for users to use DSSR for their PDB files.
1. If you have DSSR software, copy DSSR into CircularSTAR3D/tools/. After that, you will have CircularSTAR3D/tools/DSSR/x3dna-dssr.
2. If you don't have DSSR software, you can copy DSSR annotation files for your PDBs into CircularSTAR3D/STAR3D_struct_info/.
   For example, if you want to preprocess PDB 6d90, you can copy your 6d90.dssr file into CircularSTAR3D/STAR3D_struct_info/.
   

### Installation
In home directory
"chmod +x ./setup.sh && ./setup.sh"

### Preprocessing
CircularSTAR3D downloads PDB files and preprocesses them to retrieve the 
secondary structural information.

Go to the home directory and execute 
"python3 run_circularSTAR3D.py --preprocess --pdb-id [PDB ID] --chain-id [CHAIN ID]".

If you want to use your custom PDB files, copy them into CircularSTAR3D/PDB/ or LocalSTAR3D/PDB/ before running the preprocess command. CircularSTAR3D will skip downloading from PDB website if the PDB files present in the PDB folder. Please use suffix ".pdb" for PDB format and suffix ".cif" for PDBx/mmcif format files.    

Notice: Make sure the fold "tools" and file "CircularSTAR3D.jar" in the same directory. If you see "Exception in thread "main" java.lang.UnsupportedClassVersionError", try to delete the files in class/ folder and run "ant jar".

### Structural Alignment
Go to the CircularSTAR3D home directory and execute 
"python3 run_circularSTAR3D.py [options]"
More options for the program can be seen by executing 
"python3 run_circularSTAR3D.py -h"

Use "--non-rotate" to generate local alignment from LocalSTAR3D.  
Use "--print-PDB n" to generate PDB for top n alignments. The output PDB files will be named by using the same prefix of the output alingment file. For example, if the prefix is "4frn_B_3f2x_X", the output alignment file will be 4frn_B_3f2x_X.aln and the PDB file of the first alignment will be 4frn_B_3f2x_X.aln1.pdb. The default prefix is "output". The users can add path string to the prefix. The output files will be in the CircularSTAR3D root directory by default.  

Notice: Make sure the folds "PDB", "STAR3D_struct_info" and the file "CircularSTAR3D.jar" 
in the same directory. 

### Example
```
python3 run_circularSTAR3D.py --preprocess --pdb-id 3f2x --chain-id X
python3 run_circularSTAR3D.py --preprocess --pdb-id 4frn --chain-id B
python3 run_circularSTAR3D.py --print-PDB 1 --minimum-stack-size 2 --pdb-id1 4frn --chain-id1 B --pdb-id2 3f2x --chain-id2 X
python3 run_circularSTAR3D.py --non-rotate --minimum-stack-size 2 --pdb-id1 4frn --chain-id1 B --pdb-id2 3f2x --chain-id2 X
```

### ACKNOWLEDGEMENTS
CircularSTAR3D is developed for an NIH funded project (R01GM102515).

CircularSTAR3D uses java packages Commons CLI and EJML.
Commons CLI: http://commons.apache.org/proper/commons-cli/
EJML: http://ejml.org/wiki/index.php?title=Main_Page
  
### CONTACT
Xiaoli Chen is the author of CircularSTAR3D. For bug reports or comments please contact xiaoli.chen@knights.ucf.edu.
