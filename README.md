### CircularSTAR3D: a stack-based RNA 3D structural alignment tool for rotated matching

* Terms\
CircularSTAR3D is free of charge for non-commercial purposes, and it comes
with ABSOLUTELY NO WARRANTY.

	Where appropriate, please cite the following circularSTAR3D paper:
Chen et al. "CircularSTAR3D: a stack-based RNA 3D structural alignment tool for rotated matching.

	CircularSTAR3D integrates DSSR v1.5.3 (now it's called DSSR basic) in the software package. 
The use of DSSR subjects its own terms:
DSSR is free of charge for non-commercial purposes, and it comes
with ABSOLUTELY NO WARRANTY. 
Where appropriate, please cite the following DSSR paper:
Lu et al. (2015). "DSSR: an integrated software tool for
 dissecting the spatial structure of RNA." Nucleic Acids
 Res., 43(21), e142 (doi:10.1093/nar/gkv716).

* Requirement\
CircularSTAR3D is implemented by using java 1.8 and can be executed in 64-bit 
Linux. Two java packages, "commons-cli-1.2.jar" and "EJML-core-0.26.jar", 
are used in the program to support argument parsing and efficient 
matrix computation. "DSSR" is used for base pairing annotation. 
"Apache ant" is required to compile the source code.

	Apache ant: http://ant.apache.org/bindownload.cgi 
	
	For Linux system, ant can be installed directly by using command line.
Debian/Ubuntu: "apt-get install ant"
Fedora/CentOS: "yum install ant"

	The two jar files and two programs have been already included into 
the package.

* Installation\
Go to the STAR3D home directory and execute "ant jar".
(The users can also clean the compiled files by executing "ant clean")

* Preprocessing\
CircularSTAR3D downloads PDB files and preprocesses them to retrieve the 
secondary structural information.

	Go to the STAR3D home directory and execute 
"java -cp CircularSTAR3D.jar Preprocess [PDB ID] [Chain ID]".

Notice: Make sure the fold "tools" and file "CircularSTAR3D.jar" in the same directory. 

* Structural Alignment\
Go to the STAR3D home directory and execute 
"java -jar CircularSTAR3D.jar [PDB1 ID] [Chain1 ID] [PDB2 ID] [Chain2 ID]"
More options for the program can be seen by executing 
"java -jar CircularSTAR3D.jar -h"

	Notice: Make sure the folds "PDB", "STAR3D_struct_info" and the file "CircularSTAR3D.jar" 
in the same directory. 

* Example
```
java -cp CircularSTAR3D.jar Preprocess 3f2x X
java -cp CircularSTAR3D.jar Preprocess 4frn B
java -jar CircularSTAR3D.jar -s 2 -n 100 3f2x X 4frn B
```

* ACKNOWLEDGEMENTS\
CircularSTAR3D is developed for an NIH funded project (R01GM102515).

	CircularSTAR3D uses java packages Commons CLI and EJML.
Commons CLI: http://commons.apache.org/proper/commons-cli/
EJML: http://ejml.org/wiki/index.php?title=Main_Page

	Third party software DSSR was downloaded from http://forum.x3dna.org/rna-structures/.
  
* CONTACTS\
For bug reports or comments please contact shaojie.zhang@ucf.edu.
