# protTrace
@Author: Arpit Jain
@Email: arpitrb@gmail.com

### protTrace is a simulation-based framework to estimate the evolutionary traceabilities of proteins. ###

General directions to use the tool:

1. In the directory 'toy_example', there exists a program configuration file 'prog.config'. Use this file to set the paths for the various dependencies.
2. Once these dependencies are set, you can simply try running from the toy_example directory the following command:

	python ../bin/protTrace.py -i test.id -c prog.config

3. The program should start running and you can see the outputs being created in:  ../output/YEAST01111

4. One can also provide a fasta file as an input. In that case simply run:

    	python ../bin/protTrace.py -f $yourFile.fasta -c prog.config

   *** NOTE: Use short headers in the fasta file (max. 15 letters) as it will be used to create a folder in the output directory***
  
