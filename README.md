# protTrace
@Author: Arpit Jain
@Email: arpitrb@gmail.com

### protTrace is a simulation-based framework to estimate the evolutionary traceabilities of proteins. ###

################################################
Software dependencies:
################################################

1. Multiple Sequence Alignment 		- MAFFT
2. Tree reconstruction			- RAxML
3. HMMER tools				- hmmscan, hmmfetch
4. MSA format change			- ClustalW
5. For blast				- blastall
6. For protein BLAST			- blastp
7. Formatting BLAST db			- formatdb
8. Creating BLAST db			- makeblastdb
9. Plotting decay curves		- R
10. Orthologs search 			- HaMStR
11. Orthologs search (single seq)	- HaMStR-OneSeq		 
12. Likelihood distances calc.		- TreePuzzle
13. Programming languages		- Python (2.6 or higher), Perl, JAVA

################################################
Other programs / files (provided in folder 'used_files')
################################################

1. Simulating protein sequence		- REvolver
2. Script to degap MSA			- degapper.pl
3. TreePuzzle parameters file		- paramsMaxLikelihoodMapping.txt
4. Species tree for representing results- euNucWithoutRH.nw
5. Species tree's MSA			- euNucWithoutRH.nw
6. Simulation tree			- stepWiseTree.newick
7. Decay curves calculation script	- r_nonlinear_leastsquare.R
8. Plot result in pdf			- plotPdf.R
9. Mapping between OMA ids and HaMStR id- hamstrMapOmaCorrected.txt

################################################
General directions to use the tool:
################################################

1. In the directory 'toy_example', there exists a program configuration file 'prog.config'. Use this file to set the paths for the various dependencies.
2. Once these dependencies are set, you can simply try running from the toy_example directory the following command:

	python ../bin/protTrace.py -i test.id -c prog.config

3. The program should start running and you can see the outputs being created in:  ../output/YEAST01111

4. One can also provide a fasta file as an input. In that case simply run:

    	python ../bin/protTrace.py -f $yourFile.fasta -c prog.config

   *** NOTE: Use short headers in the fasta file (max. 15 letters) as it will be used to create a folder in the output directory***
  
