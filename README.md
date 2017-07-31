# protTrace
@Author: Arpit Jain
@Email: arpitrb@gmail.com

### protTrace is a simulation-based framework to estimate the evolutionary traceabilities of proteins. ###

################################################
Software dependencies:
################################################

1. Multiple Sequence Alignment (MSA)		- MAFFT
2. Tree reconstruction				- RAxML
3. HMMER tools					- hmmscan, hmmfetch
4. MSA format conversion (.aln -> .phy)		- ClustalW
5. For BLAST					- blastall
6. For protein BLAST				- blastp
7. Formatting BLAST db				- formatdb
8. Creating BLAST db				- makeblastdb
9. Plotting decay curves			- R
10. Orthologs search 				- HaMStR
11. Orthologs search (for single seq)		- HaMStR-OneSeq		 
12. Max. likelihood distances			- TreePuzzle
13. Programming languages			- Python (v2.7 or higher), Perl (v5.22.1 or higher), JAVA (v1.7 or higher)

################################################
Other programs / files (provided in folder 'used_files')
################################################

1. Simulating protein sequence evolution	- REvolver
2. Simulation tree for REvolver			- stepWiseTree.newick
3. TreePuzzle parameters file			- paramsMaxLikelihoodMapping.txt
4. Reference species tree			- speciesTree.nw
5. All vs All Max. likelihood distances Matrix	- speciesLikelihoodMatrix.txt
6. Cross reference file - HaMStR_id <tab> SpeciesName <tab> NCBI_id <tab> OMA_id	-speciesTreeMapping.txt
7. R script to fit non-linear decay curves	- r_nonlinear_leastsquare.R
8. R script to plot traceability results in pdf	- plotPdf.R
9. Tool to estimate insertions/deletions	- iqtree-24

################################################
General directions to use the tool:
################################################

1. In the directory 'toy_example', there exists a program configuration file 'prog.config'. Use this file to set the paths for the various dependencies.
2. Once these dependencies are set, you can simply try running from the toy_example directory the following command(s):

	python ../bin/protTrace.py -i test.id -c prog.config
	
	python ../bin/protTrace.py -f test.fasta -c prog.config

*** NOTE: In the first example, traceabilities will be calculated for yeast id YEAST01111 obtained from OMA. If you have OMA id for your protein of interest, simply run using this command. For cases where you do not have OMA id, run the second command. Here, you can give the protein sequence of interest. Please use short headers in the fasta file (max. 15 letters) as it will be used for naming the results folder in the output directory. Also, make sure that the query sequnce is in a single line i.e. multi-lined fasta files will not be accepted. ***

3. The program should run and you can see the output files being generated in output directories: ../output/YEAST01111 and ../output/yeast_MAPK respectively


  
  
