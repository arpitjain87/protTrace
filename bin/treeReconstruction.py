import os, sys
import subprocess
import time
import maxLikDistMatrix

# Module to reconstruct tree using RAxML

# Calculate the median of a set
def median(lst):
    even = (0 if len(lst) % 2 else 1) + 1
    half = (len(lst) - 1) / 2
    return sum(sorted(lst)[half:half + even]) / float(even)

# Rename the ortholog group sequences into oma ids
def rename_orth_file():
	fnew = open('temp_orth_%s.fa' %protein_id, 'w')
	for i in range(0, len(orth_file) - 1, 2):
		fnew.write(orth_file[i][:6] + '\n' + orth_file[i+1] + '\n')
	fnew.close()

# Perform MSA of the new renamed ortholog sequences file (MAFFT linsi)
# Perform degapping script by Ingo (degapper.pl)
# Convert the MSA format from .aln to .phy (clustalw2)
def msa_convert():
	#print 'MSA of the renamed ortholog sequences file:'
	os.system('%s temp_orth_%s.fa > temp_orth_%s.aln' %(linsi, protein_id, protein_id))
	#print 'File format conversion (.aln -> .fa):'
	s1 = open('temp_orth_%s.aln' %protein_id).read().split('\n')
	fnew = open('temp_orth_aln_%s.fa' %protein_id, 'w')
	fnew.write(s1[0] + '\n')
	for i in range(1, len(s1) - 1):
		if s1[i][0] == '>':
			fnew.write('\n' + s1[i] + '\n')
		else:
			fnew.write(s1[i])
	fnew.close()
	#print 'Apply degapping algorithm (degapper.pl):'
	os.system('perl %s -limit=0.75 -in=temp_orth_aln_%s.fa -out=temp_orth_aln_degap_%s.fa' %(degap, protein_id, protein_id))
	
	#print 'File format conversion (.fa -> .phy):'
	os.system('%s -convert -infile=temp_orth_aln_degap_%s.fa -outfile=temp_orth_%s.phy -output=PHYLIP' %(clustalw, protein_id, protein_id))
	#print 'complete..'

# Run RAxML for tree reconstruction
def run_raxml():
	#print 'Reconstructing tree with constraint topology:'
	os.system('%s -s temp_orth_%s.phy -m PROTGAMMA%s -p 12345 -n %s' %(raxml, protein_id, aaMatrix, protein_id))
	#print 'complete..'

# Rename the reconstructed tree file
def rename_raxml():
	try:
		s1 = open('RAxML_bestTree.' + protein_id).read()
		fnew = open('RAxML_bestTree.' + protein_id, 'w')
		for i in range(0, len(orth_file) - 1, 2):
			nick_name = orth_file[i][1:6]
			full_name = orth_file[i][1:]
			s1 = s1.replace(nick_name, full_name)
		fnew.write(s1)
		fnew.close()
	except IOError:
		sys.exit('RAxML tree could not be found!')

# Remove all the temp files generated
def rm_temp():
	
	os.system('rm temp_orth_%s.fa' %protein_id)
	os.system('rm temp_orth_%s.aln' %protein_id)
	os.system('rm temp_orth_%s.phy' %protein_id)
	os.system('rm temp_orth_aln_%s.fa' %protein_id)
	os.system('rm temp_orth_aln_degap_%s.fa' %protein_id)
	if os.path.exists('temp_orth_%s.phy.reduced' %protein_id):
		os.system('rm temp_orth_%s.phy.reduced' %protein_id)
	os.system('rm temp_parameters_%s.txt' %protein_id)
	os.system('rm outfile')
	os.system('rm outdist')
	os.system('rm maxLikDist_%s.txt' %protein_id)
	os.system('rm outlm.eps')

# Perform likelihood mapping on the degapped sequences (for 4 or more sequences)
def likelihoodMapping():

	fnew = open('temp_parameters_%s.txt' %protein_id, 'w')
	# Select the parameter file for puzzle i.e. fill in the alignmentFile used !!! *** CHANGE THE FILE HERE ***
	s1 = open(params_puzzle).read()
	fnew.write(s1.replace('alignmentFile', 'temp_orth_%s.phy' %protein_id))
	fnew.close()

	#print 'Running tree puzzle..'
	os.system('%s < temp_parameters_%s.txt' %(puzzle,protein_id))
	
	#print 'Puzzle run complete..'

# Perform likelihood mapping on the degapped sequences (for 3 or less sequences)
def distanceMapping():

	fnew = open('temp_parameters_%s.txt' %protein_id, 'w')
	# Select the parameter file for puzzle i.e. fill in the alignmentFile used !!! *** CHANGE THE FILE HERE ***
	s1 = open(params_puzzle).read()
	fnew.write(s1.replace('alignmentFile', 'temp_orth_%s.phy' %protein_id).replace('b', 'k\nk'))
	fnew.close()

	#print 'Running tree puzzle..'
	os.system('%s < temp_parameters_%s.txt' %(puzzle, protein_id))
	
	#print 'Puzzle run complete..'

# Calculate the scaling factor based on maximum likelihood distances
def scalingFactorMax():
	scales = []
	# Generate maximum likelihood distance file for orthologs
	### NOTE THE FILE USED HERE!!!!*******************************************!!!!
	outfile = open('outfile').read().split('\n')
	maxLikDistMatrix.main(outfile, protein_id)
	try:
		# THESE FILES HAVE TO BE CHANGED WITH NEW SPECIES TREE
		orthMaxFile = open('maxLikDist_%s.txt' %protein_id).read().split('\n')
		speciesMaxFile = open(species_maxLikMatrix).read().split('\n')
		hamstrFile = open(map_file).read().split('\n')
		for i in range(len(orthMaxFile) - 1):
			species1 = orthMaxFile[i].split('\t')[0]
			for j in range(len(hamstrFile) - 1):
				if species1 == hamstrFile[j].split('\t')[3]:
					hamstr1 = hamstrFile[j].split('\t')[0]
					break
			for k in range(i + 1, len(orthMaxFile) - 1):
				species2 = orthMaxFile[k].split('\t')[0]
				for j in range(len(hamstrFile) - 1):
					if species2 == hamstrFile[j].split('\t')[3]:
						hamstr2 = hamstrFile[j].split('\t')[0]
						break
				maxDistOrth = float(orthMaxFile[i].split('\t')[k + 1])
				for l in range(len(speciesMaxFile) - 1):
					if speciesMaxFile[l].split('\t')[0] == hamstr1:
						rowIndex = l
					elif speciesMaxFile[l].split('\t')[0] == hamstr2:
						columnIndex = l + 1
				maxDistSpecies = float(speciesMaxFile[rowIndex].split('\t')[columnIndex])
				if not maxDistSpecies == 0:
					scales.append(maxDistOrth / maxDistSpecies)
				else:
					scales.append(1.00)
	except:
		print '### ERROR: Scaling factor calculation had an error ###'
		sys.exit('Maximum likelihood files are invalid!')
	return median(scales)							

# Main module for running tree reconstruction
def main(Raxml, Linsi, Clustalw, Degap, Orthologs, AaMatrix, Protein_id, Species_tree_msa, Puzzle, Params_puzzle, Map_file, Species_maxLikMatrix, Scale_file, Tree_file, delTemp):

	global raxml, linsi, clustalw, degap, orth_file, aaMatrix, protein_id, species_tree_msa, puzzle, params_puzzle, map_file, species_maxLikMatrix, scaleFile, treeFile
	raxml = Raxml
	linsi = Linsi
	clustalw = Clustalw
	degap = Degap
	aaMatrix = AaMatrix
	protein_id = Protein_id
	species_tree_msa = Species_tree_msa
	puzzle = Puzzle
	params_puzzle = Params_puzzle
	map_file = Map_file
	species_maxLikMatrix = Species_maxLikMatrix
	orth_file = open(Orthologs).read().split('\n')
	scaleFile = Scale_file
	treeFile = Tree_file

	sf = 1.00
		
	if len(orth_file) > 7:
		if not os.path.exists(treeFile):
			try:
				rename_orth_file()
				msa_convert()
				run_raxml()
				likelihoodMapping()
				rename_raxml()
				sf = scalingFactorMax()
			except:
				print '### ERROR: Some step in the tree reconstruction was invalid!! ###'
				pass
			print 'Scaling factor: ', sf
			fnew = open(scaleFile, 'w')
			fnew.write(str(sf))
			fnew.close()
			if delTemp:
				rm_temp()			
		else:
			'WARNING: Tree file already exists!!!'
	elif len(orth_file) < 8 and len(orth_file) > 3:
		try:
			rename_orth_file()
			msa_convert()
			distanceMapping()
			sf = scalingFactorMax()
			
		except:
			pass
		print 'Scaling factor: ', sf
		fnew = open(scaleFile, 'w')
		if sf > 0:
			fnew.write(str(sf))
		else:
			fnew.write('1.00')
		fnew.close()
		if delTemp:
			rm_temp()
	
