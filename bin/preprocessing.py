#
# Module to perform preprocessing steps in protTrace workflow
#
import os, sys
import glob
import configure
import hamstr_search
import treeReconstruction
import transformAlignment
import dendropy
import subprocess
import random
import time

def Preprocessing(prot_id, querySeq, config_file):
	# Store the current working directory into a variable
	rootDir = os.getcwd()

	print 'Prot_id: ', prot_id

	# Creating instance of class configure
	# Saves the information provided by the program configuration file
	prot_config = configure.setParams(config_file)

	# Declares global variables which will be used by all methods in the module
	global cache, cache_dir, work_dir, omaIdFile, orth_file, aln_file, phy_file, id_file, proteome_file, tree_file, trans_file, hmm_file, xml_file, REvolver_output_dir, species_id, indel_file, scale_file

	# Getting information from the configuration file
	# Setting the names for the protTrace temporary and output files	
	species_id = prot_config.species
	proteome_file = 'proteome_' + prot_id
	id_file = 'ogIds_' + prot_id + '.txt'
	work_dir = prot_config.path_work_dir + '/' + prot_id
	cache_dir = prot_config.path_cache
	omaIdFile = work_dir + '/omaId.txt'
	orth_file = 'ogSeqs_' + prot_id + '.fa'
	aln_file = 'ogSeqs_' + prot_id + '.aln'
	phy_file = 'ogSeqs_' + prot_id + '.phy'
	tree_file = 'RAxML_bestTree.' + prot_id 	###	CHANGE HERE IF RAxML OUTPUT NAMES CHANGE 	###
	trans_file = 'ogSeqs_' + prot_id + '.trans'
	hmm_file = prot_id + '.hmm'
	xml_file = 'revolver_config_' + prot_id + '.xml'
	REvolver_output_dir = work_dir + '/REvolver_output/'
	indel_file = 'indel_' + prot_id
	scale_file = 'scale_' + prot_id
	delTemp = prot_config.delete_temp
	cache = prot_config.reuse_cache

	# Creates a working directory where temporary and output files will be stored
	if not os.path.exists(work_dir):
		print '##### Creating working directory:\n', work_dir
		try:
			os.mkdir(work_dir)
		except:
			sys.exit('ERROR: Working directory cannot be created!')

	# Change current working directory
	os.chdir(work_dir)

	# Parse proteome of the input species given in program configuration file
	# The proteome is extracted from the OMA database sequences file
	if prot_config.search_proteome:
		startProcessTime = time.time()
		parseOmaProteome(species_id, prot_config.path_oma_seqs, prot_config.makeblastdb, proteome_file)
		proteome_file = os.path.abspath(proteome_file)
		print '#####\tTIME TAKEN: %s mins\tSearch proteome#####' %((time.time() - startProcessTime) / 60) 	

	# Search ortholog groups for the input OMA id
	# In case of fasta sequences, first the OMA id is parsed from the OMA database
	# followed with the extraction of the OMA group
	if prot_config.search_ortholog_groups:
		startProcessTime = time.time()		
		run = findOmaOrthologs(prot_id, querySeq, prot_config.path_oma_group, prot_config.path_oma_seqs, proteome_file, prot_config.formatdb, prot_config.blastp, delTemp)
		print '#####\tTIME TAKEN: %s mins\tSearch OGs#####' %((time.time() - startProcessTime) / 60)

	# Search for the ortholog sequences for the respective OMA orthologs group
	# For all the OMA ids in the OMA group, extract sequences from OMA database sequences file
	if prot_config.search_ortholog_sequences:
		if cache and os.path.exists(orth_file):
			print 'OMA orthologs sequences file exist. Reusing it!'
		else:
			if run == 2:
				startProcessTime = time.time()
				findOmaSequences(prot_id, prot_config.path_oma_seqs, species_id, prot_config.hamstr_oma_tree_map)
				print '#####\tTIME TAKEN: %s mins\tSearch OGSeqs#####' %((time.time() - startProcessTime) / 60)
			else:
				print '##### Preparing ortholog file #####'
				fOrth = open(orth_file, 'w')
				fOrth.write('>' + species_id + '\n' + querySeq)
				fOrth.close()
			

	# Extend the ortholog set by performing a HaMStR search
	# The orthologs sequences by OMA is used as core-ortholog set
	try:
		f = open(orth_file).read().split('\n')
		# Run HaMStR search if 2 or more sequences are present. Otherwise, run HaMStROneSeq search if only 1 sequence is present
		if len(f) > 3:
			if prot_config.run_hamstr:
				startProcessTime = time.time()
				if cache and os.path.exists(cache_dir + '/' + orth_file):
					print 'Pre-HaMStR computed orthologs file found in Cache for re-use!'
					os.system('cp %s %s' %(cache_dir + '/' + orth_file, orth_file))
				else: 
					success = hamstr_search.main(prot_config.hamstr, orth_file, prot_id, prot_config.hamstr_oma_tree_map, prot_config.formatdb, prot_config.blastp, delTemp)
					if success:
						print '#####\tTIME TAKEN: %s mins\tHaMStR#####' %((time.time() - startProcessTime) / 60)
						os.system('cp %s %s' %(orth_file, cache_dir + '/' + orth_file))
					else:
						if prot_config.run_hamstrOneSeq:
							print '##### HaMStROneSeq search for orthologs #####'
							startProcessTime = time.time()
							if cache and os.path.exists(cache_dir + '/' + orth_file):
								print 'Pre-HaMStR computed orthologs file found in Cache for re-use!'
								os.system('cp %s %s' %(cache_dir + '/' + orth_file, orth_file))
							else:
								# Read the orthologs file and limit it to just the query species id and sequence
								ortholog_temp = open(orth_file).read().split('\n')
								rewrite_orth_file = open(orth_file, 'w')
								for orthLines in range(len(ortholog_temp) - 1):
									if '>' in ortholog_temp[orthLines] and species_id in ortholog_temp[orthLines]:
										rewrite_orth_file.write(ortholog_temp[orthLines] + '\n' + ortholog_temp[orthLines + 1])
										break
								rewrite_orth_file.close()
								run_hamstrOneSeq(prot_config.hamstr, os.path.abspath(orth_file), prot_config.hamstr_oma_tree_map, prot_id, prot_config.formatdb, prot_config.blastp, proteome_file, delTemp)
								print '#####\tTIME TAKEN: %s mins\tHaMStR-OneSeq#####' %((time.time() - startProcessTime) / 60)
								os.system('cp %s %s' %(orth_file, cache_dir + '/' + orth_file))

			elif prot_config.run_hamstrOneSeq:
				startProcessTime = time.time()
				print '##### HaMStROneSeq search for orthologs #####'
				if cache and os.path.exists(cache_dir + '/' + orth_file):
					print 'Pre-HaMStR computed orthologs file found in Cache for re-use!'
					os.system('cp %s %s' %(cache_dir + '/' + orth_file, orth_file))
				else:
					# Read the orthologs file and limit it to just the query species id and sequence
					ortholog_temp = open(orth_file).read().split('\n')
					rewrite_orth_file = open(orth_file, 'w')
					for orthLines in range(len(ortholog_temp) - 1):
						if '>' in ortholog_temp[orthLines] and species_id in ortholog_temp[orthLines]:
							rewrite_orth_file.write(ortholog_temp[orthLines] + '\n' + ortholog_temp[orthLines + 1])
							break
					rewrite_orth_file.close()
					run_hamstrOneSeq(prot_config.hamstr, os.path.abspath(orth_file), prot_config.hamstr_oma_tree_map, prot_id, prot_config.formatdb, prot_config.blastp, proteome_file, delTemp)
					print '#####\tTIME TAKEN: %s mins\tHaMStR-OneSeq#####' %((time.time() - startProcessTime) / 60)
					os.system('cp %s %s' %(orth_file, cache_dir + '/' + orth_file))
			
			else:
				print '#####\tCalculating traceability with ONLY OMA set\t#####'
				

		elif len(f) > 0 and len(f) < 4:
			if prot_config.run_hamstrOneSeq:
				print '##### HaMStROneSeq search for orthologs #####'
				startProcessTime = time.time()
				if cache and os.path.exists(cache_dir + '/' + orth_file):
					print 'Pre-HaMStR computed orthologs file found in Cache for re-use!'
					os.system('cp %s %s' %(cache_dir + '/' + orth_file, orth_file))
				else:
					run_hamstrOneSeq(prot_config.hamstr, os.path.abspath(orth_file), prot_config.hamstr_oma_tree_map, prot_id, prot_config.formatdb, prot_config.blastp, proteome_file, delTemp)
					print '#####\tTIME TAKEN: %s mins\tHaMStR-OneSeq#####' %((time.time() - startProcessTime) / 60)
					os.system('cp %s %s' %(orth_file, cache_dir + '/' + orth_file))

		else:
			sys.exit('ERROR: No sequence found in OMA sequences! The ortholog sequences file is empty!')
	except IOError:
		sys.exit('ERROR: Orthologs sequences file is invalid!')

	except KeyboardInterrupt:
		sys.exit('Keyboard interruption by user!!!')

	# Performs MSA on the orthologs sequences
	if prot_config.perform_msa:
		print '##### Performing MSA of the orthologs sequences #####'
		startProcessTime = time.time()
		performMSA(prot_config.msa, prot_config.clustalw)
		print '#####\tTIME TAKEN: %s mins\tMAFFT#####' %((time.time() - startProcessTime) / 60)

	# Calls tree reconstruction module which generates tree using degapped alignment
	# and also calculates the scaling factor based on maximum likelihood distance between species
	if prot_config.calculate_scaling_factor:
		print '##### Tree reconstruction and scaling factor calculation #####'
		startProcessTime = time.time()
		treeReconstruction.main(prot_config.tree_reconstruction, prot_config.msa, prot_config.clustalw, prot_config.degapping, orth_file, prot_config.aa_substitution_matrix, prot_id, prot_config.treePuzzle, prot_config.parameters_treePuzzle, prot_config.hamstr_oma_tree_map, prot_config.species_MaxLikMatrix, scale_file, tree_file, delTemp, prot_config.default_scaling_factor, cache_dir)
		print '#####\tTIME TAKEN: %s mins\tRAxML#####' %((time.time() - startProcessTime) / 60)

	# Calculate indels
	if prot_config.calculate_indel:
		if cache and os.path.exists(indel_file):
			print 'Pre-computed indel found for re-use!'
		else:

			# Transform alignment
			print '##### Transforming MSA based on indel blocks #####'
			alignmentLength = 0
			try:
				alignmentLength = transformAlignment.main(phy_file, trans_file)
			except:
				pass
			calculateIndels(tree_file, trans_file, alignmentLength, prot_config.iqtree24, prot_config.default_indel, prot_config.default_indel_distribution)

	# Domain constraint file for REvolver
	if prot_config.traceability_calculation:
		# Creates a output directory for REvolver
		if not os.path.exists(REvolver_output_dir):
			print '##### Creating REvolver output directory:\n', REvolver_output_dir
			try:
				os.mkdir(REvolver_output_dir)
			except:
				sys.exit('ERROR: REvolver output directory cannot be created!')

		print '##### Generating domain constraints for REvolver #####'
		hmmscan(prot_config.hmmscan, orth_file, prot_config.pfam_database, hmm_file, prot_id, species_id)

		# Prepare XML config file to be used as an input for REvolver
		print '##### Preparing XML configuration file for REvolver #####'
		if os.path.exists(scale_file):
			f = open(scale_file).read().split('\n')
			scaling_factor = f[0]
		else:
			print 'WARNING: Scaling factor file not found. Using default value:', prot_config.default_scaling_factor
			scaling_factor = prot_config.default_scaling_factor
			# Writing randomly generated indels into the file (Randome pick between 0.50 to 2.00)
			writeScale = open(scale_file, 'w')
			writeScale.write(scaling_factor)
			writeScale.close()

		if os.path.exists(indel_file):
			f = open(indel_file).read().split('\n')
			indel = f[0]
			p = f[1]
		else:
			print 'WARNING: Indel file not found. Using default value:', prot_config.default_indel
			indel = prot_config.default_indel
			p = prot_config.default_indel_distribution
			# Writing randomly generated indels into the file (Random pick for indel between 0.05 to 1.9; for distribution between 0.0 to 0.9)
			# The random value limit has been set as observed from the large scale indels distribution
			writeIndel = open(indel_file, 'w')
			writeIndel.write(indel + '\n' + p)
			writeIndel.close()

		prepareXML(xml_file, prot_config.pfam_database, prot_config.hmmfetch, prot_config.aa_substitution_matrix, indel, p, scaling_factor, prot_config.simulation_tree, prot_id, hmm_file, REvolver_output_dir)
	
	os.chdir(rootDir)

# HaMStROneSeq run
def run_hamstrOneSeq(hamstr, orth_file, map_file, prot_id, formatdb, blastp, proteome, delTemp):

	# Setting default paths for HaMStR-OneSeq
	hamstrOneSeq = hamstr + '/bin/oneSeq.pl'
	try:
		coreOrthDir = hamstr + '/core_orthologs'
		taxaPath = hamstr + '/genome_dir'

		for line in open(orth_file):
			if line[0] == '>':
				omaId = line.split('\n')[0][1:6]
			else:
				protSeq = line.split('\n')[0]
		
		for line in open(map_file):
			if line.split('\t')[3].split('\n')[0] == omaId:
				hamstrId = line.split('\t')[0]
				break
		print 'hamstr id: ', hamstrId
		for taxas in glob.glob(taxaPath + '/*'):
			if taxas.split('/')[-1].split('_')[0] + '_' + taxas.split('/')[-1].split('_')[1].split('@')[0] == hamstrId:
				refSpec = taxas.split('/')[-1]
				break
		
		with open(taxaPath + '/' + refSpec + '/' + refSpec + '.fa') as f:
			print 'Searching for the seqId..'
			flag = True
			for line in f:
				if f.next().split('\n')[0].replace('*', '') == protSeq.replace('*', ''):
					seqId = line.split('>')[1].split('\n')[0]
					flag = False
					break
		if flag:
			seqId = "NA"
			print 'No matching sequence found! Running BLAST search now..'
			#print 'Current working directory..'
			currentWorkDir = os.getcwd()
			#print 'Create a temporary directory..'
			tempDir = 'temp_blast_' + prot_id
			if not os.path.exists(tempDir):
				os.mkdir(tempDir)
			#print 'Change to the temporary directory..'
			os.chdir(tempDir)
			#print 'Create a temporary file with the input sequence..'
			#print 'Copy the reference proteome file into temporary directory..'
			os.system('cp -avr %s .' %(taxaPath + '/' + refSpec + '/' + refSpec + '.fa'))
			com = '%s -i %s' %(formatdb, refSpec +'.fa')
			#print 'Create blast database for the OMA sequences: ', com
			os.system(com)
			#print 'Perform BLAST search and pick up the top hit as query input ID..'
			os.system('%s -query %s -db %s -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out temp.txt' %(blastp, orth_file, refSpec + '.fa'))
			#os.system('%s -query %s -db %s -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out temp.txt' %(blastp, orth_file, proteome))
			if os.path.exists('temp.txt') and len(open('temp.txt').read().split('\n')) > 1:
				seqId = open('temp.txt').read().split('\n')[0].split('\t')[1]
			os.system('%s -query %s -db %s -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out temp.txt' %(blastp, orth_file, proteome))
			if os.path.exists('temp.txt') and len(open('temp.txt').read().split('\n')) > 1:
				omaId = open('temp.txt').read().split('\n')[0].split('\t')[1]
			#print 'Writing OMA id into file'
			oma = open(omaIdFile, 'w')
			oma.write(omaId)
			oma.close()
			#print 'Change directory to original one..'
			os.chdir(currentWorkDir)
			#print 'Remove the temporary directory..'
			if delTemp:
				os.system('rm -rf %s' %tempDir)

		print 'SeqId: %s  ...  OmaId: %s' %(seqId, omaId)
		try:
			if not seqId == "NA" and delTemp:
				command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -fasoff -local -strict -rep -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
				#command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -cleanup -fasoff -local -strict -rep -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
				print '##### Running hamstrOneSeq command: ', command
				os.system(command)

			elif not seqId == "NA" and not delTemp:
				command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -fasoff -local -strict -rep -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
				#command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -fasoff -local -strict -rep -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
				print '##### Running hamstrOneSeq command: ', command
				os.system(command)

			else:
				print 'WARNING: No sequence id found to run HaMStR-OneSeq!!!'
		except:
			pass
			print 'HaMStROneSeq did not run properly!!!'
	except IOError:
		print 'WARNING: hamstrOneSeq did not run properly!!!'


	# Read the output files generated by OneSeq run (.fa and .fa.extended)
	# Write the results back to the original orthologs file
	print 'Adding HaMStR-OneSeq results to the orthologs file..'
	try:
		orthFileDefault = open(orth_file).read().split('\n')
		speciesList = []
		speciesList.append(orthFileDefault[0][1:6])
		fnew = open(orth_file, 'w')
		fnew.write(orthFileDefault[0] + '\n' + orthFileDefault[1] + '\n')

		#faFile = coreOrthDir + prot_id + '/%s.fa' %prot_id
		extendedFile = work_dir + '/%s.extended.fa' %prot_id
		'''if os.path.exists(faFile):
			omaId = "NA"
			with open(faFile) as f:
				for line in f:
					if line[0] == '>':
						hamstrId = line.split('|')[1].split('_')[0] + '_' + line.split('|')[1].split('_')[1]
						for m in open(map_file):
							if m.split('\t')[0] == hamstrId:
								omaId = m.split('\t')[3].split('\n')[0]
								break
						if not omaId == "NA" and omaId not in speciesList:
							speciesList.append(omaId)
							#fnew.write('>' + omaId + line.split('|')[2].split('\n')[0] + '\n' + f.next().replace('*', ''))
							fnew.write('>' + omaId + '\n' + f.next().replace('*', ''))
			#fnew.write('\n')'''
		if os.path.exists(extendedFile):
			#print 'Found .extended file..'
			omaId = "NA"
			with open(extendedFile) as f:
				for line in f:
					if '>' in line:
						hamstrId = line.split('|')[1].split('_')[0] + '_' + line.split('|')[1].split('_')[1].split('@')[0]
						for m in open(map_file):
							if m.split('\t')[0] == hamstrId:
								omaId = m.split('\t')[3].split('\n')[0]
								break
						if not omaId == "NA" and omaId not in speciesList:
							speciesList.append(omaId)
							#fnew.write('>' + omaId + line.split('|')[2].split('\n')[0] + '\n' + f.next().replace('*', ''))
							fnew.write('>' + omaId + '\n' + f.next().replace('*', ''))
		fnew.close()
	except IOError:
		print 'ERROR: While writing HaMStr-OneSeq results!!!'

	if delTemp:
		if os.path.exists(coreOrthDir + '/' + prot_id):
			os.system('rm -Rf %s' %(coreOrthDir + '/' + prot_id))


# Prepares input configuration file for REvolver
def prepareXML(xml_file, pfamDB, hmmfetch, aaMatrix, indel, p, sf, simTree, prot_id, hmm_file, output_dir):
	fnew = open(xml_file, 'w')
	fnew.write('<?xml version="1.0" encoding="UTF-8" ?>\n')
	fnew.write('<configdata  xsi:schemaLocation="http://www.cibiv.at/Revolver ./input_schema.xsd" xmlns="http://www.cibiv.at/Revolver" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" >\n')
	fnew.write('\t<config>\n')
	fnew.write('\t\t<hmmdb path="%s"/>\n' %os.path.abspath(pfamDB))
	fnew.write('\t\t<hmmfetch location="%s"/>\n' %hmmfetch)
	fnew.write('\t</config>\n')
	fnew.write('\t<model>\n')
	fnew.write('\t\t<substitution name="%s"/>\n' %aaMatrix)
	fnew.write('\t\t<indel>\n')
	fnew.write('\t\t\t<insertion rate="%s">\n' %str(indel))
	fnew.write('\t\t\t\t<length distribution="geometric" p="%s"/>\n' %str(p))
	fnew.write('\t\t\t</insertion>\n')
	fnew.write('\t\t\t<deletion rate="%s">\n' %str(indel))
	fnew.write('\t\t\t\t<length distribution="geometric" p="%s"/>\n' %str(p))
	fnew.write('\t\t\t</deletion>\n')
	fnew.write('\t\t</indel>\n')
	fnew.write('\t</model>\n')
	fnew.write('\t<tree scalingFactor="%s" path="%s"  />\n' %(str(sf), os.path.abspath(simTree)))
	fnew.write('\t<root>\n')
	fnew.write('\t\t<inputSequence>\n')
	fnew.write('\t\t\t<fasta file="seq_%s.fa"/>\n' %prot_id)
	fnew.write('\t\t\t<hmmer file="%s"/>\n' %os.path.abspath(hmm_file))
	fnew.write('\t\t</inputSequence>\n')
	fnew.write('\t</root>\n')
	fnew.write('\t<output>\n')
	fnew.write('\t\t<dir path="%s" separateFastaFiles="false" trueAlignment="false" include="leaf"/>\n' %os.path.abspath(output_dir))
	fnew.write('\t</output>\n')
	fnew.write('</configdata>')
	
	fnew.close()

# Runs hmmscan and prepares the domain constraint file for REvolver
def hmmscan(hmmscan, orth_file, pfamDB, hmm_file, prot_id, species_id):
	ftemp = open('seq_%s.fa' %prot_id, 'w')
	with open(orth_file) as f:
		for line in f:
			if line[0] == '>' and line.split('\n')[0][1:6] == species_id:
				ftemp.write(line + f.next())
				break
	ftemp.close()
	os.system('%s --notextw -E 0.01 %s seq_%s.fa > %s' %(hmmscan, pfamDB, prot_id, hmm_file))
	#os.remove('tempFile.fa')
	

# Calculates indels rates
def calculateIndels(tree_file, trans, alnLength, iqtree24, def_indel, def_indel_dist):
	indel = float(def_indel)
	p = float(def_indel_dist)

	print '##### Calculating indels #####'
	try:
		trees = dendropy.TreeList.get_from_path(tree_file, "newick")
		tree_lengths = [tree.length() for tree in trees]
	except:
		pass
	
	result = ''
	try:
		command = "%s -s %s %s -tina -st MULTI" %(iqtree24, os.path.abspath(trans), os.path.abspath(tree_file))
		print 'IQ-Tree24 command: ', command
		result = subprocess.check_output(command, shell=True)
	except:
		print 'WARNING: IQTree-24 did not run properly!!!'
		pass
	
	for line in result.split('\n'):
		if line.split(':')[0] == 'mean length':
			if float(line.split(':')[1].replace(' ', '')) > 0:
				p = 1 / float(line.split(':')[1].replace(' ', ''))
				if p >= 1:
					p = 0.99
				elif p < 0.02:
					p = 0.02
			else:
				p = float(prot_config.default_indel_distribution)
		elif line.split(':')[0] == 'Parsimony score is':
			indel = (float(line.split(':')[1].replace(' ', '')) / (alnLength * tree_lengths[0])) / 2
	print 'Indel: ', indel

	fnew = open(indel_file, 'w')
	fnew.write(str(indel) + '\n' + str(p))
	fnew.close()
	

# Perform MSA of the ortholog sequences
# Convert the .aln format to .phy format
def performMSA(msa, clustalw):
	if cache:
		if os.path.exists(phy_file):
			print 'Re-using previously compiled orthologs alignment file'
		else:	
			try:
				os.system('%s %s > %s' %(msa, orth_file, aln_file))
				os.system('%s -convert -output=PHYLIP -infile=%s -outfile=%s' %(clustalw, aln_file, phy_file))
			except:
				pass
				print 'WARNING: MSA didn\'t work. Less than 2 sequences found for alignment!!!'
	else:
		try:
			os.system('%s %s > %s' %(msa, orth_file, aln_file))
			os.system('%s -convert -output=PHYLIP -infile=%s -outfile=%s' %(clustalw, aln_file, phy_file))
		except:
			pass
			print 'WARNING: MSA didn\'t work. Less than 2 sequences found for alignment!!!'

# Read OMA sequences file and parse OMA orthologs sequences
def findOmaSequences(prot_id, omaSeqs, species_id, mapFile):
	try:
		mapIds = []
		for line in open(mapFile):
			mapIds.append(line.split()[-1])

		print '##### Searching OMA orthologs sequences for %s #####' %prot_id
		fnew = open(orth_file, 'w')
		ids = open(id_file).read().split('\n')
		with open(omaSeqs) as f:
			for line in f:
				if line[0] == '>' and line.split('\n')[0][1:] in ids and line[1:6] in mapIds:
					#if line.split('\n')[0][1:6] == species_id:
						#fnew.write('>' + species_id + '\n' + f.next().replace('*', ''))
					fnew.write('>' + line[1:6] + '\n' + f.next().replace('*', '').replace('X', ''))
					#else:
					#	fnew.write(line + f.next().replace('*', ''))
		fnew.close()
	except IOError:
		sys.exit('ERROR: Cannot find OMA orthologs sequences. OMA sequence file does not exist!')
				
	
# Read OMA orthologs groups file and parses the ortholog list for input OMA id
def findOmaOrthologs(prot_id, querySeq, omaGroup, omaSeqs, proteome_file, formatdb, blastp, delTemp):
	try:
		if not querySeq == 'None':
			run = 1
			flag = False
			with open(omaSeqs) as f:
				for line in f:
					if '>' in line:
						if f.next().split('\n')[0].replace('*', "") == querySeq.split('\n')[0].replace('*', ''):
							prot_id_temp = line.split('\n')[0][1:]
							oma = open(omaIdFile, 'w')
							oma.write(prot_id_temp)
							oma.close()
							#print 'Found OMA id:', inputId
							flag = True
							#fnew.write(inputId + '\n')
							#fnew.close()
							break

			if not flag:
				print '##### WARNING: The input sequence is not present in OMA sequences. Performing a BLAST search to get the closest best hit. #####'
				currentWorkDir = os.getcwd()
				tempDir = 'temp_blast_' + prot_id
				if not os.path.exists(tempDir):
					os.mkdir(tempDir)
				os.chdir(tempDir)
				os.system('cp -avr %s proteome.fa' %(proteome_file))
				com = '%s -i proteome.fa' %(formatdb)
				os.system(com)
				tempFile = open('temp_inputSeq.fa', 'w')
				tempFile.write('>' + prot_id + '\n' + querySeq)
				tempFile.close()
				os.system('%s -query temp_inputSeq.fa -db proteome.fa -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out temp.txt' %(blastp))
				if os.path.exists('temp.txt') and len(open('temp.txt').read().split('\n')) > 1:
					prot_id_temp = open('temp.txt').read().split('\n')[0].split('\t')[1]
					oma = open(omaIdFile, 'w')
					oma.write(prot_id_temp)
					oma.close()
					flag = True
				os.chdir(currentWorkDir)
				if delTemp:
					os.system('rm -rf %s' %tempDir)

			if flag:
				run = 2
				print '##### Searching for OMA ortholog group for %s #####' %prot_id
				fnew = open(id_file, 'w')
				written = False
				for line in open(omaGroup):
					if not '#' in line and prot_id_temp in line:
						for ids in line.split('\n')[0].split('\t')[2:]:
							fnew.write(ids + '\n')
						fnew.close()
						written = True
						break
				if not written:
					run = 1
					fnew.write(prot_id_temp)
					fnew.close()
		else:
			print 'OMA ids given..'
			run = 2
			print '##### Searching for OMA ortholog group for already given OMA id %s #####' %prot_id
			oma = open(omaIdFile, 'w')
			oma.write(prot_id)
			oma.close()
			fnew = open(id_file, 'w')
			written = False
			for line in open(omaGroup):
				if not '#' in line and prot_id in line:
					for ids in line.split('\n')[0].split('\t')[2:]:
						fnew.write(ids + '\n')
					fnew.close()
					written = True
					break
			if not written:
				fnew.write(prot_id)
				fnew.close()
			
	except IOError:
		sys.exit('ERROR: Cannot find OMA orthologs id. Check OMA files given as input!')

	return run

# Read in OMA sequences file and create a new proteome file for the species id
# Performs formatdb on the proteome to be later used in reciprocal BLAST search
def parseOmaProteome(species_id, omaSeqs, makeblastdb, proteome_file):
	try:
		# Check whether the parsed proteome in already present in Cache directory
		if cache and os.path.exists(cache_dir + '/proteome_' + species_id):
			print '##### Proteome found in Cache directory #####'
			os.system('cp -avr %s %s' %(cache_dir + '/proteome_' + species_id, proteome_file))
		else:
			print '##### Parsing proteome for species %s in OMA database #####' %species_id
			fnew = open(proteome_file, 'w')
			with open(omaSeqs) as f:
				for line in f:
					if line[:6] == '>' + species_id:
						fnew.write(line + f.next())
			fnew.close()
			os.system('cp -avr %s %s' %(proteome_file, cache_dir + '/proteome_' + species_id))
		print '##### Making BLAST db of proteome to be used by the blast search #####'
		os.system('%s -in %s -input_type fasta -dbtype prot' %(makeblastdb, proteome_file))
		#proteome_file = os.path.abspath(proteome_file)

	except IOError:
		sys.exit('ERROR: Cannot create proteome for species %s. OMA sequences file does not exist!' %species_id)
