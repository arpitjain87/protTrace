import os, sys
import glob
import subprocess

# Script to run HAMSTR search to add orthologs to the core-orthologs groups derived from OMA database
# Input: 1.Path to HAMSTR folder
# 	 2.Path to orthologs group fasta file
#	 3.Protein id
#	 4.Mapping file -> HAMSTR id to OMA id
# Output: Adds orthologs found to the core-ortholog file
# NOTE: The $path and $output variables in the hamstr.pl file has been edited. 

# function to remove if pre-existing run data for a protein id is present
def remove_old_dir():
	if os.path.exists(hamstr + '/core_orthologs/' + protein_id):
		os.system('rm -Rf %s/core_orthologs/%s' %(hamstr, protein_id))
	if os.path.exists(hamstr + '/output/' + protein_id):
		os.system('rm -Rf %s/output/%s' %(hamstr, protein_id))

# function to make new directories for HAMSTR search
def make_new_dir():
	
	os.system('mkdir %s' %core_ortholog_prot_dir)
	os.system('mkdir %s/aln_dir' %core_ortholog_prot_dir)
	os.system('mkdir %s/hmm_dir' %core_ortholog_prot_dir)
	os.system('mkdir %s' %output)

# Function to read the core orthologs file and make a copy of it in the HAMSTR folder
# Edits the fasta file how in format required for the HAMSTR run
def copy_edit_core_ortholog():
	global fasta_file
	fasta_file = core_ortholog_prot_dir + '/' + protein_id + '.fa'
	s1 = open(core_ortholog).read().split('\n')
	fnew = open(fasta_file, 'w')
	flag = True ###	Check whether HaMStR has a starting sequence	###
	for i in range(0, len(s1)-1, 2):
		if s1[i][0] == '>':
			sequence = s1[i+1].replace('*', '').replace('X', '')
			s2 = open(hamstr_map_oma).read().split('\n')
			for j in range(len(s2)-1):
				if s2[j].split()[-1] == s1[i][1:6]:
					hamstrProtId = s2[j].split()[0]
					blastDir = hamstr + '/blast_dir'
					for dirs in glob.glob(blastDir + '/*'):
						if hamstrProtId in dirs:
							blastFastaFile = dirs + '/' + dirs.split('/')[-1] + '.fa'
							break
					s3 = open(blastFastaFile).read().split('\n')
					runBlast = True # If the sequence is not found in our database then, take the top blast hit #
					for k in range(len(s3) - 1):
						if sequence == s3[k].replace('*', '').replace('X', ''):
							fnew.write(s1[i][0] + protein_id + '|' + dirs.split('/')[-1] + '|' + s3[k-1][1:] + '\n' + sequence + '\n')
							flag = False
							runBlast = False
							break

					if runBlast:
						print 'No matching sequence found in hamstr blast directory! Running BLAST search now..'
						currentWorkDir = os.getcwd()
						tempDir = 'temp_blast_' + hamstrProtId
						if not os.path.exists(tempDir):
							os.mkdir(tempDir)
						os.chdir(tempDir)
						temp_query = open('temp_query.fa', 'w')
						temp_query.write(s1[i] + '\n' + s1[i+1])
						temp_query.close()
						os.system('cp -avr %s temp_proteome.fa' %(blastFastaFile))
						os.system('%s -i temp_proteome.fa' %(format_db))
						os.system('%s -query temp_query.fa -db temp_proteome.fa -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out temp_out.txt' %(blastp))

						if os.path.exists('temp_out.txt') and len(open('temp_out.txt').read().split('\n')) > 1:
							hit_id = open('temp_out.txt').read().split('\n')[0].split('\t')[1]
							fnew.write(s1[i][0] + protein_id + '|' + dirs.split('/')[-1] + '|' + hit_id + '\n' + sequence + '\n')
						os.chdir(currentWorkDir)

						if delTemp:
							os.system('rm -rf %s' %tempDir)
					break
	fnew.close()

	if flag:
		return False
	else:
		return True

# Function to run mafft and hmmbuild on the edited core ortholog fasta file
def align_hmm_run():
	global aln_file, hmm_file
	aln_file = core_ortholog_prot_dir + '/aln_dir/' + protein_id + '.aln'
	hmm_file = core_ortholog_prot_dir + '/hmm_dir/' + protein_id + '.hmm'

	f = open(fasta_file).read().split('\n')
	if len(f) > 3:
		os.system('linsi %s > %s' %(fasta_file, aln_file))
		os.system('hmmbuild %s %s' %(hmm_file, aln_file))
	else:
		os.system('hmmbuild %s %s' %(hmm_file, fasta_file))

# Function to run HAMSTR
def hamstr_run():
	species_in_tree = []
	s2 = open(hamstr_map_oma).read().split('\n')
	for i in range(len(s2) - 1):
		species_in_tree.append(s2[i].split('\t')[0])

	for dirs in glob.glob(hamstr + '/genome_dir/*/'):
		#print dirs
		#if not ".fa" in dirs and not ".sql" in dirs:
		hamstr_name = dirs.split('/')[-1]
		if hamstr_name in species_in_tree:
			file_name = dirs + '/' + dirs.split('/')[-1] + '.fa'
			command = '%s/bin/hamstr.pl -central -sequence_file=%s -taxon=misc -hmmset=%s -strict -checkCoorthologsRef -representative -outpath=%s -hit_limit=10' %(hamstr, file_name, protein_id, output)
			#command = '%s/bin/hamstr.pl -central -sequence_file=%s -taxon=misc -hmmset=%s -strict -representative -outpath=%s -hit_limit=10' %(hamstr, file_name, protein_id, output)
			print 'HaMStR run with command: ', command
			try:
				subprocess.call(command, shell=True)
			except KeyboardInterrupt as e:
				sys.exit(e)
			except:
				pass

# Function to add orthologs to the core orthologs set
def write_output():
	count = 0
	for f2 in glob.glob(output+'/*.out'):
		count += 1
	
	if count > 0:
		s1 = open(core_ortholog).read().split('\n')
		fnew = open(core_ortholog, 'w')

		s4 = open(hamstr_map_oma).read().split('\n')
		present_species = []

		for i in range(len(s4) - 1):
			present_species.append(s4[i].split('\t')[0])
		
		co = [] #Variable to store the name of the core-orthologs species (eg -> YEAST00011 would be stored in the list as 'YEAST')

		# Re-writing the original core orthologs set
		for i in range(len(s1)-1):
			if s1[i][0] == '>':
				co.append(s1[i][1:6])
			fnew.write(s1[i] + '\n')
		#print co
		# Reading the .out files and adding the sequences to the core orthologs set
		for files in glob.glob(output+'/*.out'):
			#print 'Writing file: ', files
			s2 = open(files).read().split('\n')

			if '.strict' in files.split('/')[-1]:
				hamstr_name = '_'.join(files.split('/')[-1].split('.strict.out')[0].split('_')[1:])
			else:
				hamstr_name = '_'.join(files.split('/')[-1].split('.out')[0].split('_')[1:])
			s3 = open(hamstr_map_oma).read().split('\n')
			
			oma_name = ''
			for j in range(len(s3)-1):
				if s3[j].split()[0] == hamstr_name:
					#oma_name = s3[j].split()[3] + s2[0].split('|')[3][:5]
					oma_name = s3[j].split()[-1]
					break
			#print oma_name[:5]
			if not oma_name[:5] in co and hamstr_name in present_species:
				fnew.write('>' + oma_name + '\n' + s2[0].split('|')[-1].replace('X', '').replace('*', '') + '\n')
		fnew.close()
	else:
		print 'No orthologs to be added!!'

# Edit the orthologs file to just contain species which are present in the species tree
def finalEditOrthologsSeqs():
	s1 = open(core_ortholog).read().split('\n')
	fnew = open(core_ortholog, 'w')
	
	species_in_tree = []
	s2 = open(hamstr_map_oma).read().split('\n')
	for i in range(len(s2) - 1):
		species_in_tree.append(s2[i].split('\t')[-1])

	for i in range(0, len(s1) - 2, 2):
		if s1[i][0] == '>'and s1[i][1:6] in species_in_tree:
			fnew.write(s1[i] + '\n' + s1[i+1] + '\n')
	fnew.close()

# Remove the HaMStR created directories
def removeHaMStRdirs():
	if os.path.exists(hamstr + '/core_orthologs/' + protein_id):
		os.system('rm -Rf %s/core_orthologs/%s' %(hamstr, protein_id))
	if os.path.exists(hamstr + '/output/' + protein_id):
		os.system('rm -Rf %s/output/%s' %(hamstr, protein_id))

# The main method to run HaMStR search
def main(hamstrFile, coreOrtholog, protId, hamstrMapOma, formatdb, blast, del_temp):
	print '##### Running HaMStR search #####'
	
	global hamstr, core_ortholog, protein_id, hamstr_map_oma, core_ortholog_prot_dir, output, fasta_file, aln_file, hmm_file, format_db, blastp
	global core_ortholog_prot_dir, output, delTemp
	protein_id = protId
	core_ortholog = coreOrtholog
	hamstr = hamstrFile
	core_ortholog_prot_dir = hamstr + '/core_orthologs/' + protein_id
	output = hamstr + '/output/' + protein_id
	hamstr_map_oma = hamstrMapOma
	format_db = formatdb
	blastp = blast
	delTemp = del_temp
		
	print 'HaMStR processing Step 1: Removing old directories'
	remove_old_dir()
	print 'HaMStR processing Step 2: Making new directories'
	make_new_dir()
	print 'HaMStR processing Step 3: Edit core orthologs file'
	success = copy_edit_core_ortholog()
	print 'HaMStR processing Step 4: Creating alignment and HMMER file'
	if success:
		align_hmm_run()
		print 'HaMStR processing Step 5: Running HAMSTR'
		try:
			hamstr_run()
		except KeyboardInterrupt:
			sys.exit('Keyboard interruption by user!!!')
		print 'HaMStR processing Step 6: Adding orthologs to the core set'
		write_output()
		print 'HaMStR processing Step 7: Editing the orthologs and keeping only those present in the species tree'
		finalEditOrthologsSeqs()
		print 'HaMStR processing Step 8: Removing the directories created by HaMStR'
	if delTemp:
		removeHaMStRdirs()

	return success
			
				

	
		

