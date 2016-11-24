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
# Edits the fasta file how its required for the HAMSTR run
def copy_edit_core_ortholog():
	global fasta_file
	fasta_file = core_ortholog_prot_dir + '/' + protein_id + '.fa'
	s1 = open(core_ortholog).read().split('\n')
	fnew = open(fasta_file, 'w')
	flag = True #Check whether HaMStR has a starting sequence#
	for i in range(0, len(s1)-1, 2):
		if s1[i][0] == '>':
			s2 = open(hamstr_map_oma).read().split('\n')
			for j in range(len(s2)-1):
				if s2[j].split()[3] == s1[i][1:6]:
					hamstrProtId = s2[j].split('\t')[0]
					blastDir = hamstr + '/blast_dir'
					for dirs in glob.glob(blastDir + '/*'):
						if hamstrProtId in dirs:
							blastFastaFile = dirs + '/' + dirs.split('/')[-1] + '.fa'
							break
					s3 = open(blastFastaFile).read().split('\n')
					for k in range(len(s3) - 1):
						if s1[i+1] == s3[k].replace('*', ''):
							fnew.write(s1[i][0] + protein_id + '|' + dirs.split('/')[-1] + '|' + s3[k-1][1:] + '\n' + s1[i+1] + '\n')
							flag = False
							break
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

	for dirs in glob.glob(hamstr + '/genome_dir/*'):
		#print dirs
		if not ".fa" in dirs and not ".sql" in dirs:
			hamstr_name = dirs.split('/')[-1].split('_')[0] + '_' + dirs.split('/')[-1].split('_')[1].split('@')[0]
			if hamstr_name in species_in_tree:
				file_name = dirs + '/' + dirs.split('/')[-1] + '.fa'
				command = '%s/bin/hamstr.pl -central -sequence_file=%s -taxon=misc -hmmset=%s -strict -representative -outpath=%s -hit_limit=10' %(hamstr, file_name, protein_id, output)
				print 'HaMStR run with command: ', command
			#os.system('%s/bin/hamstr -sequence_file=%s -taxon=y_rbf -hmmset=%s -strict -representative -protein -outpath=%s -hit_limit=1' %(hamstr, file_name, protein_id, output))
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

			hamstr_name = files.split('/')[-1].split('_')[1] + '_' + files.split('/')[-1].split('_')[2].split('@')[0]
			s3 = open(hamstr_map_oma).read().split('\n')
			
			oma_name = ''
			for j in range(len(s3)-1):
				if s3[j].split()[0] == hamstr_name:
					#oma_name = s3[j].split()[3] + s2[0].split('|')[3][:5]
					oma_name = s3[j].split()[3]
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
		species_in_tree.append(s2[i].split('\t')[3])

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
def main(hamstrFile, coreOrtholog, protId, hamstrMapOma, delTemp):
	print '##### Running HaMStR search #####'
	
	global hamstr, core_ortholog, protein_id, hamstr_map_oma, core_ortholog_prot_dir, output, fasta_file, aln_file, hmm_file
	global core_ortholog_prot_dir, output
	protein_id = protId
	core_ortholog = coreOrtholog
	hamstr = hamstrFile
	core_ortholog_prot_dir = hamstr + '/core_orthologs/' + protein_id
	output = hamstr + '/output/' + protein_id
	hamstr_map_oma = hamstrMapOma
		
	#print 'HAMSTR 1: Removing old directories'
	remove_old_dir()
	#print 'HAMSTR 2: Making new directories'
	make_new_dir()
	#print 'HAMSTR 3: Edit core orthologs file'
	success = copy_edit_core_ortholog()
	#print 'HAMSTR 4: Creating alignment and HMMER file'
	if success:
		align_hmm_run()
		#print 'HAMSTR 5: Running HAMSTR'
		try:
			hamstr_run()
		except KeyboardInterrupt:
			sys.exit('Keyboard interruption by user!!!')
		#print 'HAMSTR 6: Adding orthologs to the core set'
		write_output()
		#print 'HAMSTR 7: Editing the orthologs and keeping only those present in the species tree'
		finalEditOrthologsSeqs()
		#print 'HAMSTR 8: Removing the directories created by HaMStR'
	if delTemp:
		removeHaMStRdirs()

	return success
			
				

	
		

