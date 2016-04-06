# Module to read in program config file and prepare variables
# to be imported in other modules

import os, sys

class setParams:

	def __init__(self, configFile):
		config = open(configFile).read().split('\n')
		for line in config:
			if line != "" and line[0] != "#":
				if line.split(':')[0] == 'species':
					self.species = line.split(':')[1]
				if line.split(':')[0] == 'search_proteome':
					if line.split(':')[1] == 'YES':
						self.search_proteome = True
					else:
						self.search_proteome = False
				if line.split(':')[0] == 'search_ortholog_groups':
					if line.split(':')[1] == 'YES':
						self.search_ortholog_groups = True
					else:
						self.search_ortholog_groups = False
				if line.split(':')[0] == 'search_ortholog_sequences':
					if line.split(':')[1] == 'YES':
						self.search_ortholog_sequences = True
					else:
						self.search_ortholog_sequences = False
				if line.split(':')[0] == 'run_hamstr':
					if line.split(':')[1] == 'YES':
						self.run_hamstr = True
					else:
						self.run_hamstr = False
				if line.split(':')[0] == 'run_hamstrOneSeq':
					if line.split(':')[1] == 'YES':
						self.run_hamstrOneSeq = True
					else:
						self.run_hamstrOneSeq = False
				if line.split(':')[0] == 'preprocessing':
					if line.split(':')[1] == 'YES':
						self.preprocessing = True
					else:
						self.preprocessing = False
				if line.split(':')[0] == 'traceability_calculation':
					if line.split(':')[1] == 'YES':
						self.traceability_calculation = True
					else:
						self.traceability_calculation = False
				if line.split(':')[0] == 'calculate_scaling_factor':
					if line.split(':')[1] == 'YES':
						self.calculate_scaling_factor = True
					else:
						self.calculate_scaling_factor = False
				if line.split(':')[0] == 'calculate_indel':
					if line.split(':')[1] == 'YES':
						self.calculate_indel = True
					else:
						self.calculate_indel = False
				if line.split(':')[0] == 'perform_msa':
					if line.split(':')[1] == 'YES':
						self.perform_msa = True
					else:
						self.perform_msa = False
				if line.split(':')[0] == 'delete_temporary_files':
					if line.split(':')[1] == 'YES':
						self.delete_temp = True
					else:
						self.delete_temp = False
				if line.split(':')[0] == 'map_traceability_tree':
					if line.split(':')[1] == 'YES':
						self.mapTraceabilitySpeciesTree = True
					else:
						self.mapTraceabilitySpeciesTree = False
				if line.split(':')[0] == 'aa_substitution_matrix':
					self.aa_substitution_matrix = line.split(':')[1]
				if line.split(':')[0] == 'simulation_runs':
					self.simulation_runs = int(line.split(':')[1])
				if line.split(':')[0] == 'msa':
					self.msa = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'tree_reconstruction':
					self.tree_reconstruction = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'REvolver':
					self.REvolver = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'hmmfetch':
					self.hmmfetch = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'hmmscan':
					self.hmmscan = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'iqtree24':
					self.iqtree24 = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'clustalw':
					self.clustalw = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'blastall':
					self.blastall = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'blastp':
					self.blastp = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'formatdb':
					self.formatdb = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'makeblastdb':
					self.makeblastdb = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'R':
					self.R = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'hamstr':
					self.hamstr = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'hamstrOneSeq':
					self.hamstrOneSeq = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'degapping':
					self.degapping = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'treePuzzle':
					self.treePuzzle = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'parameters_treePuzzle':
					self.parameters_treePuzzle = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'species_MaxLikMatrix':
					self.species_MaxLikMatrix = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'path_work_dir':
					self.path_work_dir = os.path.abspath(line.split(':')[1])
					if not os.path.exists(self.path_work_dir):
						os.mkdir(self.path_work_dir)
				if line.split(':')[0] == 'hamstr_oma_tree_map':
					self.hamstr_oma_tree_map = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'path_oma_seqs':
					self.path_oma_seqs = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'path_oma_group':
					self.path_oma_group= os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'pfam_database':
					self.pfam_database = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'species_tree':
					self.species_tree = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'species_tree_msa':
					self.species_tree_msa = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'simulation_tree':
					self.simulation_tree = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'decay_script':
					self.decay_script = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'plot_figtree':
					self.plot_figtree = os.path.abspath(line.split(':')[1])
				if line.split(':')[0] == 'proteome':
					self.proteome = line.split(':')[1]
				if line.split(':')[0] == 'oma_ortholog_ids':
					self.oma_ortholog_ids = line.split(':')[1]
				if line.split(':')[0] == 'oma_ortholog_sequences':
					self.oma_ortholog_sequences = line.split(':')[1]
				if line.split(':')[0] == 'hamstr_ortholog_sequences':
					self.hamstr_ortholog_sequences = line.split(':')[1]
				if line.split(':')[0] == 'ortholog_group_alignment':
					self.ortholog_group_alignment = line.split(':')[1]
				if line.split(':')[0] == 'transformed_alignment':
					self.transformed_alignment = line.split(':')[1]
				if line.split(':')[0] == 'treeFile':
					self.treeFile = line.split(':')[1]
				if line.split(':')[0] == 'REvolver_config_file':
					self.REvolver_config_file = line.split(':')[1]
				if line.split(':')[0] == 'blast_result_file':
					self.blast_result_file = line.split(':')[1]
				if line.split(':')[0] == 'decay_file':
					self.decay_file = line.split(':')[1]
				if line.split(':')[0] == 'scaling_factor_file':
					self.scaling_factor_file = line.split(':')[1]
				if line.split(':')[0] == 'nexus_tree_file':
					self.nexus_tree_file = line.split(':')[1]
	

	
