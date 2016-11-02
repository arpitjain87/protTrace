import os, sys
import math

### Script to colourize the species according to the traceability values
### Reads in the generated nexus file by PTP and edits it according to our need
### Represents the traceabilities as different colours on the leaves / tip labels
### INPUT: 1. Nexus tree file
###	   2. HaMStR mapping file
###	   3. Species Id
###	   4. Species tree
###	   5. Decay result file

def calculateMaxLikDist(species1, species2):
	#print species1, species2
	speciesMaxFile = open(sp_max_file).read().split('\n')
	for j in range(len(hamstrMapFile) - 1):
		if species1 == hamstrMapFile[j].split('\t')[1]:
			#hamstr1 = hamstrMapFile[j].split('\t')[0]
			oma1 = hamstrMapFile[j].split('\t')[3]
			break
	for j in range(len(hamstrMapFile) - 1):
		if species2 == hamstrMapFile[j].split('\t')[1]:
			#hamstr2 = hamstrMapFile[j].split('\t')[0]
			#kingdom = hamstrMapFile[j].split('\t')[-1]
			oma2 = hamstrMapFile[j].split('\t')[3]
			break
	#print hamstr1, hamstr2
	#print oma1, oma2
	flag1 = True #To check if the species are present in the species ML dist matrix file. Otherwise, parse likelihood from cache directory.
	flag2 = True

	#*** CHECK HERE FOR THE SEPARATION IDENTIFIER IN SPECIES LIKELIHOOD MATRIX FILE ('\t' or '|') ***#
	sep = '\t'

	for l in range(len(speciesMaxFile) - 1):
		
		if speciesMaxFile[l].split(sep)[0] == oma1:
			rowIndex = l
			flag1 = False
		elif speciesMaxFile[l].split(sep)[0] == oma2:
			columnIndex = l
			flag2 = False

	if flag1 or flag2:
		# Checking for the likelihood score in cache directory
		if os.path.exists(cacheDir + '/' + oma1 + '_' + oma2 + '.lik'):
			return float(open(cacheDir + '/' + oma1 + '_' + oma2 + '.lik').read().split('\n')[0])
		elif os.path.exists(cacheDir + '/' + oma2 + '_' + oma1 + '.lik'):
			return float(open(cacheDir + '/' + oma2 + '_' + oma1 + '.lik').read().split('\n')[0])
		else:
			print 'No likelihood distance found between species: %s and %s. Using default likelihood distance of 1.0!' %(species1, species2)
			return 1.00		
	else:
		return float(speciesMaxFile[rowIndex].split(sep)[columnIndex])

def colourize(speciesName, nexusTreeFile):
	tree = open(nexusTreeFile).read().split('\n')
	for i in range(len(tree) - 1):
		if tree[i] == "\ttaxlabels":
			startTaxa = i + 1
		if tree[i] == ";":
			stopTaxa = i
		if tree[i] == "begin trees;":
			treeLine = i + 2
			break
	fnew = open(nexusTreeFile.replace('.nexus', '_edit.nexus'), 'w')
	for i in range(startTaxa):
		fnew.write(tree[i] + '\n')
	for i in range(startTaxa, stopTaxa):
		tempSpecies = tree[i].split()[-1]
		#print tempSpecies
		if not tempSpecies == speciesName:
			try:
				'''if tempSpecies == "crypa_4119":
					colourCode, traceValue = getColourCode(speciesName, "CRYPHONECTRIA_PARASITICA")
					traceResults.write(speciesName + '\t' + "CRYPHONECTRIA_PARASITICA" + '\t' + str(traceValue) + '\n')
				else:'''
				colourCode, traceValue = getColourCode(speciesName, tempSpecies)
				traceResults.write(speciesName + '\t' + tempSpecies + '\t' + str(traceValue) +'\n')
				fnew.write(tree[i] + '[&!color=#-' + colourCode + ']' + '\n')
			except:
				sys.exit()
				pass
			
		else:
			fnew.write(tree[i] + '\n')
	for i in range(stopTaxa, treeLine):
		fnew.write(tree[i] + '\n')
	fnew.write(speciesTree.replace('\n', '') + '\n')
	for i in range(treeLine + 1, len(tree)):
		if tree[i] == '	set branchLabels.isShown=true;':
			fnew.write("	set branchLabels.isShown=false;" + '\n')
		else:
			fnew.write(tree[i] + '\n')
	fnew.close()

def getColourCode(spName, tempName):
	mlDist = calculateMaxLikDist(spName, tempName)
	#print mlDist
	if decayRate < 0.01:
		traceability = 1
	else:
		traceability = 1 - ((decayPop * math.exp(decayRate * mlDist)) / (1 + decayPop * (math.exp(decayRate * mlDist) - 1)))
	#print traceability
	if traceability <= 1 and traceability >= 0.9:
		colCode = '16711936'
	elif traceability < 0.9 and traceability >= 0.8:
		colCode = '3604736'
	elif traceability < 0.8 and traceability >= 0.7:
		colCode = '2294016'
	elif traceability < 0.7 and traceability >= 0.6:
		colCode = '983296'
	elif traceability < 0.6 and traceability >= 0.5:
		colCode = '256'
	elif traceability < 0.5:
		colCode = '65536'

	return colCode, traceability

def main(nexusTreeFile, mapFile, protId, spTree, plotFigTree, speciesMaxLikFile, species_id, cache_dir):
	global sp_max_file, hamstrMapFile, speciesId, speciesTree, decayRate, decayPop, traceResults, cacheDir

	speciesId = species_id
	sp_max_file = speciesMaxLikFile
	decayRate = 8.5
	decayPop = 0.04
	cacheDir = cache_dir

	try:
		hamstrMapFile = open(mapFile).read().split('\n')
		speciesTree = open(spTree).read()
		decayRate = float(open('decay_summary_%s.txt_parameter' %protId).read().split('\n')[1])
		decayPop = float(open('decay_summary_%s.txt_parameter' %protId).read().split('\n')[0])

	except IOError:
		print 'ERROR: Colourizing tree encountered problem!!!'

	#traceResults = open('trace_results_%s_onlyTreeSpecies.txt' %protId, 'w')
	#traceResultsAll = open('trace_results_%s_allSpecies.txt' %protId, 'w')
	traceResults = open('trace_results_%s.txt' %protId, 'w')

	if os.path.exists(nexusTreeFile):
	#speciesId = nexusTreeFile.split('_')[1].split('.')[0][:5]
		for i in range(len(hamstrMapFile) - 1):
			if speciesId == hamstrMapFile[i].split('\t')[3]:
				speciesName = hamstrMapFile[i].split('\t')[1]
				break

		colourize(speciesName, nexusTreeFile)

		try:
			os.system('java -cp %s figtreepdf %s' %(plotFigTree, nexusTreeFile.replace('.nexus', '_edit.nexus')))
		except:
			print 'WARNING: No representation of traceabilities on tree possible.\nJAVA program figtreepdf not responding!!!'
	traceResults.close()

	### Calculating traceabilities for all species in the input species maximum likelihood distances file
	'''for i in range(len(hamstrMapFile) - 1):
		if speciesId == hamstrMapFile[i].split('\t')[3]:
			speciesName = hamstrMapFile[i].split('\t')[1]
			break
	for line in open(sp_max_file):
		tempSpecies = line.split()[0]
		for line2 in open(mapFile):
			if tempSpecies in line2:
				tempSpecies = line2.split()[1]
				break
		#print "tempSpecies:", tempSpecies
		try:
			if not tempSpecies == speciesName:
				if tempSpecies == "crypa_4119":
					colourCode, traceValue = getColourCode(speciesName, "CRYPHONECTRIA_PARASITICA")
					traceResultsAll.write(speciesName + '\t' + "CRYPHONECTRIA_PARASITICA" + '\t' + str(traceValue) + '\n')
				#else:
				colourCode, traceValue = getColourCode(speciesName, tempSpecies)
				traceResultsAll.write(speciesName + '\t' + tempSpecies + '\t' + str(traceValue) +'\n')
			else:
				traceValue = 1
				traceResultsAll.write(speciesName + '\t' + tempSpecies + '\t' + str(traceValue) + '\n')
		except:
			pass
	traceResultsAll.close()'''
	
