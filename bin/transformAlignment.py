import os, sys

# Script to transform MSA into indel blocks
# The output file is used as input for IQTree24 to calculate indel rates

def preprocessPhy(phyFile):
	f = open(phyFile).read().split('\n')
	totalSpecies = int(f[0].split()[0])
	fedit = []
	c = 0
	for line in f[1:totalSpecies + 1]:
		complete_sequence = line[11:].replace(' ','')
		addLine = c + totalSpecies + 2
		while addLine < len(f):
			complete_sequence += f[addLine][11:].replace(' ', '')
			addLine += totalSpecies + 1
		fedit.append(complete_sequence)
		c += 1
	#print fedit
	#print len(fedit)

	'''### TEST ###
	fnew = open('test.txt', 'w')
	for line in fedit:
		fnew.write(line + '\n')
	fnew.close()
	###'''

	return fedit

# Calculates the number of indel events occuring in the alignment
# Divides the complete alignment based on indel blocks and stores the column number where indel occurred
def calculateIndelBlocks(f):
	indelBlocks = []

	indelRows = []
	for i in range(len(f[0])):
		indelStart = False
		for j in range(len(f)):
			if f[j][i] == '-':
				if j not in indelRows:
					indelRows.append(j)
					indelStart = True
			else:
				if j in indelRows:
					indelRows.remove(j)
		#print indelRows
		if indelStart:
			indelBlocks.append(i)
	
	#print indelBlocks
	return indelBlocks

# Creates the transformed MSA file
def createTransformAlign(transFile, phyFile, editPhyFile, indelsPos):
	fnew = open(transFile, 'w')
	phy = open(phyFile).read().split('\n')
	totalSpecies = int(phy[0].split()[0])
	fnew.write(str(totalSpecies) + ' ' + str(len(indelsPos)) + '\n')
	c = 0
	for line in phy[1:totalSpecies + 1]:
		fnew.write(line.split()[0])
		if len(indelsPos) == 1:
			startPos = int(indelsPos[0])
			fnew.write(' ' + str(editPhyFile[c][startPos:len(editPhyFile[c])].count('-')))
		elif len(indelsPos) > 1:
			for i in range(len(indelsPos) - 1):
				startPos = int(indelsPos[i])
				stopPos = int(indelsPos[i + 1])
				#print editPhyFile[c][startPos, stopPos]
				fnew.write(' ' + str(editPhyFile[c][startPos:stopPos].count('-')))
			fnew.write(' ' + str(editPhyFile[c][stopPos:len(editPhyFile[c])].count('-')))
		else:
			fnew.write(' ' + str(editPhyFile[c].count('-')))
		fnew.write('\n')
		c += 1
	fnew.close()
		
def main(phy_file, trans_file):		
	if not os.stat(phy_file).st_size == 0:
		edit_file = preprocessPhy(phy_file)
		indelBlocksPos = calculateIndelBlocks(edit_file)
		createTransformAlign(trans_file, phy_file, edit_file, indelBlocksPos)
		return len(indelBlocksPos)
	else:
		print 'Transformed alignment cannot be created for empty phylip files!'	
	
	 
