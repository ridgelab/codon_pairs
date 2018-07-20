#! /usr/bin/env python
import os
import sys
from multiprocessing import Process, current_process, freeze_support, Pool
import argparse


def parseArgs():
	'''
	Argument parsing is done.
	Required to have an input directory and output file.
	'''
	parser = argparse.ArgumentParser(description='Find Species Matrix for Identical and co-tRNA codon pairing.')
	parser.add_argument("-t",help="Number of Cores",action="store",dest="threads",default=0,type=int, required=False)
	parser.add_argument("-id",help="Input Directory with Fasta Files",action="store", dest="inputDir", required=True)
	parser.add_argument("-o",help="Output Matrix File",action="store",dest="output", required=True)
	parser.add_argument("-oc",help="Output Characters File",action="store",dest="outputChars", required=False)
	parser.add_argument("-on",help="Output Species Names File",action="store",dest="outputNames", required=False)
	parser.add_argument("-c",help="Co-tRNA codon pairing",action="store_true",dest="co_trna", required=False)
	args = parser.parse_args()

	if not args.inputDir:
		print "You must supply an input directory with -id"
		sys.exit()
	return args

def readFile(species):
	'''
	Finds the pairing for each species
	Requires a single input file as a parameter
	Returns a dictionary of the pairing, a list of all codons used, 
	and a dictionary containing the codon and a "0" or "1" for its use.
	'''
	#print 'Reading file ' + species
	dic = {}
	allCodonsSpecies = set()
	codonsToPairingSpecies = {}
	infile = open(args.inputDir +  species, "r")
	currentGene = ""
	for line in infile:
		line = line.strip()
		if line[0] == ">":
			header = line.upper()
			geneName = header.split("GENE=")
			if len(geneName) == 1: #The gene doesn't have a name
				currentGene = ""
				continue
			currentGene = geneName[1].split("] ")[0]
			if currentGene != "":
				length = 64
				if args.co_trna:
					length = 20
				dic[currentGene] = ['?'] * length
		else:
			if currentGene == "":
				continue
			fields = line.split()
			codon = fields[0].upper()
			if codon not in codonList and codon not in aminoList:
				continue
			numPairing = fields[2]
			if float(numPairing) == 0.0:
				numPairing = '0'
			else:
				numPairing = '1'
			geneCodon = currentGene + '__' + codon
			allCodonsSpecies.add(geneCodon)
			codonsToPairingSpecies[geneCodon] = numPairing
			if args.co_trna:
				index = aminoList.index(codon)
			else:
				index = codonList.index(codon)
			dic[currentGene][index] = numPairing
	infile.close()
	return dic,allCodonsSpecies,species,codonsToPairingSpecies


def readAllFiles(fileNames):
	'''
	Processes all the files and recombines the results
	Requires a list of the file names to be passed as a parameter
	Returns a list of informative codons, a dictionary of the pairing for each species,
	and a list of all the species.
	'''
	pairing = {}
	allCodons = set()
	informativeCodons = set()
	if args.threads:
		threads = args.threads
	pool = Pool(threads)
	results = pool.map(readFile,fileNames)
	codonToFreq = {}
	codonsToPairings = {}
	for result in results:
		dic = result[0]
		allCodonsSpecies = result[1]
		species = result[2]
		codonsToPairingSpecies = result[3] 
		for codon in codonsToPairingSpecies:
			if codon in codonsToPairings:
				codonToFreq[codon] += 1
				if codonsToPairings[codon] != codonsToPairingSpecies[codon]:
					codonsToPairings[codon] = "P" #Found both a 0 and a 1, mark the codon as parsimony informative
				if codonsToPairings[codon] == "P" and codonToFreq[codon] ==4:  #Codon has 0 and 1, and also is found in minimum 4 species
					informativeCodons.add(codon)
			else:
				codonsToPairings[codon] = codonsToPairingSpecies[codon]
				codonToFreq[codon] = 1
		pairing[species] = dic
		allCodons = allCodons | allCodonsSpecies
	speciesToNumChars = {}
	speciesList,informativeCodons,pairing = removeCodons(fileNames,informativeCodons,pairing,codonToFreq,speciesToNumChars)
	return informativeCodons,pairing,speciesList

def removeCodons(speciesList,informativeCodons,pairing,codonToFreq,speciesToNumChars):
	'''
	Removes any codon that isn't found in at least 4 species.
	Removes any species that doesn't have at least 5% of the informative codons.
	Recurses until no more changes are made.
	Returns the list of species, list of informative codons, and pairing dictionary
	of species to codons.
	'''
	if len(speciesToNumChars) == 0:
		for species in speciesList:
			num = 0
			i = 0
			j = 0
			for geneCodon in informativeCodons:
				fields = geneCodon.split('__')
				gene = fields[0]
				codon = fields[1]
				if codon not in codonList and codon not in aminoList:
					continue
				if args.co_trna:
					index = aminoList.index(codon)
				else:
					index = codonList.index(codon)
				if gene in pairing[species]:
					if pairing[species][gene][index] != '?':
						num += 1
					else:
						i += 1
				else:
					j += 1
			speciesToNumChars[species] = num

	needsUpdate = False
	codonsToDelete = set()
	for species in speciesList:
		if speciesToNumChars[species] / float(len(informativeCodons)) < .05: 
			speciesList.remove(species)
			for geneCodon in informativeCodons:
				gene = geneCodon.split('__')[0]
				if gene in pairing[species]:
					if geneCodon in codonToFreq:
						codonToFreq[geneCodon] -= 1
						if codonToFreq[geneCodon] < 4:
							codonsToDelete.add(geneCodon)
							needsUpdate = True

	for geneCodon in codonsToDelete:
		informativeCodons.remove(geneCodon)
		fields = geneCodon.split('__')
		gene = fields[0]
		codon = fields[1]
		for species in speciesList:
			if gene in pairing[species]:
				if args.co_trna:
					index = aminoList.index(codon)
				else:
					index = codonList.index(codon)
				if pairing[species][gene][index] != '?':
					speciesToNumChars[species] -= 1

	if needsUpdate:
		removeCodons(speciesList,informativeCodons,pairing,codonToFreq,speciesToNumChars)

	return speciesList,informativeCodons,pairing
	
def readInputFiles(args):
	'''
	Prepares the input files to be read.
	Writes out the results to the files.
	'''
	outfile = open(args.output,"w")
	outfileChars = open(args.outputChars,"w")
	if args.outputNames:
		outfileNames = open(args.outputNames,"w")
	fileNames = []
	for species in os.listdir(args.inputDir):
		fileNames.append(species)
	informativeCodons,pairing,speciesList = readAllFiles(fileNames)
	
	if args.outputChars:
		for codon in informativeCodons:
			outfileChars.write(codon + "\n")
	
	outfile.write("xread\n")
	outfile.write(str(len(informativeCodons)) + " " + str(len(speciesList)) + "\n")
	count = 0
	for species in speciesList:
		if args.outputNames:
			outfile.write("Species_" + str(count) + "\t")
			outfileNames.write("Species_" + str(count) + "\t" + species + "\n")
			count += 1
		else:
			outfile.write(species + "\t")
		for item in informativeCodons:
			fields = item.split('__')
			gene = fields[0]
			codon = fields[1]
			if gene in pairing[species]:
				if args.co_trna:
					index = aminoList.index(codon)
				else:
					index = codonList.index(codon)
				outfile.write(pairing[species][gene][index])
			else:
				outfile.write('?')
		outfile.write('\n')
	outfile.write(";\n")
	outfile.write("proc /;\ncomments 0\n;")
	outfile.close()
	outfileChars.close()

codonList = ['ATT','ATC','ATA','CTT','CTC','CTA','CTG','TTA','TTG','GTT','GTC','GTA','GTG','TTT','TTC','ATG','TGT','TGC','GCT','GCC','GCA','GCG','GGT','GGC','GGA','GGG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','TCT','TCC','TCA','TCG','AGT','AGC','TAT','TAC','TGG','CAA','CAG','AAT','AAC','CAT','CAC','GAA','GAG','GAT','GAC','AAA','AAG','CGT','CGC','CGA','CGG','AGA','AGG','TAA','TAG','TGA']
aminoList = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

if __name__ =='__main__':
	'''
	Main.
	'''
	args = parseArgs()
	if args.inputDir[-1] != '/':
		args.inputDir = args.inputDir + '/'
	readInputFiles(args)
