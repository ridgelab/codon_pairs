#! /usr/bin/env python

import sys
import argparse
from multiprocessing import Process, current_process, freeze_support, Pool
import re
import os


def makeAllPossibleCodons(rna):
	'''
	Input is a flag to specify if the sequence is DNA or RNA.
	Returns a set of all 64 possible codons.
	'''

	from itertools import product
	codons = product("ACGT",repeat=3)
	if rna:
		codons = product("ACGU",repeat=3)
	codonsComb = set()
	for c in codons:
		codonsComb.add("".join(c))
	return codonsComb
def parseArgs():
	'''
	Argument parsing is done.
	Required to have an input file.
	'''
	parser = argparse.ArgumentParser(description='Find Identical and co-tRNA codon pairing.')
	parser.add_argument("-t",help="Number of Cores",action="store",dest="threads",default=0,type=int, required=False)
	parser.add_argument("-i",help="Input Fasta Files",nargs='*',action="store", dest="input", required=False)
	parser.add_argument("-id",help="Input Directory with Fasta Files",action="store", dest="inputDir", required=False)
	parser.add_argument("-od",help="Output Directory",action="store",dest="output", required=False)
	parser.add_argument("-f",help="Ribosome Footprint",action="store",dest="footprint", type=int, default=9, required=False)
	parser.add_argument("-c",help="Co-tRNA codon pairing",action="store_true",dest="co_trna", required=False)
	parser.add_argument("-rna",help="Flag for RNA sequences",action="store_true",dest="rna", required=False)
	args = parser.parse_args()

	if not args.input and not args.inputDir:
		print "You must supply an input file with either -i or -id"
		sys.exit()
	return args


def getPairs(seq):
	footprint = args.footprint
	pairs = dict()
	codons = []
	if args.co_trna:
		from Bio.Seq import Seq
		from Bio.Alphabet import generic_dna
		from Bio.Alphabet import generic_rna
		if args.rna:
			rna = Seq(seq,generic_rna)
			aa = str(rna.translate())
			codons = re.findall(".",aa)
		else:
			seq = Seq(seq,generic_dna)
			aa = str(seq.translate())
			codons = re.findall(".",aa)
	else:
		codons = re.findall("...",seq)
	lastFound = dict() #key= codon, value= position of last found codon with pairing
	for x in xrange(len(codons)):
		curCodon = codons[x]
		if not curCodon in lastFound or (x - lastFound[curCodon] >= footprint):
			if not curCodon in pairs:
				pairs[curCodon] = [0.0,0.0] #Total times codon is found, total codon pairings
			pairs[curCodon][0] +=1
			lastFound[curCodon] =x
			continue
		pairs[curCodon][1] += 1.0 
		pairs[curCodon][0] +=1
		##pairs[curCodon][1] += (1.0 / len(codons)) #normalizes for length of sequence
		#pairs[curCodon][1] += 1.0 
		lastFound[curCodon] = x

	return pairs


def readOneFile(inputFile):
	'''
	Reads one input file that is supplied as a parameter.
	'''
	input = ""
	header = ""
	sequence = ""
	codonPairs = []
	try:
		if inputFile[-3:] =='.gz':
			import gzip
			input = gzip.open(inputFile,'r')
		else:
			input = open(inputFile,'r')
		for line in input:
			if line[0] =='>':
				if sequence !="":
					codonPairs.append((header,getPairs(sequence)))
				header = line
				sequence = ""
				continue
			sequence +=line.upper().strip()
		if sequence != "":
			codonPairs.append((header,getPairs(sequence)))
			sequence = ""

	except Exception: #If the input file is malformatted, do not stop the program.
		input.close()
		print 'Could not open file ' + str(inputFile)
		return {}
	input.close()
	outputFile = sys.stdout
	if args.output:
		outputFile = open(args.output + inputFile.split("/")[-1].split(".gz")[0],'w')
	total = dict()
	numGenes = 0
	for header,pairs in codonPairs:
		numGenes +=1
		outputFile.write(header)
		if len(pairs) == 0:
			outputFile.write("NONE\n")
		for key,value in pairs.items():
			if not key in total:
				total[key] = [0,0]
			total[key][0] +=value[0]
			total[key][1] +=value[1]
			outputFile.write(key + "\t" +str(value[0]) + "\t" + str(value[1]) + "\n")
	outputFile.write(">Pairing\tTotal_Instances\tTotal_Pairs\n")
	for key,value in total.items():
		outputFile.write(key + "\t" +str(value[0]) + "\t" + str(value[1]) + "\n")
	outputFile.write(">Pairing_normalized_for_number_of_genes\tTotal_Instances\tTotal_Pairs\n")
	normalized = dict()
	for key,value in total.items():
		outputFile.write(key + "\t" +str(float(value[0]) / numGenes) + "\t" +str(float(value[1]) / numGenes) + "\n")
		#normalized[key] = [(float(value[0]) / numGenes),((float(value[1])/float(value[0]))/numGenes)]
		normalized[key] = [(float(value[0]) / numGenes),(float(value[1])/float(value[0]))]
	return normalized
	


def readInputFiles(args):
	'''
	Requires the system arguments to be passed to the function.
	'''
	threads = args.threads
	if threads ==0:
		pool = Pool()
	else:
		pool = Pool(threads)	
	allInputFiles = []
	allSets = set()
	fileToSet = {}
	if args.input:
		allInputFiles = args.input
	elif args.inputDir:
		allFasta = []
		path = args.inputDir
		allFasta = os.listdir(path)
		if path[-1] != '/':
			path += '/'
		allInputFiles = [path +i for i in allFasta]
	if len(allInputFiles) < 1:
		print "At least one input file is required"
		sys.exit()
	n = pool.map(readOneFile,allInputFiles)
	normal = {}
	for normalized in n:
		for key,value in normalized.items():
			if not key in normal:
				normal[key] = [0,0]
			normal[key][0] += normalized[key][0]
			normal[key][1] += normalized[key][1]
	output = open("everythingCombined",'w')
	for key,value in normal.items():
		if any([False if char in ["A","T","C","G"] else True for char in key]):
			continue	
		output.write(key + "\t" + str(value[0]/len(n)) + "\t" + str(value[1]/len(n)) + "\n")
	output.close()

if __name__ =='__main__':
	'''
	Main.
	'''

	args = parseArgs()
	if args.output:
		if args.output[-1] != '/':
			args.output = args.output + '/'
		if not os.path.isdir(args.output):
			os.makedirs(args.output)
	codonsComb = makeAllPossibleCodons(args.rna)
	readInputFiles(args)







