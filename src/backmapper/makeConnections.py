#!/usr/bin/python
"""A program for creating tab-separated-values files representing predicted contact residues on a PDB model.
"""
import sys, os, os.path
from re import split
import csv
from optparse import OptionParser
from itertools import count, izip, repeat

import Backmap
import Bio.PDB as P

import glob
def chooseMapFiles(bundleFile, nameList):
	mapFiles = list()
	for oneName in nameList:
		pattern = os.path.join(bundleFile,"*"+oneName+"*")
		candidates = glob.glob(pattern)
		if len(candidates) == 0:
			raise Exception('No mapfile matching "' + pattern + '".')
		elif len(candidates) == 1:
			mapFiles.append(candidates[0])
		else:
			raise Exception('Pattern "' + pattern +'" is ambiguous.')
	
	return mapFiles
	
def readDistFile(filelike):
	"""
	Expects a distfile to have this header:
	"#resID1	chain1	resID2	chain2	atom1	atom2	distance\n"
	
	"""
	retval = dict()
	for line in filelike:
		if line.startswith("#"):
			continue
		resID1,chain1,resID2,chain2,atom1,atom2,distance = line.strip().split()
		retval[(resID1, chain1, resID2, chain2)] = ( atom1, atom2, float(distance))
	return retval



def main():
	parser = OptionParser()
	parser.set_usage(parser.get_prog_name() + " [options] <mapbundle> <MIFILE or '-' for STDIN> <FirstMap> [<list of subsequent maps>] \n")
	options, args = parser.parse_args()
	if(len(args) < 3):
		parser.print_help()
		sys.exit()

	bundleFile = args[0]
	distfile = os.path.join(bundleFile,"distances")
	dsspFile = os.path.join(bundleFile,bundleFile.split(".")[0] + ".dssp")
	MIfile = args[1]
	mapFileLists = args[2:]
	mapFiles=list()
	for fileOrFileList in mapFileLists:
		mapFiles.extend(split(",| ",fileOrFileList))
	mapFiles = chooseMapFiles(bundleFile,mapFiles)

	#read top couplets
	with (sys.stdin if MIfile == "-" else open(MIfile)) as reader:
		Nlines = [l for l in reader]
		topcups = [[int(j) for j in i.split()[:2]] for i in Nlines if not i.startswith("#")]

	#read Distance file if it exists and we can open it.
	resDist = dict()
	try:
		with open(distfile) as distfile_file:
			resDist = readDistFile(distfile_file)
	except Exception as e:
		sys.stderr.write("Could not open distance file: " + e.message)
		sys.exit(1)

	#create a mapping of HMMColumnID -> Strand
	

	#read all specified backmaps
	#if options.strands is None:
	
	allBackmaps = [Backmap.read_backmap(open(i)) for i in mapFiles]
	for i in allBackmaps:
		i.remove_inserts()
	
	strandMap = reduce(lambda x,y:x+y , [[i for i in repeat(j.chain,len(j))] for j in allBackmaps])																#strandMap[colID] gives the strand for that column
	
	completeMap = reduce(lambda x,y:x+y, allBackmaps)
	
	#read DSSP dictionary
	SASdict = dict()
	if dsspFile is not None and os.path.exists(dsspFile):
		dsspDict = P.make_dssp_dict(dsspFile)
		SASdict = dict([((a[0],a[1][1]), b[2]) for a, b in dsspDict[0].iteritems()])

	print """#hmmI	hmmJ	pdbInum	pdbIstrd	pdbIatom	pdbJnum	pdbJstrd	pdbJatom	iSAS	jSAS	dist"""
	for i,n in izip(topcups,count(1)):
		
		hmmI = i[0]
		hmmJ = i[1]
		
		ls = strandMap[hmmI]
		leftRes = completeMap[hmmI][2]
		if leftRes is None:
			leftRes = "?"
		else:
			leftRes = Backmap.format_resid(leftRes)
		
		rs = strandMap[hmmJ]
		rightRes = completeMap[hmmJ][2]
		if rightRes is None:
			rightRes = "?"
		else:
			rightRes = Backmap.format_resid(rightRes)
		
		leftSAS = SASdict.get((ls,leftRes),-1)
		rightSAS = SASdict.get((rs,rightRes),-1)
		atom1,atom2,distance = resDist.get((leftRes,ls,rightRes,rs),("C","C",-1.0))		#default to C-Alpha if no distance included
		
		print """%i	%i	%s	%s	%s	%s	%s	%s	%.2f	%.2f	%f""" % (hmmI, hmmJ, leftRes, ls, atom1, rightRes, rs, atom2, leftSAS, rightSAS, distance)

if __name__ == "__main__":
	main()
