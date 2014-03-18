import os
#from backmapper.readbackmap import *

#def backMapFileStruct(fileList):
#	fileList.sort()
#	fileNamesByStrand = dict()
#	for i in fileList:
#		strand = i.split(".")[0]
#		domListForThisStrand = fileNamesByStrand.setdefault(strand,list())
#		domListForThisStrand.append(i)
	
#	for i in fileNamesByStrand.values():
#		i.sort()

#def backMapStruct(bundDir, fileStruct):
#	retDict = dict()
#	for k in fileStruct:
#		strandDict = dict()
#		retDict[k] = strandDict
#		for aMapFile in fileStruct[k]:
#			strandDict[aMapFile] = readBackmapList(os.path.join(bundDir,aMapFile))
#	return retDict

#class MapBundle:
#	def __init__(self,bundDir):
#		fileList = os.listdir(bundDir)
#		
#		ents = [i for i in fileList if i.endswith("ent")]
#		maps = [i for i in fileList if i.endswith("map")]
#		faas = [i for i in fileList if i.endswith("faa")]
#		dssps = [i for i in fileList if i.endswith("dssp")]
#		#the lists should only have one element.
#		self.pdbFile = ents[0]
#		self.mapFiles = backMapFileStruct(maps)
#		self.faaFile = faas[0]
#		self.dsspFile = dssps[0]

#		self.backMaps = backMapStruct(bundDir,self.mapFiles)

def create_bundle(pdbid,suffix='.mapbundle',parentdir=''):
	bundleDir = os.path.join(parentdir,pdbid.lower() + suffix)
	try:
		os.makedirs(bundleDir)
	except:
		pass
	try:
		os.makedirs(os.path.join(bundleDir,'Chains'))
	except:
		pass
	try:
		os.makedirs(os.path.join(bundleDir,'Domains'))
	except:
		pass
	
	return bundleDir
	