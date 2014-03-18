'''
Created on Apr 26, 2013

@author: lunt
'''

import Bio.PDB as PDB

import os
from urllib import urlretrieve as download

__seqres_url = "http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=FASTA&compression=NO&structureId=%s"
__dssp_url = "http://www.pdb.org/pdb/files/%s.dssp"
def download_PDB_and_FASTA(modelname,targetdir='.'):
	plist = PDB.PDBList()
	
	try:
		modelfile = plist.retrieve_pdb_file(modelname,pdir=targetdir)#this handles directory creation
	except:
		raise Exception("Failed to Retrieve PDB file")
	
	targetdir = os.path.dirname(modelfile)
	
	url = __seqres_url % modelname.upper()
	fastafile = os.path.join(targetdir,modelname+".faa") # this makes the assumption that modelfile is in fact the filename of the model.
	download(url,fastafile)
	
	#download the DSSP output
	dsspUrl = __dssp_url % modelname.lower()
	dsspSaveFile = os.path.join(targetdir,modelname+".dssp")
	download(dsspUrl,dsspSaveFile)
	
	return modelfile, fastafile