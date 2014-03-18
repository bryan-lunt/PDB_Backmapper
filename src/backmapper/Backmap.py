'''
Created on May 26, 2011

@author: lunt
'''

import os
import re

class Backmap(list):
	'''
	This Object stores backmaps, hopefully in a properly Pythonic way.
	'''


	def __init__(self,pdbid="",model="",chain="",name="",seqrec=None,struct=None):
		'''
		Constructor
		'''
		self.pdbid = pdbid
		self.model = model
		self.chain = chain
		self.name = name
		self.seqrec = seqrec
		self.struct = struct

	def add_res(self,hmmCol,hitChar,oneRes):
		self.append((hmmCol,hitChar,oneRes))
	
	
	def __add__(self,other):
		retval = Backmap(",".join([self.pdbid,other.pdbid]), ",".join([self.model,other.model]), ",".join([self.chain,other.chain]), ",".join([self.name,other.name]),)
		retval.extend(self)
		startint = self[-1][0]#There are no inserts afterward, so this is always a hit state.
		foo = other.copy()
		foo.column_offset(startint + 1)
		retval.extend(foo)

	def column_offset(self,anint):
		for i in range(len(self)):
			foo = self[i]
			if foo[0] is not None:
				self[i] = (foo[0]+anint,foo[1],foo[2])
	
	def num_columns(self):
		return len([i for i in self if i[0] is not None])
	
	def remove_inserts(self):
		for i in range(len(self)-1,-1,-1):
			if self[i][0] is None:
				self.pop(i)
	
	def remove_gaps(self):
		for i in range(len(self)-1,-1,-1):
			if self[i][2] is "-":
				self.pop(i)

	def write(self,outfile):
		outfile.write("#HMM col\tPDB res ID\tHmmHit Residue\n")
		for hmmCol,hitCharacter,one_residue in self.__iter__():
			if hmmCol is None:
				hmmCol = ''
			
			formattedID=''
			if one_residue is not None:
				formattedID = format_resid(one_residue)
			
			outfile.write("%s\t%s\t%s\n" % (hmmCol,formattedID,hitCharacter))
		#end of loop
		#return

def format_resid(residue):
	hetid, idnum, insertid = residue.id
	idnum = str(idnum)
	insertid = insertid.strip()
	return "%s%s"%(idnum,insertid)

class FakeResidue(object):
	def __init__(self,het,idnum,insert):
		self.id = (het,idnum,insert)

def read_backmap(filehandle):
	
	strand, number, name = os.path.basename(filehandle.name).split(".")[:3]#don't need the .map
	
	retmap = Backmap(chain=strand,name=name)
	
	for line in filehandle:
		if line.startswith("#"):
			continue
		hmmCol,pdbResID,resLetter = line.split("\t")[:3]
		try:
			hmmColint = int(hmmCol)
		except:
			hmmColint = None
			
		pdbFakeRes = None
		try:
			pdbResIDmatch = re.match('(^\d+)([^\d]?$)',pdbResID)
			pdbResIDint = int(pdbResIDmatch.group(1))
			pdbResIDinsert = pdbResIDmatch.group(2)
			pdbFakeRes = FakeResidue(' ',pdbResIDint,pdbResIDinsert)
		except:
			pass
		retmap.add_res(hmmColint,resLetter.strip(),pdbFakeRes)
		
	return retmap
