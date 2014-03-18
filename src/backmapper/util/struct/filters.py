'''
Created on Apr 24, 2013

@author: lunt
'''

#TODO: Totally replace this with Bio.PDB.PDBIO.Select based filtering.

from copy import deepcopy

import Bio.PDB as P

def _dict_alternate_haskey(aDict,alternateKeys):
	for aKey in alternateKeys:
		if aDict.has_key(aKey):
			return aKey
	raise KeyError(alternateKeys)

def _dict_alternate(aDict, alternateKeys):#function to specify a list of alternatives for a dictionary get.
	for aKey in alternateKeys:
		if aDict.has_key(aKey):
			return aDict[aKey] #if it is found
	raise KeyError(alternateKeys) #if not, consistent with normal dictionary subscripting

class StructFilter(object):
	def __call__(self,someEntity):
		raise NotImplementedError("Abstract Base Class")

class NoFilter(StructFilter):
	"""
	Does Nothing, useful in some situations where you have code that expects a filter.
	This allows you to remove a lot of oddball case checking.
	"""
	def __call__(self,someEntity):
		return someEntity

class AAOnly(StructFilter):
	def __call__(self,someEntity):
		"""
		Destructive
		"""
		retval = deepcopy(someEntity) #Copy doesn't make a deep copy of P.PDB.Residue.DisorderedResidue, it just gets the currently selected disordered residue
		all_res_list = P.Selection.unfold_entities(retval, 'R')
		
		for oneres in all_res_list:
			if not P.Polypeptide.is_aa(oneres):#Skip non-amino acid residues.
				theParent = oneres.parent
				theParent.detach_child(oneres.id)

		return retval

class AtomFilter(StructFilter):
	
	def __init__(self,prefatoms=['CA'],remove=None,keep=None,view=False):
		"""
		@param prefatoms: A list of strings, giving alternative atoms in order of preference (ex ['CB', 'CA']), the first atom that is present in the residue is taken.
		@type prefatoms: list of str
		
		@param remove: A list of strings for residues that are always removed from each residue
		@type remove: list
		
		@param keep: A list of strings for residues that are always kept from each residue
		@type keep: list
		"""
		self.prefatoms = prefatoms
		self.remove = remove
		self.keep = keep
		self.view = view
		
	def _filter_single_residue(self,oneres):
		"""
		Warning, this is a destructive function.
		
		It should work with DisorderedEntities now?
		"""
		
		assert isinstance(oneres,P.Residue.Residue)
		
		toRemove = set([i.id for i in oneres.child_list]) #The list of atoms we want to remove from this residue.
		try:
			keepone = _dict_alternate_haskey(oneres.child_dict,self.prefatoms) #remove the ones we don't want to remove from the remove list.
			toRemove.remove( keepone )
		except:
			pass
			#TODO: Add informative stderr statement.
			
		
		if self.keep is not None:
			toRemove.difference_update(self.keep)

		
		for oneID in toRemove:
			oneres.detach_child(oneID)
		#return

	def __call__(self,someEntity):
		"""
		Destructive
		"""
		retval = deepcopy(someEntity) #Copy doesn't make a deep copy of P.PDB.Residue.DisorderedResidue, it just gets the currently selected disordered residue
		all_res_list = P.Selection.unfold_entities(retval, 'R')
		
		for oneres in all_res_list:
			if not P.Polypeptide.is_aa(oneres):#Skip non-amino acid residues.
				continue
			if isinstance(oneres,P.Residue.DisorderedResidue):#Disordered residue
				for subres in oneres.disordered_get_list():
					self._filter_single_residue(subres)
			else:#Normal Residue
				self._filter_single_residue(oneres)
			#never makes it here.

		
		return retval

class RemoveEmptyResidues(StructFilter):
	def __call__(self,someEntity):
		retval = deepcopy(someEntity)
		all_res_list = P.Selection.unfold_entities(retval,'R')
		
		for oneres in all_res_list:
			if len(oneres.child_list) == 0:
				oneres.get_parent().detach_child(oneres.id)
		return retval