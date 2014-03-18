'''
Created on Apr 24, 2013

@author: lunt
'''

from itertools import chain as iterchain

import Bio.PDB as P
import scipy as S
import scipy.spatial.distance as SPD

def dists_from_coords(coords):
	"""
	Each row of coords should be the coordinate vector for one object.
	
	@param coords: An Nxd array of N p-dimensional coordinates.
	
	@return: An NxN matrix with at least the upper triangular section giving the distances between elements i and j.
	"""
	#The original version Used the law of cosines, M{c^2 = a^2 + b^2 - 2(a)(b)cos(C)} to find the distances. (a,b,c are edge lengths, and C is the angle opposite c)
	#This is faster
	return SPD.squareform(SPD.pdist(coords))

def all_atom_min_distances(reslist):
	"""
	Foobar...
	"""
	N = len(reslist)
	other_output_dict = dict()
	closest_atom = dict()
	other_output_dict['closest_atom'] = closest_atom
	
	#each residue might have a variable number of atoms, just to make things hard.
	#It might not even be standard based on which amino acid it is, for example with unresolved atoms, having or hiding hydrogen, or wahtever. >.<
	num_atoms_per_res = [len(oneres) for oneres in reslist]
	atom_offset = S.cumsum([0] + num_atoms_per_res) #beginning indexes for each residue, an extra end element gives the total length to make it easy to use, see below
	
	#get one flattened list of atoms
	allatoms = list(iterchain(*[[atom for atom in res] for res in reslist]))#We can trust this more than P.Selection.unfold... This keeps atoms in the same order.
	
	#create a single flat list of atom coordinates
	allcoords = S.array([atom.coord for atom in allatoms]) #and turn it into an Nx3 matrix, N is the number of residues
	
	atom_D = dists_from_coords(allcoords)#use the law of cosines to quickly calculate all the distances.
	
	final_D = S.zeros((N,N))
	
	for i in range(N):
		for j in range(i+1,N):
			submatrix = atom_D[atom_offset[i]:atom_offset[i+1],atom_offset[j]:atom_offset[j+1]]
			minspot = submatrix.argmin()
			min_dist = submatrix.flat[minspot]
			
			final_D[i,j] = min_dist
			final_D[j,i] = min_dist
			
			i_atom = allatoms[ atom_offset[i] + minspot / submatrix.shape[1] ]
			j_atom = allatoms[ atom_offset[j] + minspot % submatrix.shape[1] ]
			
			closest_atom[i,j] = (i_atom,j_atom)
			closest_atom[j,i] = (j_atom,i_atom)
	
	#end nested for

	return final_D, other_output_dict

class DistanceMap(object):
	def __init__(self,some_entity):
		self.residue_list = P.Selection.unfold_entities(some_entity,'R')
		self.D, self.otheroutput = all_atom_min_distances(self.residue_list)
		self.closest_pairs = self.otheroutput['closest_atom']

	def write(self,filehandle,formatstring='%(iid)s\t%(jid)s\t%(dist).4f\n',header='',footer=''):
		"""
		The keys that can be used in formatstring are :
			i	: Iteration parameter.
			j	: Iteration parameter, greater than i
			dist : Distance between residues.
			iid : The ID of the ith residue.
			jid : The ID of the jth residue.
			ichain : Which chain is the i residue on?
			jchain :
			iatom : Which atom was the closest from the ith residue?
			jatom : Which atom was the closest from the jth residue?
		
		"""
		M = len(self.residue_list)
		
		allIDs = [(str(i.id[1]) + i.id[2]).strip() for i in self.residue_list]
		allChains = [i.parent.id for i in self.residue_list]
		
		filehandle.write(header)
		
		for i in range(M):
			for j in range(i+1,M):
				left_atom, right_atom = self.closest_pairs[i,j]
				thedict = {'i':i,'j':j,'dist':self.D[i,j],'iid':allIDs[i],'jid':allIDs[j],
						'ichain':allChains[i],'jchain':allChains[j],'iatom':left_atom.id,'jatom':right_atom.id}
				filehandle.write(formatstring % thedict)
		
		filehandle.write(footer)
		#return
	
	#end of DistanceMap