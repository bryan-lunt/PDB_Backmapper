'''
Created on Jan 21, 2013

@author: lunt
'''

import scipy as S

import Bio.PDB as __P
import Bio.PDB.StructureBuilder as __SB
from Bio.PDB.Polypeptide import PPBuilder as __PPBuilder

try:
	from Bio.PDB.Polypeptide import to_one_letter_code as __PP_one_letter_code
except:
	from Bio.Data.SCOPData import protein_letters_3to1 as __PP_one_letter_code

from Bio.SeqRecord import SeqRecord as __SeqRecord
from Bio.Seq import Seq as __Seq

import Bio.pairwise2 as __PW

#a sortof singleton
__pdb_parser = None
__pp_builder = None

from sys import stderr as __stderr

#we need a workaround because of the assertion in line 69 of Bio.PDB.Atom.__init__
#it asserts that the element name is all upper case... stupid, what about Se!?
class Workaround_Builder(__SB.StructureBuilder):
	def __init__(self):
		super(Workaround_Builder,self).__init__()
	
	def init_atom(self, name, coord, b_factor, occupancy, altloc, fullname, serial_number=None, element=None):
		myElement = None if element is None else element.upper()
		return super(Workaround_Builder,self).init_atom(name, coord, b_factor, occupancy, altloc, fullname, serial_number, element=myElement)
		
		
def get_pdb_parser():
	global __pdb_parser
	if __pdb_parser is None:
		__pdb_parser = __P.PDBParser(get_header=True,structure_builder=Workaround_Builder())
	return __pdb_parser

def get_peptide_builder():
	global __pp_builder
	if __pp_builder is None:
		__pp_builder = __PPBuilder()
	return __pp_builder

def parse_structure(STRUCTFILENAME,name='X'):
	return get_pdb_parser().get_structure(name,STRUCTFILENAME)

def build_polypeptides(one_chain,aa_only=0):	
	return get_peptide_builder().build_peptides(one_chain,aa_only=False)


__default_sequence_header_format = "%(PDBID)s:%(CHAINID)s|PDBID|CHAIN|SEQUENCE%(METHOD)s"

#end basic helpers
#Sequence from headers
def read_all_seqres_from_structfile(structfile,pdbid="",headerfmt=__default_sequence_header_format):
	"""Expects an opened handle to the pdbXXXX.ent file to be read
	Returns a dict of ChainID -> Bio.SeqRecord.SeqRecord objects representing the seqres entries of the pdb file.
	"""
	
	lines = list()
	#get all SEQRES entries
	for line in structfile:
		if line.startswith("SEQRES"):
			lines.append(line)
	
	#TODO: See if you need to worry about MODSEQ
	
	#build sequences for them. I know I could have done it all in the above loop...
	#Keep this structure because we might need to deal with MODSEQ etc in a later version.
	resStrings = dict()#maps chainID -> list of three character 
	for line in lines:
		chainID = line[11]
		ress = line[19:].split()#A list of residues in three letter format
		resStrings.setdefault(chainID,'')
		resStrings[chainID] += ''.join(map(__PP_one_letter_code.get, ress))
	
	#building it all at once will make it use a protein alphabet
	returnDict = dict()
	for chainID, thesequence in resStrings.iteritems():
		sequenceID = headerfmt % {'PDBID':pdbid, 'CHAINID':chainID , 'METHOD':'_SEQRES'}
		returnDict[chainID] = __SeqRecord(__Seq(thesequence),id=sequenceID,description='')
	
	return returnDict



#Sequence from residues
def build_chain_seqrec(one_chain,headerfmt=__default_sequence_header_format):
	#TODO make it get description, organism, etc from the PDB structure header.
	pdbid, chainID = one_chain.get_full_id()[::2]
	sequenceID = headerfmt % {'PDBID':pdbid, 'CHAINID':chainID , 'METHOD':'_ATOM'}
	sequence = reduce(lambda x,y:x+y,   map(lambda x:x.get_sequence(),build_polypeptides(one_chain))  )
	return __SeqRecord(sequence,id=sequenceID,description='')

def build_all_chain_seqs(one_entity):
	chains = __P.Selection.unfold_entities(one_entity,'C')
	retval = dict([(onechain.id, build_chain_seqrec(onechain)) for onechain in chains])
	return retval


def align_chain_to_seq(sequence,chain,verbose=False):
	#Build Polypeptides from the chains
	polypeptides = build_polypeptides(chain)
	
	#Can't be broken out into another function, because we need seq_lens
	contiguous_seqs = [single_pp.get_sequence().tostring() for single_pp in polypeptides]
	ATOM_joined_seq = ''.join(contiguous_seqs)
	
	seq_lens = [0] + [len(single_pp) for single_pp in polypeptides]


	#Figuring all of this out took days...
	#I am so tired of dealing with mapping various numberings around
	#I wish Biopython, especially Bio.pairwise2 had better documentation
	breaks = set(S.cumsum(seq_lens) )#TODO : Tear hair out GYAAAAA
	
	nogaps = lambda x,y: -2000 -200*y #There really should not be inserts with respect to the database sequence.

	def specificgaps(x,y):
		if x in breaks:#very minor penalty for gaps at breaks in the PDB structure, consider using 0
			return (0 -y)
		else:
			return (-2000 -200*y)#strongly discourage gaps anywhere else.
	
	alignments = __PW.align.globalxc(sequence.seq.tostring(),ATOM_joined_seq,nogaps,specificgaps)
	
	if verbose:
		#some output?
		for a in alignments:
			__stderr.write( __PW.format_alignment(*a) )
			__stderr.write('\n')

	return alignments

def aligned_seq_to_res_list(aligned_chain_seq,joined_polypep):
	"""
	Expand a polypeptide or list of Bio.PDB.Residue objects to include None objects for gaps.
	Useful with zip(a,b) to create lists of ('master seq letter', Residue) tuples.
	
	@param aligned_chain_seq: The polypeptide sequence after alignment.
	
	@param polypep: T Bio.PDB.Polypeptide object whence aligned_chain_seq was extracted.
	
	@return: A new representation of the alignment. [None, Res, Res, None, None, ...., etc ]
	
	"""
	anit = joined_polypep.__iter__()
	
	#This nice function generates our list.
	def resornot(letter):
		if letter == '-':
			return None
		else:
			oneres = anit.next()
			return oneres
	
	retlist = [resornot(i) for i in aligned_chain_seq]
	return retlist

__insertionAlphabet = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'
def renumber_chain_to_align(seq_aligned,aligned_res_list,renumber_offset=1):
	"""
	Takes a list of residues as from C{aligned_seq_to_res_list} and the master sequence to which it was aligned (needed to determine inserts, numbering, etc.) And renumbers those residues appropriately.
	
	@precondition: aligned_res_list only contains amino-acids or other covalantly bonded residues in a chain. IE, not water or ligands.
	@postcondition: Residues that were in aligned_res_list are renumbered according to seq_aligned numberings.
	
	@param seq_aligned: The master sequence to which the residues were aligned.
	
	@param aligned_res_list: A list as prepared by C{aligned_seq_to_res_list}.
	
	@param renumber_offset: What number to begin with when renumbering.
	"""
	
	newid = renumber_offset-1
	newinsertid = 0
	for seq_letter, one_residue in zip(seq_aligned,aligned_res_list):
		if seq_letter == '-':#insert w.r.t. master
			newinsertid += 1
		else:				#match or gap w.r.t. master
			newinsertid = 0
			newid += 1

		#TODO: Find a better scheme for handling overly large inserts.
		if newinsertid >= len(__insertionAlphabet):
			raise Exception("Too large an insert!")

		if one_residue != None:
			one_residue.id = (one_residue.id[0],newid,__insertionAlphabet[newinsertid])

	#return

def align_and_renumber_chain(masterseq,chain,renumber_offset=1):
	"""
	Convenience function. Does what it says.
	Uses renumber_chain_to_align, aligned_seq_to_res_list, align_chain_to_seq
	"""
	joinedpolypep = reduce(lambda x,y:x+y, build_polypeptides(chain))
	alignments = align_chain_to_seq(masterseq,chain,verbose=True)
	if len(alignments) > 1:
		raise Exception( "Multiple valid Aligments!" )
	seq_aligned, chain_aligned_seq = alignments[0][:2]
	chain_aligned = aligned_seq_to_res_list(chain_aligned_seq,joinedpolypep)
	renumber_chain_to_align(seq_aligned,chain_aligned,renumber_offset)
	