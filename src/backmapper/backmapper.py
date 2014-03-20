#!/usr/bin/python
"""
@author: lunt
"""
import sys, os
from optparse import OptionParser
from itertools import count
from collections import defaultdict

import Bio.SeqIO as SeqIO
import Bio.PDB as PDB

import util.struct as structU
import util.struct.fetcher as structF
import util.struct.contacts as structC

from util.pypfam_scan import pfam_scan

from . import *
from Backmap import *
	
		
def make_backmap(HMMSeqRes,chain_res_list,name=None):
	"""
	
	@param HMMSeqRes: A profile HMM hit to the SEQRES sequence for this chain.
	@type HMMSeqRes: Bio.SeqRecord
	
	@param chain_res_list: A list of [Bio.PDB.Residue or None] entries, as from C{dcaBioUtils.struct.aligned_seq_to_res_list()}
	@type chain_res_list: list
	
	@param name: The name of the backmap, for future use?
	@type name: str
	
	@return: A list of (HMMCol or None, hitChar or '-', Bio.PDB.Residue or None)
	"""

	res_iter = chain_res_list.__iter__()
	
	#RETURN VALUE
	finalBackMap = Backmap(seqrec=HMMSeqRes,name=name)
	
	#skip sufficient residues in chain_res_list to get to the proper starting position of the HMM match
	#This is needed (instead of using slicing) because there are sometimes (are there?) residues in the PDB that are not in the seqres.
	numToPop = HMMSeqRes.annotations['start'] - 1
	while numToPop > 0:
		first_pair = res_iter.next()
		if first_pair[0] != '-':#we only count match states, of course.
			numToPop -= 1
	
	#Begin preparing the backmap
	
	
	hmmPos = count()
	# HMMCol, hitCharacter, Residue
	for hitChar in HMMSeqRes:
		hmmCol = None
		oneResidue = None
		
		#A gap or a match, not an insert.
		#hmmCol remains None for 
		if not hitChar.islower():
			hmmCol = hmmPos.next()
		else:#It's an insert and does not correspond to any HMM column.
			hmmCol = None #Redundant but added for readability. 
		
		#A non-gap (hit or insert), meaning a residue in the PDB file. Those get added to backmaps.
		#if it was a gap in the HMM, then it corresponds to no residue, (But we still output the HMM column ID, handled above.)
		if hitChar != '-':
			seqResChar, oneResidue = res_iter.next()
			assert seqResChar.upper() == hitChar.upper(), "mismatch between seqres and PBD structure"
		
		finalBackMap.append((hmmCol,hitChar,oneResidue))

	return finalBackMap

def pdb_seqres_seqname_toChain(seqname):
	return seqname.split("|")[0].split(":")[1]

def is_no_aa_chain(chain):
	"""
	Test if a chain contains no amino acids.
	"""
	return all([(not PDB.is_aa(r)) for r in chain])

def align_every_chain(structure,fastafile):
	"""
	Reads the fasta file and the structure, Figures out which sequence in the fasta file corresponds to which chain in the structure.
	
	uses dcaBioUtils.struct.align_chain_to_seq and dcaBioUtils.struct.aligned_seq_to_res_list to create reslists for every chain.
	
	
	Returns a Dict : chainID -> [list of alignedSeqs]
	where alignedSeqs is [(seqres_after_alignment_char, Residue) ...]
	
	"""
	#Get all chains from the structure, crate a dictionary CHAINID->Chain
	ALLCHAINS = dict([(i.id, i) for i in structure.get_chains()])
	
	#Read the sequences from the fasta file, create a dictionary CHAINID->Seq
	ALLSEQS = dict([(pdb_seqres_seqname_toChain(i.id),i) for i in SeqIO.parse(fastafile,"fasta")])

	retdict = dict()
	for chainID, thisChainSeqRes in ALLSEQS.iteritems():
		try:
			thechain = ALLCHAINS[chainID]
			if is_no_aa_chain(thechain):
				continue#Skip DNA and RNA chains.
		except:
			raise Exception("Chain ID %s was present in the FASTA file, but not in the structure." % chainID)
		
		
		alignments = structU.align_chain_to_seq(thisChainSeqRes,thechain,verbose=True)
		if len(alignments) > 1:
			raise Exception( "Multiple valid Aligments!" )
		chain_aligned_seqres, chain_aligned_seq = alignments[0][:2]
		
		joinedpolypep = reduce(lambda x,y:x+y, structU.build_polypeptides(thechain))
		chain_aligned = structU.aligned_seq_to_res_list(chain_aligned_seq,joinedpolypep)
		
		retdict[chainID] = zip(chain_aligned_seqres,chain_aligned)
	
	return retdict


def bundle_from_file(pdbfileName,structID,bundleDIR):
	if structID == None:
		structID = os.path.basename(pdbfileName).rsplit(".",1)[0]
	
	structfilename = os.path.join(bundleDIR,"pdb"+structID+".ent")
	import shutil #not always necessary so not in general import
	shutil.copy(pdbfileName,structfilename)
	fastafile = os.path.join(bundleDIR,structID+".faa")

	#dump the SEQRES headers since a fasta file is not available.
	#TODO: What if there are no SEQRES headers? shouldn't you get a sequence from the residues?
	STRUCTFILEHANDLE = open(structfilename)
	sequenceList = structU.read_all_seqres_from_structfile(STRUCTFILEHANDLE,structID).values()
	STRUCTFILEHANDLE.close()
	
	FASTAOUTFILE = open(fastafile,'w')
	SeqIO.write(sequenceList,FASTAOUTFILE,'fasta')
	FASTAOUTFILE.close()
	
	return structfilename, fastafile


def main():
	usage = "Usage: %prog [options] <HmmFile> <PDBID>"
	parser = OptionParser(usage=usage)
	parser.add_option("-d","--dl",dest="download",help="Just download the files and create the mapbundle, but do not create backmaps.",default=False,action="store_true")
	parser.add_option("-f",dest="file",help="Open an extant pdbfile. (PDBID is still used to name the output.)",type="string",default=None)
	options, args = parser.parse_args()

	if not ((options.file and len(args) >= 1) or (len(args) >= 2)):
		parser.print_help()
		parser.exit(1)
	
	#put args into sensible names
	structID = args[1].lower()
	pfamDBfile = args[0]
	pfamDBfile = os.path.join(pfamDBfile, "Pfam-A.hmm")

	#placeholder, later the user will be given the option to specify a bundledir
	
	#TODO update this to have a nice interface for creating backmaps for a file on the file system.!!
	bundleDIR = create_bundle(structID)
	
	if options.file is not None:
		structfile, fastafile = bundle_from_file(options.file,structID,bundleDIR)
	else:
		print "Beginning Download"
		structfile, fastafile = structF.download_PDB_and_FASTA(structID,bundleDIR)
		if options.download:
			sys.exit(0)

	print "Starting pfam_scan.pl"
	pfam_scan_proc = pfam_scan(fastafile, pfamDBfile, save_results=os.path.join(bundleDIR, "hmmscan"))
	pfam_scan_proc.start()

	#parse the structure file, and create some "mappings"	
	aparser = PDB.PDBParser()
	structure = aparser.get_structure("X", structfile)

	ALLCHAINIDS = [i.id for i in structure.get_chains()]

	
	#By this point, the bundle is created, the .faa file was either downloaded, or read from the structure file itself.
	#We can trust it to contain the sequences we need.
	seqResToResidue = align_every_chain(structure,fastafile)
	
	#get the results of pfam_scan.pl
	if pfam_scan_proc.wait() is not 0:
		print "HMMSEARCH FAILED: "
		print pfam_scan_proc.return_stderr()
		sys.exit()
	else:
		print "HMMSearch Success"
	
	pfam_scan_proc.get_results()
	pfam_scan_proc.filter_clan_overlap()
	
	
	import pdb
	pdb.set_trace()

	
	#finally create the actual backmaps. UGH!
	backmaps = defaultdict(list)
	for one_query_result in pfam_scan_proc.query_results:
		chainID = pdb_seqres_seqname_toChain(one_query_result.id)
		for one_hit in one_query_result.hits:
			hit_length = one_hit.seq_len
			for one_hsp in one_hit.hsps:
				start, end = one_hsp.hit_range
				oneSeqRecord = "-"*start + one_hsp.query + "-"*(hit_length-end)
				oneSeqRecord.annotations['start'] = one_hsp.query_start+1
				
				backmaps[chainID].append(make_backmap(oneSeqRecord, seqResToResidue[chainID], name=one_hit.id))
	


	#output backmaps
	for mapChainID, curChainMaps in backmaps.iteritems():#output maps to files
		#domNo = 0#domains come pre-sorted by hmmpfam, we will count up, because we want to order them
		for abmap, domNo in zip(curChainMaps,count(1)):
			
			mapfilename = '%s.%i.%s.map' % (mapChainID,domNo,abmap.name)
			
			outfile = open(os.path.join(bundleDIR,mapfilename),"w")
			abmap.write(outfile)
			outfile.close()

	#Distance Map File
	dMap = structC.DistanceMap(structure[0])
	distfileHeader = "#resID1	chain1	resID2	chain2	atom1	atom2	distance\n"
	distfileFmt = '%(iid)s\t%(ichain)s\t%(jid)s\t%(jchain)s\t%(iatom)s\t%(jatom)s\t%(dist).4f\n'
	
	distoutfile = open(os.path.join(bundleDIR,"distances"),"w")
	dMap.write(distoutfile,formatstring=distfileFmt,header=distfileHeader)
	distoutfile.close()
	
	#output renumbered chains and renumbered domains.
	saver = PDB.PDBIO()
	saver.set_structure(structure)
	#Save renumbered chains
	for chainID, seqResMapping in seqResToResidue.iteritems():
		chain_aligned_seqres = ''.join([i[0] for i in seqResMapping])
		aligned_residues     = [i[1] for i in seqResMapping]
		
		structU.renumber_chain_to_align(chain_aligned_seqres,aligned_residues,0)
		
		class __MySelector(PDB.Select):
			def accept_chain(self,chain):
				return chain.id == chainID
		
		onechainOutFile = open(os.path.join(bundleDIR,'Chains',chainID+'.pdb'),'w')
		saver.save(onechainOutFile,__MySelector())
		onechainOutFile.close()
	
	#Save renumbered domains
	for mapChainID, curChainMaps in backmaps.iteritems():
		for abmap, domNo in zip(curChainMaps,count(1)):
			hmmSeq = ''.join(['-' if i[0] is None else 'X' for i in abmap])
			reslist = [i[2] for i in abmap]
			resset = set(reslist)
			
			structU.renumber_chain_to_align(hmmSeq,reslist,0)
			
			class __MySelector(PDB.Select):
				#def accept_chain(self,chain):
				#	return chain.id == mapChainID
				def accept_residue(self,residue):
					return residue in resset
			
			mapfilename = '%s.%i.%s.pdb' % (mapChainID,domNo,abmap.name)
			oneDomOutFile = open(os.path.join(bundleDIR,'Domains',mapfilename),'w')
			saver.save(oneDomOutFile,__MySelector())
			oneDomOutFile.close()
			
			

if __name__ == "__main__":
	main()
