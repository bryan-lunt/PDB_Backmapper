'''
Created on Mar 18, 2014

@author: Bryan Lunt <blunt@cs.ucsd.edu>
'''

import Bio.SearchIO as SIO
import Bio.SearchIO.HmmerIO as HIO

import collections
from re import split as _re_split

from itertools import count
import os
import tempfile
import subprocess as sp

def _HMMIO_joiner(domain_table_results, text_results):
	"""
	Function to update into text_results information that is only available in domain_table_results. (Specifically, the length of the HMM.)
	
	 
	"""
	#filter the textHits so they will match the domain table output, (which does not give any output for sequences with no hit.)
	filtered_textHits = filter(lambda x: len(x.hits) > 0, text_results)
	
	assert len(domain_table_results) == len(filtered_textHits), "Bryan Lunt is an idiot who makes too many assumptions. Please complain to him and include your input."
	for domQR, textQR in zip(domain_table_results, filtered_textHits):
		assert len(domQR.hits) == len(textQR.hits), "Bryan Lunt is an idiot who makes too many assumptions. Please complain to him and include your input."
		for domHit, textHit in zip(domQR.hits, textQR.hits):
			textHit.seq_len = domHit.seq_len

	return filtered_textHits

class pfam_datfile(object):
	"""
	Class to parse and otherwise interact with a Pfam-A.hmm.dat file.
	
	This file is a regular language, so it's pretty easy.
	"""
	
	converters = {'GA': lambda val: tuple(map(
									lambda x: float(x.strip()),
									filter(lambda x:x!='', val.split(';'))
									)),
				'ML': lambda val:int(val.strip())}
	
	class PfamFamily(collections.defaultdict):
		
		def __init__(self):
			super(pfam_datfile.PfamFamily, self).__init__(lambda : None)
			self.__dict__ = self
		
		def __hash__(self):
			return hash(self.get('AC', ""))
		
		def __cmp__(self,other):
			return cmp(self.get('AC', ""), other.get('AC',''))
		
		def __str__(self):
			return self.get('ID',"") + ": " + self.get('AC', "") + ": " + self.get('DE','')
		
		def __repr__(self):
			return "PfamFamily(" + self.get('ID',"") + ", " + self.get('AC', "") + ", " + self.get('DE','') + ")"
	
	def __init__(self):
		self.clans = collections.defaultdict(set) #A dictionary of sets, clan-name -> Pfam Accession
		self.accs = dict() #A dictionary of accessions -> PfamFamily
		self.ids = dict() #A dictionary of IDs -> PfamFamily
		self.as_list = list()
	
	def __add_family(self, newFam):
		if newFam.get('CL', None) is not None:
			self.clans[newFam.CL].add(newFam)
		self.accs[newFam.AC] = newFam
		self.ids[newFam.ID] = newFam
		self.as_list.append(newFam)
	
	def load(self, datfile):
		"""
		Finite state machine with two states : idle -> reading_attributes (with a self-loop) -> idle
		"""
		state = 0 # zero for idle, 1 for reading
		current = None
		for line in datfile:
			line = line.strip()
			if state == 0:
				if not line == "# STOCKHOLM 1.0":
					print line
					raise Exception("Expected beginning of entry, got:" + line)
				current = pfam_datfile.PfamFamily()
				state = 1
			elif state == 1:
				if line == "//":
					self.__add_family(current)
					current = None
					state = 0
				else: #read one property
					if not line.startswith("#=GF "):
						raise Exception("Problem parsing, expected GF entry.")
					attribute_name, attribute_value = _re_split("\s*", line[5:], 1)
					
					
					converted_value = pfam_datfile.converters.get(attribute_name, lambda x:x).__call__(attribute_value)
					current[attribute_name] = converted_value
		
		if state != 0:
			raise Exception("Reached end of file while still expecting input")
				
			

class pfam_scan(object):
	"""
	To use this class:
	instantitate it.
	call start()
	call wait()
	call get_results()
	call filter_clan_overlap() if you want.
	get the results from pfam_scan.query_results
	"""
	
	
	def __init__(self, fastafile, pfamDB_file, pfamDB_dat_file=None, save_results=None, hmmscan_path="hmmscan"):
		'''
		Constructor.
		
		@param fastafile: The FASTA format file containing unaligned sequences to scan with pfam_scan.pl
		@type fastafile: str
		
		@param pfamDB_file: The path to "Pfam-A.hmm" for hmmscan to use.
		@type pfamDB_file: str
		
		@param pfamDB_dat_file: The path to "Pfam-A.hmm.dat" for hmmscan to use, if left as None, it will be inferred automatically from pfamDB_file .
		@type pfamDB_dat_file: str
		
		@param save_results: A file prefix to use to save the output of hmmscan. If None, a temporary directory is created to hold this as long as needed.
		@type save_results: str
		
		@param hmmscan_path: The path (or partial path, since PATH resolution is used) to the "hmmscan" utility, usually "hmmscan"
		@type hmmscan_path: str
		
		'''
		self.fastafile = fastafile
		self.pfamDB_file = pfamDB_file
		if pfamDB_dat_file is None:
			pfamDB_dat_file = pfamDB_file + ".dat"
		self.pfamDB_dat_file = pfamDB_dat_file
		self.pfamDB_dat = None #We'll parse this later
		self.pfamScanPath = hmmscan_path
		
		
		if save_results is None:
			save_results = os.path.join(tempfile.mkdtemp(),"hmmscan_results")
		self.save_results = save_results
		
		self.hmmproc = None
		self.returncode = None
		self.byDom = None
		self.bySeq = None
		self.output = None
		self.error = None
		
		self.query_results = list()
	
	def start(self):
		"""
		Start a background process for scanning.
		
		"""
		cmdline = [self.pfamScanPath] #The program
		cmdline += ["-o", self.save_results + ".out", "--domtblout", self.save_results + ".domtblout"]
		cmdline += ['--cut_ga'] #See PfamScan/Bio/Pfam/Scan/PfamScan.pm 
		cmdline += [self.pfamDB_file] #Where Pfam-A.hmm is stored
		cmdline += [self.fastafile] #which Fasta file to invoke on
		
		#Create the invocation, we want to be able to read STDOUT and STDERR back in.
		self.hmmproc = sp.Popen(cmdline,stdout=sp.PIPE,stderr=sp.PIPE)
		
		#now parse the dat file simultaneously
		self.pfamDB_dat = pfam_datfile()
		with open(self.pfamDB_dat_file) as DATFILE:
			self.pfamDB_dat.load(DATFILE)
	
	def wait(self):
		self.output,self.error = self.hmmproc.communicate()
		self.returncode = self.hmmproc.returncode
		return self.returncode
	
	def poll(self):
		self.returncode = self.hmmproc.poll()
		return self.returncode
	
	def return_stderr(self):
		return self.error
	
	def get_results(self):
		"""
		Parse the results of "hmmscan" into self.query_results
		"""
		if self.returncode != 0:
			raise Exception( "You can only get results if the scan succeeded!")
		
		with open(self.save_results + ".domtblout") as DOMTBLFILE:
			with open(self.save_results + ".out") as TEXTFILE:
				domTblHits = [i for i in HIO.Hmmer3DomtabHmmhitParser(DOMTBLFILE)]
				textHits = [i for i in HIO.Hmmer3TextParser(TEXTFILE)]
		
		#filter useless hits from textresults (below reporting threshold per-domain)
		textHits = [i.hit_filter(lambda x:len(x.hsps) > 0) for i in textHits]
		
		#Copies information unavailable in the HMMer text format. Specifically, the length of the HMM.
		_HMMIO_joiner(domTblHits, textHits)
		
		#return the full textHits, which may have some queries with no Hits
		self.query_results = textHits
		#TODO: test our assumption that HMMER only gives one fragment per hit?
	
	def filter_clan_overlap(self):
		"""
		Remove overlapping hits that are in the same clan.
		"""
		self.query_results = map(lambda x:_single_filter_clan_overlap(x,self.pfamDB_dat), self.query_results)
	
	
def _single_filter_clan_overlap(single_query_result, pfam_dat_object):
	"""
	Removes overlapping hits _within_ clans.
	
	@param single_query_result: A query result to filter.
	@type single_query_result: Bio.SearchIO.QueryResult
	
	@param pfam_dat_object: A representation of the "Pfam-A.hmm.dat" file, used to separate clan overlap.
	@type pfam_dat_object:
	
	@return: The original query with overlapping HSPS from the same clan removed.
	@rtype: Bio.SearchIO.QueryResult
	"""
	sorted_hsps = single_query_result.hsps
	sorted_hsps.sort(key=lambda x:x.evalue)
	
	accepted_hsps = list()
	
	clan_hsps = collections.defaultdict(list) #Clan -> list. we have already sorted all HSPS, so these sub-lists are also sorted.
	
	#accept any hits that are not members of a clan
	#sort all clan members according to clan
	#BEGIN SORT
	for one_hsp in sorted_hsps:
		clan = None
		datEntry = pfam_dat_object.ids.get(one_hsp.hit_id,None)
		if datEntry is not None:
			clan = datEntry.get('CL',None)
		
		if clan is None:
			accepted_hsps.append(one_hsp)
		else:
			clan_hsps[clan].append(one_hsp)
	#END SORT
	
	#anonymous inner function to test overlap
	def test_query_overlap(left,right):
		return not (left.query_end <= right.query_start or left.query_start >= right.query_end)
	
	#For each clan, accept non-overlapping hits greedily by evalue.
	for one_clan in clan_hsps.itervalues():
		this_clan_accepted = list()
		
		for clan_member in one_clan:
			if all(map(lambda x:not test_query_overlap(clan_member, x), this_clan_accepted)):
				this_clan_accepted.append(clan_member)

		accepted_hsps.extend(this_clan_accepted)
	
	#filter out HSPs that are not in the accepted set.
	return single_query_result.hsp_filter(lambda x:x in accepted_hsps)

