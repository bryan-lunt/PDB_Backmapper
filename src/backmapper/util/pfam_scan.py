'''
Created on Jun 28, 2010

@author: lunt
'''
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess as sp
import json


class pfam_scan(object):
	'''
	This object wraps an invocation of the pfam_scan.pl utility using the "subprocess" library.
	
	It should be replaced with Bio.SearchIO.HmmerIO . Its advantage is that by using the actual PFAM code, we get exactly their results...
	
	'''

	def __init__(self,pfamDB_dir,fastafile,pfamScanPath="pfam_scan.pl"):
		'''
		Constructor.
		
		@param pfamDB_dir: The directory containing the files "Pfam-A.hmm*" and "Pfam-A.hmm.dat" for pfam_scan.pl to use.
		@type pfamDB_dir: str
		
		@param fastafile: The FASTA format file containing unaligned sequences to scan with pfam_scan.pl
		@type fastafile: str
		
		@param pfamScanPath: The path (or partial path, since PATH resolution is used) to the PFMA scan utility, usually "pfam_scan.pl"
		@type pfamScanPath: str
		
		'''
		self.pfamDBfile = pfamDB_dir
		self.fastafile = fastafile
		self.pfamScanPath = pfamScanPath
		self.hmmproc = None
		self.returncode = None
		self.byDom = None
		self.bySeq = None
		self.output = None
		self.error = None

	def start(self):
		"""
		Start a background process for scanning.
		
		"""
		cmdline = [self.pfamScanPath] #The program
		cmdline += ["-dir",self.pfamDBfile] #Where Pfam-A.hmm is stored
		cmdline += ["-fasta",self.fastafile] #which Fasta file to invoke on
		cmdline += ["-json","pretty"] #We want pretty-printed JSON format output.
		
		#Create the invocation, we want to be able to read STDOUT and STDERR back in.
		self.hmmproc = sp.Popen(cmdline,stdout=sp.PIPE,stderr=sp.PIPE)

	def wait(self):
		self.output,self.error = self.hmmproc.communicate()
		self.returncode = self.hmmproc.returncode
		return self.returncode

	def poll(self):
		self.returncode = self.hmmproc.poll()
		return self.returncode

	def return_stderr(self):
		return self.error

	def parse_output(self):
		allHits = json.loads(self.output)
		
		self.byDom = dict()
		self.bySeq = dict()
		
		for aHit in allHits:
			hmmName = str(aHit["name"])
			name = str(aHit["seq"]["name"])
			alignedSeq=str(aHit["align"][-1].split()[1])
			seqStart=int(aHit["seq"]["from"])
			seqEnd=int(aHit["seq"]["to"])
			hmmStart=int(aHit["hmm"]["from"])-1#convert to 0-index
			alignedSeq = "-"*hmmStart + alignedSeq + "-"*(int(aHit["model_length"])-int(aHit["hmm"]["to"])) #pad the alignment with the necessary number of gaps, because HMMER3 does not do global alignment.
			aSeq = Seq(alignedSeq)#to eliminate unicode strings
			aSeqRec = SeqRecord(aSeq,name=name+"/%i-%i"%(seqStart,seqEnd))
			aSeqRec.annotations["start"] = seqStart
			aSeqRec.annotations["end"] = seqEnd
			aSeqRec.annotations["region"] = (seqStart,seqEnd)
			
			byDomList = self.byDom.setdefault(hmmName,list())
			byDomList.append((name,aSeqRec))
			bySeqList = self.bySeq.setdefault(name,list())
			bySeqList.append((hmmName,aSeqRec))