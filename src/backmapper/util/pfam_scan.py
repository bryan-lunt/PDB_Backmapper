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
	This object wraps an invocation of the pfam_scan.pl utility.
	'''

	def __init__(self,pfamDBfile,fastafile,pfamScanPath="pfam_scan.pl"):
		'''
		Constructor
		'''
		self.pfamDBfile = pfamDBfile
		self.fastafile = fastafile
		self.pfamScanPath = pfamScanPath
		self.hmmproc = None
		self.returncode = None
		self.byDom = None
		self.bySeq = None
		self.output = None
		self.error = None

	def start(self):
		self.hmmproc = sp.Popen([self.pfamScanPath,"-dir",self.pfamDBfile,"-fasta",self.fastafile,"-json","pretty"],stdout=sp.PIPE,stderr=sp.PIPE)

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