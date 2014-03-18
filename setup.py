#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

import os
def find_all_scripts(topdir):
	allfiles = list()
	for dirpath,subdirs,filenames in os.walk(topdir):
		for onefile in filenames:
			allfiles.append(os.path.join(dirpath,onefile))
	return allfiles

setup(name="PDB_Backmapper",
	version="4.1",
	description="A tool for mapping HMMER HMM columns to exact residues in PDB.",
	author="Bryan Lunt",
	author_email="lunt@cs.ucsd.edu",
	url="http://matisse.ucsd.edu/~lunt/amaranth",
	package_dir={'':'src',"backmapper":"src/backmapper", "backmapper.util":"src/backmapper/util", "backmapper.util.struct":"src/backmapper/util/struct"},
	packages=["backmapper", "backmapper.util", "backmapper.util.struct"],
	scripts=find_all_scripts("scripts"),
	requires=["biopython (>=1.62)", "numpy (>=1.8.0)", "scipy (>=0.11.0)"],
	classifiers=["Intended Audience :: Bioinformaticians"]
)
