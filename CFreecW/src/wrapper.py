# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-07-16 13:40:24
# @Last modified by:   jsgounot
# @Last Modified time: 2019-07-16 18:42:29

import shutil, os, glob, copy
from collections import defaultdict
from Bio import SeqIO
import processing

# -------------------------------
# Reference preparation [makeRef]
# -------------------------------

def prepare_reference(reference, samtools) :
	outdir = reference + ".split4freec"
	if not os.path.isdir(outdir) : os.mkdir(outdir)
	split_reference(reference, outdir)
	chrlenfile = index_reference(reference, samtools)

	print ("makeRef done ...")
	print ("ADD THIS TO YOUR CONFIGURATION FILE :")
	print ("chrFiles =", outdir)
	print ("chrLenFile =", chrlenfile)

def index_reference(reference, samtools) :
	cmdline = "%s faidx %s" %(samtools, reference)
	processing.run_cmdline(cmdline)
	return reference + ".fai"

def split_reference(reference, outdir) :
	fdata = SeqIO.parse(reference, "fasta")
	for record in fdata :
		name = record.name
		outfile = os.path.join(outdir, name + ".fa")
		SeqIO.write([record], outfile, "fasta")

# --------------------------
# Running CFreec [runCFreec]
# --------------------------

class ParsingError(Exception) :
	pass

def parse_cfile(cfile, debug) :
	# Parse the configuration file
	infos = defaultdict(dict)
	cat = None

	for line in processing.readfile(cfile) :

		if debug : print ("DEBUG", line)

		if line.startswith("[") and line.endswith("]") :
			cat = line[1:-1]
			continue

		if cat is None : raise ParsingError("[Category] not found in your configuration file")
		key, value = line.split("=")
		key, value = key.strip(), value.strip()
		infos[cat][key] = value

	return infos

def write_cfile(cdata, outfile) :
	order = ["general", "sample", "control", "target"]

	with open(outfile, "w") as f :
		for cat in order :
			f.write("[%s]\n" %(cat))
			for key, value in cdata[cat].items() :
				f.write("%s = %s\n" %(str(key), str(value)))

def modify_cdata(cdata, ifile, samtools, outputDir) :
	cdata = copy.deepcopy(cdata)
	cdata["general"]["samtools"] = samtools
	cdata["sample"]["mateFile"] = ifile
	cdata["general"]["outputDir"] = outputDir
	return cdata

def runCFreec_sample(ifile, cdata, cfreec, samtools, index_bam, suffix) :
	outdir = ifile + ".freec" + suffix
	if not os.path.isdir(outdir) : os.mkdir(outdir)
	cfile = os.path.join(outdir, "configuration.freec.txt")
	cdata = modify_cdata(cdata, ifile, samtools, outdir)
	write_cfile(cdata, cfile)

	if index_bam :
		cmdline = "%s index %s" %(samtools, ifile)
		processing.run_cmdline(cmdline)

	cmdline = "%s -conf %s" %(cfreec, cfile)
	print ("Run command line : '%s'" %(cmdline))
	processing.run_cmdline(cmdline)

def runCFreec(input_files, configuration, cfreec, samtools, suffix="", index_bam=False, ncore=1, debug=False) :
	ifiles = glob.glob(input_files)
	print ("Work with %i files with %i cores" %(len(ifiles), ncore))

	try : cdata = parse_cfile(configuration, debug)
	except Exception as e : raise ParsingError("An error occured during the configuration file parsing : %s" %(str(e)))

	args = (processing.FunArgs(runCFreec_sample, ifile, cdata, cfreec, samtools, index_bam, suffix) for ifile in ifiles)
	processing.mproc(args, ncore)

