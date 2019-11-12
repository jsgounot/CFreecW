# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-07-16 16:10:41
# @Last modified by:   jsgounot
# @Last Modified time: 2019-11-12 19:14:21

import glob, os
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio import SeqIO
from pysegs import Segment, SegList
import processing

class ParsingError(Exception) :
	pass

class CNVRegion() :

	def __init__(self, start, end, CN, CNK) :
		self.seg = Segment(int(start), int(end))
		self.CN = CN
		self.CNK = CNK

def get_max(row, cnvs, min_coverage, ploidy) :
	column_names = ["CNK", "CN", "coverage"]
	contig = row["chrom"]
	segment = Segment(int(row["chromStart"]), int(row["chromEnd"]))

	if not cnvs[contig] :
		return pd.Series(("Expected", ploidy, np.nan), index=column_names)

	fun = lambda cnv : segment.overlapp_count(cnv.seg, prc=True)
	besthit = max((cnv for cnv in cnvs[contig]), key=fun)
	coverage = segment.overlapp_count(besthit.seg, prc=True)

	if coverage < min_coverage :
		return pd.Series(("Expected", ploidy, np.nan), index=column_names)

	return pd.Series((besthit.CNK, besthit.CN, coverage), index=column_names)

def add_telo_dist(annos, fasta) :
	fdata = SeqIO.parse(fasta, "fasta")
	fdata = {record.name : len(record.seq) for record in fdata}
	fun = lambda row : min(row["chromStart"], fdata[row["chrom"]] - row["chromEnd"])
	annos["min_telo_dist"] = annos.apply(fun, axis=1)
	return annos

def compare2bed(cnv_file, annos, min_coverage, ploidy) :
	print ("Work with result file : %s" %(cnv_file))

	cnvs = parse_result(cnv_file)
	annos = annos.copy(deep=True)

	fun = lambda row : get_max(row, cnvs, min_coverage, ploidy)
	res = annos.apply(fun, axis=1)
	annos = pd.concat((annos, res), axis=1)

	outfile = cnv_file + ".compared.tsv"
	annos.to_csv(outfile, sep="\t")

def parse_result(cnv_file) :
	# columns = ["contig", "start", "end", "cn", "cnk"]
	data = defaultdict(list)
	
	for line in processing.readfile(cnv_file) :
		values = line.split()
		values = [value.strip() for value in values]
		data[values[0]].append(CNVRegion(* values[1:]))

	return data

def parse_bedfile(bedfile) :
	columns = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", 
	"thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]

	df = pd.read_csv(bedfile, sep="\t", names=columns)
	df = df.dropna(axis=1, how="all")

	mandatories = ("chrom", "chromStart", "chromEnd")
	if not all(name in df.columns for name in mandatories) :
		raise ParsingError("Missing required fields in your bedfile (chrom, chromStart and chromEnd), see : https://www.ensembl.org/info/website/upload/bed.html")

	return df

def compare2beds(cnv_files, bedfile, ploidy, min_coverage=50, fasta=None, ncore=1) :
	fnames = glob.glob(cnv_files)
	annos = parse_bedfile(bedfile)
	if fasta is not None : annos = add_telo_dist(annos, fasta)

	args = (processing.FunArgs(compare2bed, fname, annos, min_coverage, ploidy) for fname in fnames)
	processing.mproc(args, ncore)

def read_compared(compared_file) :
	name = os.path.basename(compared_file)[:-22]
	df = pd.read_csv(compared_file, sep="\t", index_col=0)
	df["sample"] = name
	return df

def make_matrix(compared_files, outfile, complete=True, ncore=1) :
	fnames = glob.glob(compared_files)
	args = (processing.FunArgs(read_compared, fname) for fname in fnames)
	df = processing.mproc(args, ncore)
	df = pd.concat(df)

	if not complete :
		genes = set(df[df["CNK"].isin(("gain", "loss"))]["name"].unique())
		df = df[df["name"].isin(genes)]

	df = pd.pivot_table(df, index="name", columns="sample", values="CN")
	df.to_csv(outfile, sep="\t")