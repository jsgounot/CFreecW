# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-07-16 13:06:53
# @Last modified by:   jsgounot
# @Last Modified time: 2019-07-16 13:35:59

import argparse
from Bio import SeqIO, SeqUtils
import pandas as pd

def extract_gc_content(fasta, window, step) :
    fdata = SeqIO.parse(fasta, "fasta")
    for record in fdata :
        sequence = record.seq
        for i in range(0, len(sequence), step) :
            subseq = sequence[i:i+window]
            yield SeqUtils.GC(subseq)

def fasta_analyzer(fasta, window, step) : 
    print ("Extract GC content statistics (window : %i, step : %i) ...\n" %(window, step))
    values = (value for value in extract_gc_content(fasta, window, step))
    values = pd.Series(values, dtype=float)
    print (values.describe())