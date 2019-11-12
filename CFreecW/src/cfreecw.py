# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-07-16 12:39:18
# @Last modified by:   jsgounot
# @Last Modified time: 2019-11-12 19:11:03

import fire
import fanalyzer
import wrapper
import ranalyzer

class FManager(object) :

    def fastaGC(self, fasta, window=20000, step=10000) :
        """
        Extract basic GC content informations from your reference
        
        Arguments:
            fasta {str} -- Fasta file path
            window {int} -- Sliding window to use
            step {int} -- Sliding step to use
        """

        fanalyzer.fasta_analyzer(fasta, window, step)

    def makeRef(self, reference, samtools) :
        """        
        Prepare your reference sequence by :
        - Indexing it using samtools
        - Split your reference chromosomes for CFreec
        
        Arguments:
            reference {str} -- Your reference sequence path (fasta)
            samtools {samtools} -- Samtools path
        """

        wrapper.prepare_reference(reference, samtools)

    def runCFreec(self, input_files, configuration, cfreec, samtools, suffix="", index_bam=False, ncore=1, debug=False) :
        """        
        Run control freec for your different samples
        
        Arguments:
            input_files {str} -- Path of your different samples (use ' if you select multiple files)
            configuration {str} -- Path of your configuration file
            cfreec {str} -- Path to cfreec executable
            samtools {str} -- Path to samtools executable
        
        Keyword Arguments:
            suffix {str} -- Suffix to add to each sample created directory (default: {""})
            index_bam {bool} -- Wether we index your inputfile in case it is a bamfile (default: {False})
            ncore {number} -- Number of core to use (default: {1})
        """

        wrapper.runCFreec(input_files, configuration, cfreec, samtools, suffix, index_bam, ncore)

    def compare2bed(self, cnv_files, bedfile, ploidy, min_coverage=50, fasta=None, ncore=1) :
        """        
        Compare your Control FREEC results to a bedfile
        
        Arguments:
            cnv_files {str} -- Control FREEC result files (use ' if you select multiple files)
            bedfile {str} -- File to your bedfile
            ploidy {int} -- Ploidy of your strain
        
        Keyword Arguments:
            min_coverage {number} -- Minimum coverage for a region to be selected (percent of your feature) (default: {50})
            fasta {str} -- Path to your fasta file if you want to add the minimum telomoeric distance to your output files
            ncore {int} -- Number of core to use (default: {1})
        """

        ranalyzer.compare2beds(cnv_files, bedfile, ploidy, min_coverage, fasta, ncore)

    def make_matrix(self, compared_files, outfile, complete=True, ncore=1) :
        """
        Produce a matrix based on compared files
        
        Arguments:
            compared_files {[str]} -- Compared files produced by the compare2bed function
            outfile {[str]} -- Outfile for the matrix, a tsv file

        Keyword Arguments:
            complete {bool} -- Either if you want that the matrix contains all the genes (True, default) or just genes which have a loss or a gain in at least one strain (False)
            ncore {number} -- Number of core to use (default: {1})
        """

        ranalyzer.make_matrix(compared_files, outfile, complete, ncore)

if __name__ == '__main__' :
    fManager = FManager()
    fire.Fire(fManager)