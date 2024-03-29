# CFreecW : Python Control-Freec Wrapper 

[Control-Freec](http://boevalab.com/FREEC/) is a tool aiming to identify region with copy number variations. 
With this wrapper, I tried to simplify how to produce and analyze CFreec results for simple whole genome data.

## Dependancies

Softwares dependancies :

* [Control-Freec](http://boevalab.com/FREEC/index.html#downloads) (v. >= 11.5). Warning, previous version might not work at all !
* [samtools](http://www.htslib.org/doc/samtools.html)

For the python dependancies :

* Python v. >= 3.6
* [Fire](https://github.com/google/python-fire) : Google CLI
* [Pandas](https://pandas.pydata.org/) : DataFrame analysis
* [Biopython](https://biopython.org/) : Fasta managing
* [PySegs](https://github.com/jsgounot/PySegs) : Bed files comparison

## CFreec configuration file

First you have to think about the [arguments needed by CFreec in order to run correctly](http://boevalab.com/FREEC/tutorial.html#CONFIG). All arguments are splitted in different categories within CFreec : General, sample, control, target. You will find in this repository a default configuration file that you will have to modify. Arguments `samtools` and `mateFile` will be automaticaly added for each sample. For arguments `chrLenFile`, `chrFiles` and the two `ExpectedGC` arguments, see functions below. Take special attention to other arguments, especially `mateOrientation` since wrong value will produce errors.

Note that if the ploidy or the format file (for example, paired end and single end mapping) is not the same between the sample, you will have to run them independantly with different configuration files. Be sure also that you use exactly the same reference than the one use for generating your bamfiles, and that the contig names are the same between you alignment files and your reference.

## CFreecW functions

All functions should be well documented. Use `python cfreecw --help` or for example `python cfreecw.py fastaGC --help` for more informations about available functions.

### Before running Control-Freec :

Two functions are available to make the preparation of your data easier and help you to complete your configuration file :

* `fastaGC` : Compute GC statistics from your reference for the min and max ExpectedGC arguments.
* `makeRef` : Will prepare your reference for CFreec : Indexing the reference and split the contigs in different files (both are printed at the end of the command execution)

### Running Control-Freec :

When your configuration file is completed, you can run Control-Freec using this command : 

`python cfreecw.py runCFreec 'path/to/your/*.bam' config.txt freec samtools --ncore 4`

Other arguments are also available, see `help`. If this command works fine, you should have control freec results for each bamfile in each `path/to/your/bamfile.bam.freec` directory. I recommand to first start with only one sample with ncore set to 1. 

**If the error `FREEC was not able to extract reads` occured** : First off, you can try to run the command line without the wrapper to look at were it failed (the command line should be printed). Then check your `mateOrientation` and `inputFormat` parameters, and finally check that the reference contigs name in your `chrFiles` and `chrLenFile` are exactly the same compared to the one in your bamfile. All these problems can lead to the same error.

### Comparison of your results with a bedfile

Once your data are generated, you can compare regions with CNV (.bam_CNVs files) with annotations formatted within a bedfile using this command line :

`python cfreecw.py compare2bed '/path/to/your/bam.CNVs' /path/to/your/bed ploidy`

Similar to the previous step, you can use multiple cores but try first with only one. This command will compare regions found to your CNVs file to your annotations. By default a threshold of 50% if applied as the minimum coverage requiered for a feature to be determined as impacted by a region, you can change this behavior using the `min_coverage` parameter. You can also add a `fasta` which will add a new column with the distance of each feature with the closest telomere.

Resulting file can be found in the same directory where the related CFreec results file has been given, with the suffix `compared.tsv`.

### CNV matrix

You can concatenate the copy number values of your resulting file into one single matrix using the `make_matrix` function :

`python cfreecw.py make_matrix '/path/to/your/*.bam.compared.tsv' outfile --complete True --ncore 4`

Complete option allows your matrix to contain either all the genes (if false) or only genes which have a gain or a loss in at least one strain. This matrix can be therefore used with [saturation function](https://github.com/jsgounot/PyPan#saturation-curve) available within the PyPan librairy. By default all gene name should be unique, however if multiple genes with the same name are found, an average value for each gene will be calculated.