-------------
-------------
DivA 1.0
M. Lisandra Zepeda Mendoza & Rute R. da Fonseca
-------------
-------------


-------------
DESCRIPTION
-------------

Set of python scripts designed to detect non-homologous and very Divergent regions in protein sequence Alignments. DivA was tested with python 2.7

DivA makes no assumptions on evolutionary models, and it is ideal for detecting incorrectly annotated segments within individual gene sequences. DivA is a python script that is a binary decision making method that inapplies a sliding-window approach to estimates four divergence-based parameters and defines their outlier values according to automatically defined thresholds that can be optionally modified. DivA then classifies the windows of a sequence of an alignment as very divergent (potentially non-homologous) if it presents a combination of outlier values for the four parameters. The windows classified as very divergent can optionally be masked in the alignment.  This allows DivA to discard a minimum amount of sequence information compared to other currently available methods that remove entire sequences or blocks of a multiple sequence alignment. One important application of DivA is in the detection of incorrect automatic gene annotated sequences, which can have confounding effects in comparative genomics and phylogenomics analyses.


-------------
INSTALLATION
-------------

DivA is a python script that does not need any sort of compilation. It was developed in Python 2.7.3 and uses the following modules which should be already installed in the user's system:

- numpy
- function AlignIO from module Bio
- re
- os
- sys
- argparse

Make sure to put the bin in your path, where the blosum62.txt should also be placed; alternatively place the blosum62.txt or another distance matrix of preference on the same directory where DivA is going to me used.


------
USAGE
------

usage: DivA.py [-h] [--mask] [--printAllwindows] [-w W] [-g G] [-p P] [-zp ZP]
               [-d D] [-zd ZD] [-o O] [-m M]
               alnNamesFile

Identify very divergent potentially non-homologous windows in a protein
multiple sequence alignment.

positional arguments:
  alnNamesFile       A txt file with the file name(s) of the MSA(s) on which
                     to perform the method

optional arguments:
  -h, --help         show this help message and exit
  --mask             Flag for the output of an alignment with the wrong
                     windows masked with XXs [default not set]
  --printAllwindows  Flag for the output of a file with the parameter values
                     and start and end positions of all the windows in the
                     MSA(s) [default not set]
  -w W               The size of the sliding window [default 12]
  -g G               Maximum gap content in a window to be considered [default
                     0.6]
  -p P               The number of standard deviations from the mean of the
                     alpha parameter to use as threshold [default 1]
  -zp ZP             The number of standard deviations from the mean of the
                     Zalpha parameter to use as threshold [default 2]
  -d D               The number of standard deviations from the mean of the
                     beta parameter to use as threshold [default 2]
  -zd ZD             The number of standard deviations from the mean of the
                     Zbeta parameter to use as threshold [default 2]
  -o O               Output basename prefix [default "out"]
  -m M               The amino acid distance matrix [default "blosum62.txt"]





#Example:


 1. Create a file with the names/paths of the alignments to be analyzed. The final thresholds will be calculated using all those alignemnts.

 2. Run DivA:

python DivA_RF.py ListOfAlignments.txt #Basic default DivA run

python DivA.py -h # Will display the help

python DivA.py ListOfAlignments.txt -o DivaOutput --mask --printAllwindows # The outputs will have the prefix "DivaOutput" and alignments with the wrong windows masked wll begenerated, as well as an etra output file containing all the windows with the four parameter values and start and end positions.

python DivA.py ListOfAlignments.txt -o DivaOutput -p 2 # The number of standard deviations form the mean of the alpha parameter is changed to 2 and the outputs will have the prefix "DivaOutput"


#Example files in the 'Test' directory
The Test.aln file corresponds to ortholog alignment 14518.fasta from Jarvis et al.


-----
CITE
-----

Zepeda Mendoza ML, Nygaard S, and da Fonseca R (2014)  "DivA: detection of non-homologous and very Divergent regions in protein sequence Alignments"

--------
CONTACT
--------

For any enquiries correspondence is sent to rute.r.da.fonseca@gmail.com



