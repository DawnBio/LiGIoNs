# LiGIoNs
LiGIoNs is a profile Hidden Markov Model (pHMM) based method capable of predicting Ligand-Gated Ion Channels (LGICs). The method consists of a library of 10 pHMMs, each corresponding to a single LGIC subfamily. In addition, 14 Pfam pHMMs are used to further annotate and correctly classify unknown protein sequences into one of the 10 LGIC subfamilies.

Avgi E. Apostolakou, Katerina C. Nastou, Georgios N. Petichakis, Zoi I. Litou and Vassiliki A. Iconomidou  
LiGIoNs: Α Computational Method for the Detection and Classification of Ligand-Gated Ion Channels (In Preparation)

The method is available online at http://bioinformatics.biol.uoa.gr/ligions/.  (soon)

# Description
This repository contains a python script PredictLGIC.py and the necessary files to run LiGIoNs locally, a python script profileMaker.py that can be used to create the pHMMs and the predicted LGICs that resulted from the application of LiGIoNs on all UniProtKB reference proteomes (proteome_predicted_LGICs.gz).

# Requirements
*These versions were used during production and testing*  
Python - 3.7.10  
Biopython - 1.78  
Numpy - 1.19.2  
[HMMER](http://hmmer.org/) - 3.3.2

*For recreating the pHMMs the following are also needed:*  
[CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/) - 4.7  
[Clustal Omega](http://www.clustal.org/omega/) - 1.2.4  

# Usage
First make sure that PredictLGIC.py, LGICslib (4 files allLGIChmms.h3*) and PfamLGICslib (pfamHMM directory), if Pfam annotation is wanted, are in the working directory. The protein sequence(s) of interest should be in a file in FASTA format.  

*Run the script against one or more protein sequences in FASTA format:*  
`python PredictLGIC.py -i my_input_filename`  
The output is two files, one with the full results (.lgic) and one with the brief results (.lbrief). By default these files are named after the input file.  
*To specify output filename:*  
`python PredictLGIC.py -i my_input_filename -o my_output_filename`  
