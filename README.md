# LiGIoNs
LiGIoNs is a profile Hidden Markov Model (pHMM) based method capable of predicting Ligand-Gated Ion Channels (LGICs). The method consists of a library of 10 pHMMs, each corresponding to a single LGIC subfamily. In addition, 14 Pfam pHMMs are used to further annotate and correctly classify unknown protein sequences into one of the 10 LGIC subfamilies.

*Apostolakou AE, Nastou KC, Petichakis GN, Litou ZI, Iconomidou VA.*   
**LiGIoNs: A computational method for the detection and classification of ligand-gated ion channels.**   
Biochim Biophys Acta Biomembr. 2022 Sep 1;1864(9):183956. doi: 10.1016/j.bbamem.2022.183956. Epub 2022 May 14. PMID: 35577076.

The method is available online at http://bioinformatics.biol.uoa.gr/ligions/.  

# Contents
*This repository contains:* 
- **PredictLGIC.py** and the necessary files to run LiGIoNs locally
- **profileMaker.py** that can be used to create the pHMMs and the resulting directories
- **negative test sets** used for evaluation
- **predicted LGICs** from the application of LiGIoNs on all **UniProt reference proteomes**

# Required software
*These following software tool versions were used during production and testing:*  
Python - 3.7.10   
[NumPy](https://numpy.org/) - 1.19.2  
[Biopython](https://biopython.org/) - 1.78  
[HMMER](http://hmmer.org/) - 3.3.2

*For recreating the pHMMs the following are also needed:*  
[CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/) - 4.7  
[Clustal Omega](http://www.clustal.org/omega/) - 1.2.4  

*(Operating system used: Ubuntu - 18.04)*
# Usage
First make sure that **PredictLGIC.py**, **LGICslib** (4 files allLGIChmms.h3*) and **PfamLGICslib** (pfamHMM directory), if Pfam annotation is wanted, are in the working directory. The protein sequence(s) of interest should be in a file in **FASTA** format.  

*Run the script against one or more protein sequences in FASTA format:*  
`python PredictLGIC.py -i my_input_filename`  
The output is two files, one with the full results (**.lgic**) and one with the brief results (**.lbrief**). By default these files are named after the input file.  
*To specify output filename:*  
`python PredictLGIC.py -i my_input_filename -o my_output_filename`  

## Recreating the pHMMs
In the directory **LGIC_fasta** are 10 multiFASTA files, one per LGIC subfamily. The script **profileMaker.py** uses **CD-HIT**, **Clustal Omega** and **HMMER** to create **LGICslib**, a pHMM library containing 1 pHMM per subfamily.

# Application in reference proteomes
The archive **proteome_predicted_LGICs.gz** contains the predicted LGICs from the application of LiGIoNs on all UniProt reference proteomes. Each file in the archive corresponds to a single proteome, is named according to the UniProt **Proteome ID** and has 3 columns, the LGIC subfamily name, the score and the protein information. 
