# LiGIoNs
LiGIoNs is a profile Hidden Markov Model based method capable of predicting Ligand-Gated Ion Channels (LGICs). The method consists of a library of 10 pHMMs, each corresponding to a single LGIC subfamily. In addition, 14 Pfam pHMMs are used to further annotate and correctly classify unknown protein sequences into one of the 10 LGIC subfamilies.

Avgi E. Apostolakou, Katerina C. Nastou, Georgios N. Petichakis, Zoi I. Litou and Vassiliki A. Iconomidou  
LiGIoNs: Α Computational Method for the Detection and Classification of Ligand-Gated Ion Channels (In Preparation)

The method is available online at http://bioinformatics.biol.uoa.gr/ligions/.  (soon)

# Description
This repository contains a python script PredictLGIC.py and the necessary files to run LiGIoNs locally, a python script profileMaker.py that can be used to create the HMM profiles and the predicted LGICs that resulted from the application of LiGIoNs on all UniProtKB reference proteomes (proteome_predicted_LGICs.gz).
