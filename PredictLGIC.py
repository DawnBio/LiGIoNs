#!/usr/bin/env python3
##
## PredictLGIC
##
## Authors:  
##           Final version         Avgi E. Apostolakou (Dawn)
##           Original version      Georgios N. Petichakis
##
## Input:    FASTA file
## Output:   Results file
##           Brief results file (Optional)
##
## Web tool:
##           http://thalis.biol.uoa.gr/ligions/     
##

import re
import os, glob,subprocess
from Bio import SeqIO
from Bio import SearchIO
import argparse

#Parse Command Line arguments
#-i, --input   -> Required FASTA file with sequence(s) of interest
#-o, --output  -> Optional result files name, default uses input file name
parser = argparse.ArgumentParser(description='FASTA file required as INPUT.')
parser.add_argument('-i','--input', help='Input file name',required=True) #required input file argument
parser.add_argument('-o','--output', help='Output file name') #optional output file name
args = parser.parse_args()
proteinInputFile = args.input #get input file name (or path if not in working directory)
proteinInputFileWithoutPath =  proteinInputFile.split("/")[-1] #get name of file without path (e.g. 'path/to/file')
proteinInputFileWithoutExt = proteinInputFileWithoutPath.split(".")[0] #get name of file without exstension (e.g. '.fa')

#Need to change if HMMER installed elsewhere
hmmscan_path = "/usr/local/bin/hmmscan" 
hmmsearch_path = "/usr/local/bin/hmmsearch"

#pHMM libraries
LGICslib = "allLGIChmms" #Required for LGIC prediction, name of LGICslib composed of 4 files (.h3*)
pfamPath = "pfamHMM/" #Optional, required for Pfam annotation

#Results directory
path = 'LGIC_res/' #Output directory name (or path to the directory)
if not os.path.isdir(path): #check if directory for results does not exist
    os.makedirs(path) #create directory

#Set output file names
if args.output: #Check if (optional) output file name is provided
    outfile = "%s/%s.lgic" % (path,args.output) #file to write detailed results for all sequence
    outfileBrief = "%s/%s.lbrief" % (path,args.output) #file to write brief results for all sequences
else: #otherwise use the input file name
    outfile = "%s/%s.lgic" % (path,proteinInputFileWithoutExt) #file to write detailed results for all sequence
    outfileBrief = "%s/%s.lbrief" % (path,proteinInputFileWithoutExt) #file to write brief results for all sequences


#Sequence score thresholds per LGIC subfamily used for LGIC prediction
thresh_fam = {}
thresh_fam['5HT3']= 177.5
thresh_fam['ASICs']= 117.5
thresh_fam['ENaCs']= 90
thresh_fam['GABAA']= 117.5
thresh_fam['iGluR']= 147.5
thresh_fam['Glycine']= 152.5
thresh_fam['IP3R']= 295
thresh_fam['nAchR']= 62.5
thresh_fam['P2XR']= 125
thresh_fam['ZACs']= 142.5

#Helper function that prints a FASTA sequence to file
def proteinToFasta(prtn,file): #prtn=FASTA seq, file=output file
    outfile = file
    handle = open(outfile, "w")
    SeqIO.write(prtn, handle, "fasta") #use SeqIO function to write FASTA to file
    handle.close()

#Helper function to print a sequence 60 characters per line
def seqOut(outHandler,seq,ty): #outHandler=file handle, seq=sequence to print, ty= type of sequence for label
    out = "%s\t%-5s\n" % (ty, len(seq)) #first line with sequence label and length
    outHandler.write(out) #print first line
    x=0 #needed for sequence<60
    for x in range(0,int(len(seq)/60)): #for every 60 aa in sequence
        out = "  %5s %s %-5s\n" % (x*60+1,seq[x*60:(x+1)*60],(x+1)*60) #get these 60 aa
        outHandler.write(out)
    out = "  %5s %-60s %-5s\n" %  ((x+1)*60+1,seq[(x+1)*60:],len(seq)) #get the remaining 60 aa or till sequence end
    outHandler.write(out)

#Helper function to get the maxScore from hmmsearch results
def getBestScore(infileHH): #infileHH=output of hmmsearch
    hmmSearchOut = list(infileHH.readlines())
    maxScore = -100 #initialize maxScore
    for x in hmmSearchOut:
        if x.strip().startswith("=="): #find lines where each domain details start
            test = x.strip().split() #to get score for domain
            if maxScore < float(test[4]):
                maxScore = float(test[4])
    return 	maxScore #return the highest domain score

#Helper function to get the coordinates for the aligned sequence region
def getRegions(infileHH,bestScore): #infileHH=output of hmmsearch, bestScore=highest domain score
    bss = str(bestScore) #bestScore into string for use in pattern
    lines = infileHH.readlines()
    for x in range(0,len(lines)):
        line = lines[x].rstrip() #remove spaces from end of string
        #make pattern to find line with highest scoring domain information and capture start and end of alignment
        #   2 !   42.8   1.2   2.1e-15   2.1e-15       2      28 .]     239     265 ..     238     265 .. 0.98
        pattern1 = "\s+\d+\s+[!\?]\s+%s\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)" % bss 
        m = re.search(pattern1,line,re.U) #look for the pattern
        if(m): #if a match is found
            aliStart = m.group(1) #captured Start
            aliEnd = m.group(2) #captured End
            return (aliStart,aliEnd) 
    return (0,0) #if no match is found return (0,0)

#Main function responsible for the prediction
def predict(infile,outfile): #infile=FASTA input, outfile=detailed results
    hhl = open(infile,'r')
    ptr1 = list(SeqIO.parse(hhl, "fasta")) #make a list of FASTAs
    hhl.close()
    nSeq = 0 #sequence counter
    with open(outfile,'w') as ho:
        for x in ptr1: #go through all sequences
            #initializing necessary variables
            inside = 0 #to check when in subfamily block of hmmscan results
            indom = 0 #to check when in domain block of hmmscan results
            di = 0 #counter for domains
            isLGIC = 0
            allhits = [] #all hits from hmmscan
            dom_ann = [] #annotation for each domain of hit for results file
            que_S = [] #query start
            que_E = [] #query end
            que_E.append(0) #necessary start
            
            if not x.id or x.id.isspace(): #in case there is FASTA ID
                x.id = "Sequence%d" % (nSeq) #set ID to SequenceN
            print("\nSequence ", x.id)
            fileName = infile.split("/")[-1] #get name of file without path (e.g. 'path/to/file')
            seq_path= "%s%s.single" % (path,fileName) #path for sequence file
            proteinToFasta(x,seq_path) #save sequence FASTA to file
            ho.write("ID %s\n" % x.id) #Write Results: ID line
            search_path = "%s%s.tmp.txt" % (path,fileName) #temporary storage of hmmscan and hmmsearch results
            #Use hmmscan to scan query sequence against LGICslib (library of pHMMs one per subfamily)
            try:
                run_hmmscan = subprocess.run([hmmscan_path, "-o", search_path, LGICslib, seq_path])
                if run_hmmscan.returncode != 0:
                    raise Exception()
            except Exception  as ex: #print error message and exit if hmmscan fails
                print("hmmscan could not proceed, ensure HMMER is installed and the path is correct.\n More details:\n")
                print(ex)
                exit()
            for record in SearchIO.parse(search_path, 'hmmer3-text'): #in the context of hmmscan record=query protein and hit=pHMM
                hits = record.hits 
                num_hits = len(hits)
                if num_hits > 0: # if there are any hits 
                    for i in range(0,num_hits):
                        allhits.append(hits[i].id +'('+str(hits[i].bitscore)+')') #keep all hits and the corresponding scores for Results
                        if isLGIC==0 and hits[i].bitscore >= thresh_fam[hits[i].id]: #if another higher scoring pHMM was not already assigned and the current pHMM score is above threshold
                            best = hits[i].id #name of pHMM/LGIC subfamily
                            doms = hits[i].domain_obs_num
                            isLGIC = 1
                            alisc = [[] for i in range(doms)] #to fill with alignment lines from hmmscan results file later
                            for a in range (0, doms):
                                que_S.append(hits[i].hsps[a].query_start)
                                que_E.append(hits[i].hsps[a].query_end)
                                dom_ann.append("TP   %4s   %15s   %5s   %5s   %5s\n" % ("DOM",best+"_"+str(a+1),str(hits[i].hsps[a].query_start+1),str(hits[i].hsps[a].query_end),str(hits[i].hsps[a].bitscore))) #create annotation for domains, start and end based on query sequence and domain score
            #if query protein is predicted as an LGIC proceed to annotation 
            if isLGIC: 
                print("is predicted as LGIC in subfamily ", best)
                with open(search_path, 'r') as reader: #open hmmscan result file to extract alignment information
                    for line in reader.readlines():
                        if re.search("\A>> "+best, line): #find the correct subfamily/hit results
                            inside = 1
                            indom = 0
                            continue
                        elif re.search("\A>> ", line) or re.search("\AInternal pipeline", line): #when entering next hit or end of results
                            inside = 0
                            indom = 0
                            continue
                        if re.search("  == domain", line) and inside == 1: #when we enter the domain line of the correct subfamily/hit
                            indom = 1
                            di = di+1
                            continue
                        if indom ==1: 
                            alisc[di-1].append(line) #keep the alignment and annotation (per 5 lines, 0=pHMM, 1=alignment annotation, 2=query sequence, 3=posterior probability, 4=newline)
                annoF = '' #alignment annotation from hmmscan, is processed so it matches gapless query sequence
                #Quote from HMMER User's Guide (http://eddylab.org/software/hmmer/Userguide.pdf) "The midline indicates matches between the query profile and target sequence. A letter indicates an exact match to the profile consensus. A + indicates that this residue has a positive log-odds emission score, a “conservative substitution” given what the profile expects at that position."
                seqF = '' #gapless query sequence, lowercase=position not aligned to pHMM, uppercase=position aligned to pHMM
                for d in range(0, doms): #to join sequence segments before, inbetween and after domains
                    seq=''
                    anno=''
                    for h in range(0, len(alisc[d])-4,5): #every 5 lines (=1 segment) of alignment results
                        a = alisc[d][h+2] #query sequence
                        xx = re.search("\A.+?\s+\d+\s",a) #find which position the sequence starts
                        xx = int(xx.end())
                        a = re.sub("\A.+?\s+\d+\s",'',a)
                        a = re.sub("\s\d+\s*\n",'',a)
                        seq = seq + a
                        a = alisc[d][h+1] #alignment annotation
                        a = a[xx:] #remove the extra gaps before annotation begins
                        a = re.sub("\n",'',a)
                        anno = anno + a
                    t=0 #adjusts the index once anno starts changing
                    for a in range(0, len(seq)): #find every gap ("-") in the query sequence and remove corresponding annotation character (" ")
                        if seq[a] == "-":
                            anno = anno[:a-t] +anno[a-t+1:]
                            t=t+1
                    anno = "-"*(que_S[d]-que_E[d]) + anno #add gaps ("-") to beginning of annotation matching not aligned segment length
                    annoF = annoF + anno
                    seq = re.sub("-",'',seq) #remove all gaps from query sequence
                    seq = x.seq[que_E[d]:que_S[d]].lower() + seq #add unaligned sequence
                    seqF = seqF + seq
                seqF = seqF + x.seq[que_E[doms]:].lower() #add sequence end
                annoF = annoF + "-"*(len(x.seq)-que_E[doms]) #add annotation end
                out1= "RS %s,LGIC,%s\n" % (x.id.replace(",","|"),best)
                ho.write(out1) #Write Results: RS line
                ho.write("AH "+','.join(allhits)+"\n") #Write Results: AH line
                ho.write(''.join(dom_ann)) #Write Results: TP lines
                #search against Pfam LGIC related profiles (PfamLGICslib)
                if os.path.isdir(pfamPath):
                    for k in glob.glob(pfamPath+"*.hmm"):
                        tmpHmmOut = "%s.tmp.txt" % infile
                        try:
                            run_hmmsearch = subprocess.run([hmmsearch_path, '-o', tmpHmmOut, k, seq_path])
                            if run_hmmsearch.returncode != 0:
                                raise Exception()
                        except Exception  as ex: #print error message and exit if hmmsearch fails
                            print("hmmsearch could not proceed, ensure HMMER is installed and the path is correct.\n More details:\n")
                            print(ex)
                        hh = open(tmpHmmOut,'r')
                        if "No hits detected" not in hh.read():
                            hh.seek(0)
                            bsProf=  getBestScore(hh)
                            if float(bsProf) > 30: #if the highest scoring domain is above threshold
                                hh.seek(0)
                                (st,en) = getRegions(hh,bsProf) #get the domain start and end positions
                                prf = k.split("/")[-1].split(".")[0] #get profile name from HMM 
                                if st!=0 and en!=0:
                                    out = "TP   %4s   %15s   %5s   %5s   %5s\n" % ("PFAM",prf,st,en,bsProf)
                                    ho.write(out) #Write Results: TP lines
                else:
                    print("PfamLGICslib not found at ", pfamPath)
                seqOut(ho,seqF,"SQ") #Write Results: SQ lines
                seqOut(ho,annoF,"AN") #Write Results: AN lines
                ho.write("//end\n")
            else: #if query protein is not predicted as LGIC
                print("is not predicted as LGIC")
                out1= "RS %s,NOT LGIC,\n" % (x.id.replace(",","|"))
                ho.write(out1)
                seqOut(ho,x.seq,"SQ")
                ho.write("//end\n")

#Helper function to create Brief Results file
def readResultsGetLGICsBrief(filename,outfileBrief):
    counter = 0
    briefResults = ""
    with open(filename,'r') as fh:
        readResultsLines = list(fh.readlines())
        for x in range(0,len(readResultsLines)):
            if re.search("\ARS.+,LGIC", readResultsLines[x]):
                counter += 1
                split2 = readResultsLines[x].split(",")
                briefResults = "%s%-12s\t%12s\n" % (briefResults,split2[0].rstrip(),split2[-1].rstrip())
        lid = filename.split(".")[-2].split("/")[-1]
        with open(outfileBrief, 'w') as brR:
            brR.write("ID %s\n" % lid)
            brR.write("CN %s\n" % counter)
            if counter > 0:
                brR.write(briefResults)


predict(proteinInputFile,outfile) #Main function execution
readResultsGetLGICsBrief(outfile,outfileBrief) #Optional, can be commented out if Brief Results are not needed
