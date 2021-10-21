import subprocess
import os

families = ["5HT3", "ASICs", "ENaCs", "GABAA", "iGluR", "Glycine", "IP3R", "nAchR", "P2XR", "ZACs"]
paths = ["LGIC_cdhit/", "LGIC_clustal/", "LGIC_hmmer/"]
fasta = "LGIC_fasta/"

cdhit = "cdhit"
clustal = "clustalo"
hmmer = "/usr/local/bin/"

for path in paths:
    if not os.path.isdir(path): #check if directory for results does not exist
        os.makedirs(path) #create directory

for fam in families:
    print(fam)
    
    run_cdhit = subprocess.run([cdhit, "-i", fasta+fam+".fasta", "-o", paths[0]+"cd40"+fam, "-c", "0.4", "-n", "2"])
    print("The exit code for CDHIT was: %d" % run_cdhit.returncode)
    
    run_clustal = subprocess.run([clustal, "-i", paths[0]+"cd40"+fam, "-o", paths[1]+fam+".clu", "--outfmt", "clu"])
    print("The exit code for ClustalO was: %d" % run_clustal.returncode)
    
    run_hmmbuild = subprocess.run([hmmer+"hmmbuild", "-n", fam, "-o", paths[2]+fam+".out", paths[2]+fam+".hmm", paths[1]+fam+".clu"])
    print("The exit code for hmmbuild was: %d" % run_hmmbuild.returncode)

with open('allLGIChmms', 'w') as outfile:
    for fam in families:
        with open(paths[2]+fam+".hmm") as infile:
            outfile.write(infile.read())

run_hmmpress = subprocess.run([hmmer+"hmmpress", 'allLGIChmms'])
print("The exit code for hmmpress was: %d" % run_hmmpress.returncode)
