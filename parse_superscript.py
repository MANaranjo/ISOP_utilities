import sys, os, argparse
import shutil
from Bio import SeqIO

parser = argparse.ArgumentParser(description="")
parser.add_argument("-d", "--directory", default="./", help="Directory output from superscript.py")
parser.add_argument("-m", "--missing", default=0, help="Maximum number of missing species per sequence that will be accepted")
parser.add_argument("-c", "--clipkit", default=False, action="store_true", help="When set, the output job file will contain commands to launch clipkit on individual ")
parser.add_argument("-q", "--query_file", default = "something.fasta", help="Query file used to search the sequences. This is only used to process the names of the files by using the file extension")


args = parser.parse_args()

dirinput = os.path.abspath(args.directory)

def estimate_total (my_dir):
	mylist = []
	for i in os.listdir(my_dir):
		if os.path.isdir(i) == True and i.find("parsed_seqs") == -1 and i.find("selected_seqs"):
			mylist.append(i)
	return(len(mylist))

def pick_longest(fastafile):
	fa_dict = SeqIO.parse(fastafile, "fasta")
	longest = ""
	for i in fa_dict:
		if len(i.seq) > len(longest):
			longest = str(i.seq)
	return(longest)

def parse_superscript(dirinput):
	s = args.query[args.query.find("."):]+"."
	if dirinput[-1] != "/": 
		dirinput = dirinput + "/"
	outdir = dirinput+"parsed_seqs/"
	if os.path.exists(outdir) == True:
		shutil.rmtree(outdir)
	os.makedirs(outdir)
	for d in os.listdir(dirinput):
		if os.path.isdir(d) == True:
			subdir = dirinput+d+"/"
			for f in os.listdir(subdir):
				if f.find("filtered_contigs.fa") > -1:
					longest = pick_longest(subdir+f)
					a = f.split(s)
					samplename = d 
					seqname = a[1].replace(samplename+"_", "")
					myfasta = outdir + seqname.split(".")[0] + ".fasta"
					
					if os.path.exists(myfasta) == True:
						outfile = open(myfasta, "a")
						outfile.write(">"+samplename+"\n")
						outfile.write(longest)
						outfile.write("\n")
						outfile.close()
					else:
						outfile = open(myfasta, "w")
						outfile.write(">"+samplename+"\n")
						outfile.write(longest)
						outfile.write("\n")
						outfile.close()
						continue
	return(outdir)
					
def extract_seqs(dirinput, outdir, cutoff):
	outdir2 = dirinput+"/selected_seqs/"
	if os.path.exists(outdir2) == True:
		shutil.rmtree(outdir2)
	os.makedirs(outdir2)
	for f in os.listdir(outdir):
		fasta = SeqIO.parse(outdir+f, "fasta")
		faslist = []
		for i in fasta:
			faslist.append(i.id)
		if len(faslist) >= cutoff:
			shutil.copyfile(outdir+f, outdir2+f)
	return(outdir2)
			
def launch_mafft(outdir2, outfile):
	for f in os.listdir(outdir2):
		newname = f[:f.rfind(".")]+".mafft"
		outfile.write("mafft --auto "+outdir2+f+" > "+outdir2+newname+"\n")

def launch_clipkit(outdir2, outfile):
	for f in os.listdir(outdir2):
		newname = f[:f.rfind(".")]+"_clipkit.fasta"
		mafftname = f[:f.rfind(".")]+".mafft"
		outfile.write("clipkit "+outdir2+mafftname+" -o "+outdir2+newname+"\n")

cutoff = estimate_total(dirinput) - int(args.missing)
outdir = parse_superscript(dirinput)
outdir2 = extract_seqs(dirinput, outdir, cutoff)
outfile = open(dirinput+"/alignment_job.txt", "w")
launch_mafft(outdir2, outfile)
if args.clipkit != False:
	launch_clipkit(outdir2, outfile)
outfile.close()
