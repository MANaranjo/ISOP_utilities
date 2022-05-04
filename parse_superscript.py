import sys, os
import shutil
from Bio import SeqIO

dirinput = sys.argv[1]

def pick_longest(fastafile):
	fa_dict = SeqIO.parse(fastafile, "fasta")
	longest = ""
	for i in fa_dict:
		if len(i.seq) > len(longest):
			longest = str(i.seq)
	return(longest)

def parse_superscript(dirinput):
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
					a = f.split(".txt.")
					samplename = d 
					seqname = a[1].replace(samplename+"_", "")
					myfasta = outdir + seqname.split(".")[0] + ".fasta"
					
					if os.path.exists(myfasta) == True:
						outfile = open(myfasta, "a")
						outfile.write(">"+samplename+"\n")
						outfile.write(longest)
						outfile.write("/n")
						outfile.close()
					else:
						outfile = open(myfasta, "w")
						outfile.write(">"+samplename+"\n")
						outfile.write(longest)
						outfile.write("/n")
						outfile.close()
						continue
					
parse_superscript(dirinput)
