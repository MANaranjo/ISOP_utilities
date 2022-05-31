import sys, os, random
from Bio import AlignIO
from Bio.Nexus import Nexus

fastalist = sys.argv[1:]

outpref = "concatenated"	

concatenated_align_list = []
for element in fastalist:
	newfasta = open(element+".new.fa", "w")
	for line in open(element):
		if line[0] == '>':
			newfasta.write(line[:line.find('|')]+"\n")
		else:
			newfasta.write(line.replace("*", "-").replace("?","-"))
	newfasta.close()
	AlignIO.convert(element+".new.fa", "fasta", element+".nex", "nexus", 'protein')
	concatenated_align_list.append(element+".nex")

data = [(fname, Nexus.Nexus(fname)) for fname in concatenated_align_list]

combined = Nexus.combine(data)
combined.write_nexus_data(filename=open(outpref+'.nex', 'w'))
combined.export_fasta(filename=outpref+'.fa')

os.system("clipkit "+outpref+".fa -o "+outpref+"_clipkit.fa \n")
os.system("sed 's/?/-/g' "+outpref+"_clipkit.fa > tmp ; mv tmp "+outpref+"_clipkit.fa")

for i in fastalist:
	if os.path.isfile(i+".new.fa") == True:
		os.remove(i+".new.fa")
		os.remove(i+".nex")


