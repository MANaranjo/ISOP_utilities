import sys, argparse
from Bio import AlignIO, SeqIO

'''
This script wil take a multifasta file as input and will generate a multifasta file containing those sequences whose headers satisfy the conditions set by the user.
'''

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', required=True, help='Fasta file to parse. Required.')
parser.add_argument('-d', '--destination', default=False, help="File to write output. If ommited, it will write the results to the standard output")
parser.add_argument("-i", "--include", nargs="+", default=[], help="All headers that include any of these strings will be included")
parser.add_argument("-e", "--exclude", nargs="+", default=[], help="All headers that include any of these strings will be included")
parser.add_argument("-I", "--include_strict", nargs="+", default=[], help="All headers exactly equal to any of these strings will be included")
parser.add_argument("-E", "--exclude_strict", nargs="+", default=[], help="All headers that are exactly equal to any of these strings will be excluded")
parser.add_argument("-n", "--include_file", default=False, help="File including strings to search, one per line. All headers including any of these strings will be included")
parser.add_argument("-x", "--exclude_file", default=False, help="File including strings to search, one per line. All headers including any of these strings will be excluded")
parser.add_argument("-N", "--include_file_strict", default=False, help="File including strings to search, one per line. All headers including exactly any of these strings will be included")
parser.add_argument("-X", "--exclude_file_strict", default=False, help="File including strings to exclude, one per line. All headers including exactly any of these strings will be excluded")
parser.add_argument("-s", "--slow_mode", default=False, action="store_true", help="Uses a slower mode to parse the fasta file that does not depend on Biopython. Use it if you are having problems with the way Biopython parses fasta headers. Be careful, if your file is really big (In the order of Gb) it might run out of memory")

args = parser.parse_args()
include_all = False
if len(args.include)+len(args.include_strict) == 0 and args.include_file == False and args.include_file_strict == False:
	include_all = True

def fasta_to_dict(myfasta, mode):
	ogfasta = {}
	if mode == False:
		ogfasta = SeqIO.index(myfasta, "fasta")
	else:
		name = False
		for line in open(myfasta):
			if line[0] == ">":
				name = line[1:-1]
				ogfasta[name] = ""
			elif len(line) > 1:
				ogfasta[name] = ogfasta[name]+line[:-1]
			else:
				continue
	return (ogfasta)

ogfasta = fasta_to_dict(args.fasta, args.slow_mode)

includelist = []

if len(args.include) > 0:
	for i in args.include:
		includelist.append(i)
if args.include_file != False:
	f = open(args.include_file)
	for line in f:
		includelist.append(line[:-1])
		
excludelist = []
if len(args.exclude) > 0:
	for i in args.exclude:
		excludelist.append(i)
if args.exclude_file != False:
	f = open(args.exclude_file)
	for line in f:
		excludelist.append(line[:-1])

strict_include = []
if len(args.include_strict) > 0:
	for i in args.include_strict:
		strict_include.append(i)
if args.include_file_strict != False:
	f = open(args.include_file_strict)
	for line in f:
		strict_include.append(line[:-1].split(" ")[0])

strict_exclude = []
if len(args.exclude_strict) > 0:
	for i in args.exclude_include:
		strict_exclude.append(i)
if args.exclude_file_strict != False:
	f = open(args.exclude_file_strict)
	for line in f:
		strict_exclude.append(line[:-1].split(" ")[0])

to_output = []
not_output = []

for seq in ogfasta:
	if include_all == True:
		to_output.append(seq)
	else:
		for name in includelist:
			if args.slow_mode == False:
				if ogfasta[seq].description.find(name) > -1:
					to_output.append(seq)
			else:
				if seq.find(name) > -1:
					to_output.append(seq)
		for iname in strict_include:
			if iname not in to_output and seq.id == iname:
				to_output.append(seq)
				
if args.destination != False:
	destin = open(args.destination, "w")

for seqid in to_output:
	for x in excludelist:
		if ogfasta[seqid].description.find(x) > -1:
			not_output.append(seqid)
	else:
		for xx in strict_exclude:
			if seqid == xx:
				not_output.append(seqid)
for seqid in to_output:
	if seqid in not_output: continue
	else:
		if args.destination == False:
			print(">"+seqid)
			if args.slow_mode == False:
				print(ogfasta[seqid].seq+"\n")
			else:
				print(ogfasta[seqid]+"\n")
			
		else:
			if args.slow_mode == False:
				destin.write(">"+ogfasta[seqid].description+"\n")			
				destin.write(str(ogfasta[seqid].seq))
			else:
				destin.write(">"+seqid+"\n")
				destin.write(ogfasta[seqid])
			destin.write("\n\n")

