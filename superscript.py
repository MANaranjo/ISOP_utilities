from Bio import SeqIO
import sys, argparse, os, re
import numpy
import gzip, bz2, tarfile
import shutil
from Bio import SeqIO

parser = argparse.ArgumentParser(description="")
parser.add_argument("-T", "--trimming", default=False, choices=['fastp','trimmomatic'])
parser.add_argument('-l', '--libraries', required=True, nargs='+', help="Fastq libraries to use. Required.")
parser.add_argument("-x", "--trimmomatic_commands", default="")
parser.add_argument("-c", "--conda", default="./", help="bin directory for your conda installation.")
parser.add_argument("-d", "--directory", default="./superscript_output", help="Directory to send all the output files")
parser.add_argument("-t", "--table", default=False, help="File containing a conversion table for the libraries. The format must be a csv file separated by ';' with the first column containing the names of the libraries and the second containing the name of the species. Ensure the names of the species are unique, adding some identifier if needed.")
parser.add_argument("-a", "--aTRAM", default=False, help="Directory containing aTRAM")
parser.add_argument("--atram_program", default="spades", choices=["trinity", "spades"])
parser.add_argument("-p", "--HybPiper", default=False)
parser.add_argument("-b", "--bait", default=False, nargs="+")
parser.add_argument("-z", "--compression", default="", choices=["gzip", "bzip"])

args = parser.parse_args()
here = os.path.abspath("./")

if args.conda[-1] != "/":
	args.conda = args.conda + "/"

directory = os.path.abspath(args.directory)
if directory[-1] != "/":
	directory = directory + "/"

aTRAM = args.aTRAM
if args.aTRAM != False:
	aTRAM = os.path.abspath(aTRAM)
	if aTRAM[-1] != "/":
		aTRAM = aTRAM + "/"

HybPiper = args.HybPiper
if args.HybPiper != False:
	HybPiper = os.path.abspath(HybPiper)
	if HybPiper[-1] != "/":
		HybPiper = HybPiper + "/"
bait = args.bait
if aTRAM != False or HybPiper != False:
	if bait == False:
		raise NameError("You need a bait file in order to use aTRAM or HybPiper!")
	else:
		for n in args.bait:
			n = os.path.abspath(n)
			if n[-1] != "/":
				n = n + "/"

def get_fasta_type(fasta_list):
	amino_dict = dict()
	for fasta in fastalist:
		amino = False
		my_fasta = SeqIO.parse(fasta, "fasta")
		for my_seq in my_fasta:
			if amino == True:
				break
			else:
				if len(my_seq.seq) != my_seq.seq.count("A")+my_seq.seq.count("a")+my_seq.seq.count("G")+my_seq.seq.count("g")+my_seq.seq.count("T")+my_seq.seq.count("t")+my_seq.seq.count("C")+my_seq.seq.count("c"):
					amino = True
		amino_dict[fasta] = amino
	return(amino_dict)

########################################################################

def remove_false_files(filelist): #We remove empty files that otherwise would make everything crash
	clean_list = []
	for i in filelist: 
		if os.path.getsize(i) < 10000: #I probably need a better way to judge it
			continue
		else:
			clean_list.append(i)
	return clean_list

def get_mean_read_len(fastqfile, sample_size, compressed_dict):
	mean_read_dict = {}
	sampling = []
	if (compressed_dict[fastqfile] != "no-compression" and fastqfile[:fastqfile.rfind(".")+1][-3:] == "fa") or (compressed_dict[fastqfile] != "no-compression" and  fastqfile[:fastqfile.rfind(".")+1][-6:] == ".fasta"):
		fastafile = SeqIO.parse(fastqfile, "fasta")
		for i in fastafile:
			sampling.append(len(i))
			if len(sampling) >= sample_size:
				break
	elif fastqfile[-3:] == "fa" or fastqfile[-6:] == ".fasta" :
			fastafile = SeqIO.parse(fastqfile, "fasta")
			for i in fastafile:
				sampling.append(len(i))
				if len(sampling) >= sample_size:
					break
	else:
		switch = False
		if compressed_dict[fastqfile] == "gzip":
			open_fastqfile = gzip.open(fastqfile, 'r')
		elif compressed_dict[fastqfile] == "bz2":
			open_fastqfile = bz2.BZ2file(fastqfile, "r")
		elif compressed_dict[fastqfile] == "tar":
			open_fastqfile = tarfile.open(fastqfile, "r")
		else:
			open_fastqfile = open(fastqfile)
		for line in open_fastqfile:
			if switch == True:
				sampling.append(len(line)-1)
				switch = False
				if len(sampling) >= sample_size:
					break
			if type(line) is str:
				if line[0] == "+":
					switch = True
			elif type(line) == bytes:
				try:
					bline = line.decode('cp437')
				except UnicodeDecodeError :
					continue
				else:
					if bline.find("+") < 5:
						switch = True
			else: continue
	return [int(numpy.mean(sampling)), numpy.std(sampling)]

def compression_parse(fastq):
	compressed_dict = {}
	for i in fastq:
		if i[i.rfind("gz"):] == "gz" or i[i.rfind("gzip"):] == "gzip":
			compressed_dict[i] = "gzip"
		elif i[i.rfind("bzip2"):] == "bz2" or i[i.rfind("bzip2"):] == "bzip2" > -1:
			compressed_dict[i] = "bz2"
		elif i[i.rfind("tar"):] == "tar":
			compressed_dict[i] = "tar"
		else:
			compressed_dict[i] = "no-compression"
	return compressed_dict

def do_line(line, switch, phred64dict, counter):
	breakswitch = False
	if switch == False:
		if line[0] == "+":
			switch = True
			counter = counter + 1
	else:
		if switch == True:
			if line.find("Z") > -1:
				phred64dict[element] = "64"
				breakswitch = True
			else:
				switch = False
	return (switch, counter, breakswitch)

def phred_parse (fastqlist, sample_size):
	phred64dict = {}
	counter = 0
	for element in fastqlist:
		switch = False
		if element.split('.')[-1] == 'gz' or element.split('.')[-1] == 'gzip':
			for line in gzip.open(element, 'r'):
				bline = line.decode()[2:-3]
				if counter > sample_size: break
				elif len(bline) < 3: continue
				elif line[0] == "@": continue
				else:
					switch, counter, breakswitch = do_line(bline, switch, phred64dict, counter)
					if breakswitch == True:
						break
		elif element.split('.')[-1] == 'bz2' or element.split('.')[-1] == 'bzip2':
			for line in bz2.BZ2file(element, "r"):
				bline = line.decode()[2:-3]
				if counter > sample_size: break
				elif len(bline) < 3: continue
				elif line[0] == "@": continue
				else:
					switch, counter, breakswitch = do_line(line, switch, phred64dict, counter)
					if breakswitch == True:
						break
		elif element.split('.')[-1] == 'tar':
			for line in tarfile.open(element, "r"):
				bline = line.decode()[2:-3]
				if counter > sample_size: break
				elif len(bline) < 3: continue
				elif line[0] == "@": continue
				else:
					switch, counter, breakswitch = do_line(line, switch, phred64dict, counter)
					if breakswitch == True:
						break
		else:
			for line in open(element):
				if counter > sample_size: break
				elif len(line) < 3: continue
				elif line[0] == "@": continue
				else:
					switch, counter, breakswitch = do_line(line, switch, phred64dict, counter)
					if breakswitch == True:
						break

	for element in fastqlist:
		if element not in phred64dict:
			phred64dict[element] = "33"
	return phred64dict

def hypo_dict_parse(fastqlist):
	hypo_dict = {}
	for element in fastqlist:
		if element in hypo_dict:
			continue
		for m in re.finditer('2', element):
			hypothetical = element[:m.start()]+"1"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				break
		for m in re.finditer('R', element):
			hypothetical = element[:m.start()]+"F"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				break
		for m in re.finditer('rev', element):
			hypothetical = element[:m.start()]+"fwd"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				break
		for m in re.finditer('rev', element):
			hypothetical = element[:m.start()]+"fw"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				break
		for m in re.finditer('R2', element):
			hypothetical = element[:m.start()]+"R1"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				break
	return hypo_dict

def format_parse(fastq):
	format_dict = {}
	for i in fastq:
		if i[-3:] == "fa" or i[-6:] == ".fasta" :
			format_dict[i] = "fasta"
		else:
			format_dict[i] = "fastq"
	return format_dict

def type_parse(fastq, hypo_dict, mean_read_dict):
	type_dict = {}
	library_dict = {}
	library_size_dict = {}
	for i in fastq:
		if i in type_dict: continue
		if i in hypo_dict:
			type_dict[i] = [1, hypo_dict[i]]
			type_dict[hypo_dict[i]] = [2, i]
		elif int(mean_read_dict[i][0]) > 1500:
			type_dict[i] = ["pb", "no_partner"]
		else: type_dict[i] = ["s", "no_partner"]
	for i in fastq:
		if type_dict[i] == 1 or type_dict[i] == 2:
			library_size_dict[i] = (os.stat(i).st_size)*2
		else:
			library_size_dict[i] = (os.stat(i).st_size)
	return type_dict, library_size_dict

# this function creates the prepared_libraries.txt file
def preparation(initial_fastq, sample_size, directory):
	os.mkdir(directory)
	library_report = directory+"libraries.txt"
	fastq = remove_false_files(initial_fastq)
	mean_read_dict = {}
	compressed_dict = compression_parse(fastq)
	toreturn = []
	for i in fastq:
		mean_read_dict[i] = get_mean_read_len(i, sample_size, compressed_dict)

	hypo_dict = hypo_dict_parse(fastq)
	type_dict, library_size_dict = type_parse(fastq, hypo_dict, mean_read_dict)


	phred64dict = phred_parse(fastq, sample_size)
	format_dict = format_parse(fastq)

	report = open(library_report, 'w')
	for i in fastq:
		print (i + "\t" + str(mean_read_dict[i][0]) + "\t" + str(mean_read_dict[i][1]) + "\t" + str(library_size_dict[i]) + "\t" + str(phred64dict[i]) + "\t" + str(type_dict[i][0]) + "\t" + str(type_dict[i][1]) + "\t" + format_dict[i] + "\t" + compressed_dict[i]+"\n")
		report.write(i + "\t" + str(mean_read_dict[i][0]) + "\t" + str(mean_read_dict[i][1]) + "\t" + str(library_size_dict[i]) + "\t" + str(phred64dict[i]) + "\t" + str(type_dict[i][0]) + "\t" + str(type_dict[i][1]) + "\t" + format_dict[i] + "\t" + compressed_dict[i]+"\n")
		toreturn.append((i + "\t" + str(mean_read_dict[i][0]) + "\t" + str(mean_read_dict[i][1]) + "\t" + str(library_size_dict[i]) + "\t" + str(phred64dict[i]) + "\t" + str(type_dict[i][0]) + "\t" + str(type_dict[i][1]) + "\t" + format_dict[i] + "\t" + compressed_dict[i]+"\n"))
	report.close()
	return(toreturn, library_report)

def trimmomatic (paired_list, single_list, conversion, directory, path):
	paired_list = []
	single_list = []
	pacbio_list = []
	backstring = ''
	trimmo_exec = ''
	output_list = []
	if len(commands) == 0:
		commands = 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

	for i in os.listdir(path):
		if i[-4:] == '.jar' and i.find('trimmomatic') > -1:
			trimmo_exec = path + i

	if "All_adapters.fa" not in os.listdir(path+"adapters"):
		allada = open(path+"adapters/All_adapters.fa", "w")
		for i in os.listdir(path+"adapters"):
			for line in i:
				allada.write(line[:-1])
			allada.write("\n")
		allada.close()

	for i in open(library_file):
		chunk = i.split()
		if chunk[5] == "1":
			paired_list.append([chunk[0], chunk[6], chunk[4]])
		elif chunk[5] == "2": continue
		elif chunk[5] == "s":
			single_list.append(chunk[0])
		elif chunk[5] == "pb":
			pacbio_list.append(chunk[0])
		else: continue
	
	output_string = ''
	to_remove = []

	for i in paired_list:
		if i[0][0] == '.': i[0] = i[0][1:]
		if i[1][0] == '.': i[1] = i[1][1:]
		output_string = output_string + "java -jar " + trimmo_exec + " PE -phred" + str(i[2]) + " " + i[0] + " " + i[1] + " "+output+"parsed_paired_"+ \
		i[0][i[0].rfind("/")+1:] + " " + output+"parsed_unpaired_"+i[0][i[0].rfind("/")+1:] + " " +output+"parsed_paired_"+i[1][i[1].rfind("/")+1:]\
		 + " " + output+"parsed_unpaired_"+i[1][i[1].rfind("/")+1:] + " " +"ILLUMINACLIP:" + path+"adapters/All_adapters.fa" + ":2:30:10 "+commands+"\n"
		if remove_originals == True:
			to_remove.append(i[0])
			to_remove.append(i[1])

	for i in single_list:
		if i[0] == '.': i = i[1:]
		output_string = output_string + "java -jar " + trimmo_exec + " SE -phred" + str(i[2]) + " " + i[0] + " " + i[1] + " " +output+i[0][:i[0].rfind("/")+1]+"parsed_"+\
		" ILLUMINACLIP:" + path+"adapters/All_adapters.fa" + ":2:30:10 "+commands+"\n"
		if remove_originals == True:
			to_remove.append(i[0])

	for i in to_remove:
		output_string = output_string + "rm " + i + "\n"
	return(output_string, output_list)

def fastp (paired_list, single_list, conversion, directory, path):
	output_string = '\necho "Configuring commands for fastp..."'+'\n\n'
	output_list = []
	for i in paired_list:
		if i[0].find('./') == 0 : i[0] = i[0][2:]
		if i[1].find('./') == 0 : i[1] = i[1][2:]
		print (i[0], i[1])
		destiny_list = i[0][i[0].rfind("/")+1:].replace(".", " ").split()
		destiny = destiny_list[0]
		if conversion != False and destiny in conversion:
			destiny = conversion[destiny]
		o1, o2 = i[0][i[0].rfind("/")+1:], i[1][i[1].rfind("/")+1:]
		o = [directory+destiny, o1[:o1.rfind(".")]+"parsed_paired.fq", o2[:o2.rfind(".")]+"parsed_paired.fq"]		
		os.mkdir(directory+destiny)
		output_string = output_string + path + "bin/fastp -i " + os.path.abspath(i[0]) + " -I " + os.path.abspath(i[1]) + " -o " + o[0]+"/"+o[1].split("/")[-1]+".gz" + " -O " + o[0]+"/"+o[2].split("/")[-1]+".gz" + " --detect_adapter_for_pe\n"
		output_list.append(o)

	for i in single_list:
		if i[0] == '.': i = i[1:]
		if conversion != False and destiny in conversion:
			destiny = conversion[destiny]
		output_string = output_string + path + "bin/fastp -i " + os.path.abspath(i[0]) + " -I " + os.path.abspath(i[1]) + " -o "+directory+destiny+"/"+destiny+"_parsed_single"+"\n"
	return(output_string, output_list)

def no_trimming(paired_list, single_list, conversion, directory, path):
	output_string = ''
	output_list = []
	suf = ""
	if args.compression == "gzip" and args.trimming == False:
		comp = " --gzip "
		suf = ".gz"
	for i in paired_list:
		if i[0].find('./') == 0 : i[0] = i[0][2:]
		if i[1].find('./') == 0 : i[1] = i[1][2:]
		destiny_list = i[0][i[0].rfind("/")+1:].replace(".", " ").split()
		destiny = destiny_list[0]
		if conversion != False and destiny in conversion:
			destiny = conversion[destiny]
		o1, o2 = i[0][i[0].rfind("/")+1:], i[1][i[1].rfind("/")+1:]
		o = [directory+destiny, o1[:o1.rfind(".")]+".fq"+suf, o2[:o2.rfind(".")]+".fq"+suf]
		os.mkdir(directory+destiny)
		output_string = output_string + "ln -s " + os.path.abspath(i[1]) + " " + o[0]+"/"+o[1] + "\nln -s " + os.path.abspath(i[1]) + " " + o[0]+"/"+o[2][o[2].rfind("/")+1:] + "\n"
		output_list.append(o)

	for i in single_list:
		if i[0] == '.': i = i[1:]
		if conversion != False and destiny in conversion:
			destiny = conversion[destiny]
		output_string = output_string + "ln -s " + os.path.abspath(i[0]) + " " + o[0]+"/"+o[1] + "\n"
	return(output_string, output_list)

def atram(output_list, bait, aTRAM):
	output_string = '\necho "Configuring commands for aTRAM..."'+'\n\n'
	suf = ""
	if args.trimming == "fastp":
		comp = " --gzip "
		suf = ".gz"
	else: comp = ""
	for o in output_list:
		ministring1 = aTRAM+"atram_preprocessor.py -b "+o[0]+o[0][o[0].rfind("/"):]+" --end-1 "+o[0]+"/"+o[1]+suf+" --end-2 "+o[0]+"/"+o[2]+suf+comp+" -t "+ o[0]+"/tmp\n"
		ministring2 = ""
		for n in bait:
			N = os.path.abspath(n)
			if os.path.isdir(o[0]+"/tmp") == False:
				os.mkdir(o[0]+"/tmp")
			ministring2 = ministring2 + aTRAM+"atram.py -b "+o[0]+o[0][o[0].rfind("/"):]+" -Q "+N+" -a " + args.atram_program + " -o "+o[1][:o[1].rfind(".")]+"_"+n[n.rfind("/")+1:]+" -t "+ o[0]+"/tmp\n"
		output_string = output_string + ministring1 + ministring2
	return (output_string)
	
def hybpiper(o, bait, HybPiper):
	pass
########################################################################

def arrange_files(prep, table, path, directory, aTRAM, HybPiper, bait):
	outfile = open(directory+"super_script.job", "w")
	paired_list = []
	single_list = []
	pacbio_list = []
	backstring = ''
	for i in prep:
		chunk = i.split()
		if chunk[5] == "1":
			paired_list.append([chunk[0], chunk[6], chunk[4]])
		elif chunk[5] == "2": continue
		elif chunk[5] == "s":
			single_list.append(chunk[0])
		elif chunk[5] == "pb":
			pacbio_list.append(chunk[0])
		else: continue
	output_string = ''
	if args.trimming == False:
		ministring, output_list = no_trimming(paired_list, single_list, args.table, directory, path)
		output_string = output_string + ministring
	if args.trimming == "fastp":
		ministring, output_list = fastp(paired_list, single_list, args.table, directory, path)
		output_string = output_string + ministring
	if args.trimming == "trimmomatic":
		ministring, output_list = output_string + trimmomatic(paired_list, single_list, args.table, directory, path)
		output_string = output_string + ministring
	if aTRAM != False:
		atram_string = atram(output_list, bait, aTRAM)
		output_string = output_string + atram_string
	outfile.write(output_string)	
		
def parse_conversion_table (table):
	ctable = open(table)
	cdict = {}
	for line in ctable:
		chunk = line[:-1].split(";")
		ctable[chunk[0]] == ctable[chunk[1]]
	return (cdict)

prep, output_report = preparation(args.libraries, 100, directory)
arrange_files(prep, args.table, args.conda, directory, aTRAM, HybPiper, bait)
