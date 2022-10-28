import sys, argparse, os, re

parser = argparse.ArgumentParser(description="")
parser.add_argument("-a", "--alignment", required=True)
parser.add_argument("-t", "--tree", required=True)
parser.add_argument("-s", "--species_tree", required=True)
parser.add_argument("-r", "--root", default=False)
parser.add_argument("-j", "--job_file", default="evaluate_trees.job")

args = parser.parse_args()

def make_job_file(alig, tree, spptree, root, jobfile):
	name = alig[alig.rfind("/")+1:alig.rfind(".")]

	jobfile.write("echo '>PARSIMONY INFORMATIVE SITES' > ./PhyKIT_tests/" + name+"\n")
	jobfile.write("echo '#' >> ./PhyKIT_tests/" + name+"\n" )
	jobfile.write("phykit parsimony_informative_sites " + alig + " >> ./PhyKIT_tests/" + name+"\n\n")
	
	jobfile.write("echo '>TREENESS' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("echo '#This measures the ration between phylogenetic signal and noise. Higher treeness is desirable'  >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("phykit treeness " + tree + " >> ./PhyKIT_tests/" + name+"\n\n")
	
	jobfile.write("echo '>RELATIVE COMPOSITION VARIABILITY' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("echo '#This measures the variability in the compositoin of the alignment. Lower values are desirable.' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("phykit relative_composition_variability " + alig + " >> ./PhyKIT_tests/" + name+"\n\n")
	
	jobfile.write("echo '>TREENESS/RCV' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("echo '#This combines treeness and RCV. Higer values are desirable.' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("phykit treeness_over_rcv -a " + alig + " -t " + tree + " >> ./PhyKIT_tests/" + name+"\n\n")
	
	jobfile.write("echo '>BIPARTITION SUPPORT STATS' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("echo '#Higher values are desirable.'  >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("phykit bipartition_support_stats " + tree + " >> ./PhyKIT_tests/" + name+"\n\n")
	
	if root != False:
		jobfile.write("echo '>DEGREE OF VIOLATION OVER THE MOLECULAR CLOCK' >> ./PhyKIT_tests/" + name+"\n")
		jobfile.write("echo '#Measures how much the tree deviates from a molecular clock model. Lower values are desirable.' >> ./PhyKIT_tests/" + name+"\n")
		jobfile.write("phykit dvmc  -t " + tree + " -r " + root +" >> ./PhyKIT_tests/" + name+"\n\n")
	
	jobfile.write("echo '>SATURATION' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("echo '#This measure saturation in the tree and alignment. Varies from 0 to 1, and lower values are desirable' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("phykit saturation -a " + alig + " -t " + tree + " >> ./PhyKIT_tests/" + name+"\n\n")
	
	jobfile.write("echo '>ROBINSON-FOULDS DISTANCE TO SPECIES TREE' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("echo '#' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("phykit robinson_foulds_distance " + tree + " " + spptree + " >> ./PhyKIT_tests/" + name+"\n\n")
	
	jobfile.write("echo '>EVOLUTIONARY RATE' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("echo '#' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("phykit evolutionary_rate " + tree + " >> ./PhyKIT_tests/" + name+"\n\n")
	
	jobfile.write("echo '>PAIRWISE IDENTITY' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("echo '#This is used as an approximation to the evolutionary rate based on the alignment' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("phykit pairwise_identity " + alig + " >> ./PhyKIT_tests/" + name+"\n\n")
	
	jobfile.write("echo '>INTERNAL BRANCH STATS' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("echo '#' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("phykit internal_branch_stats " + tree + " >> ./PhyKIT_tests/" + name+"\n\n")
	
	jobfile.write("echo '>TERMINAL BRANCH STATS' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("echo '#' >> ./PhyKIT_tests/" + name+"\n")
	jobfile.write("phykit terminal_branch_stats " + tree + " >> ./PhyKIT_tests/" + name+"\n\n")
	
if os.path.isfile(args.alignment) == True and os.path.isfile(args.tree) == True:
	jobfile = open(args.job_file, "w")
	jobfile.write("mkdir ./PhyKIT_tests\n")
	make_job_file(args.alignment, args.tree, args.species_tree, args.root, jobfile)

elif os.path.isdir(args.alignment) == True and os.path.isdir(args.tree) == True:
	jobfile = open(args.job_file, "w")
	jobfile.write("mkdir ./PhyKIT_tests\n")
	ali2tree = {}
	for i in os.listdir(args.alignment):
		name = i.split(".")[0]
		for e in os.listdir(args.tree):
			if e.find(name) > -1:
				ali2tree[os.path.abspath(args.alignment+i)] = os.path.abspath(args.tree+e)
				break
	for entry in ali2tree:
		make_job_file(entry, ali2tree[entry], args.species_tree, args.root, jobfile)

else:
	print ("Alignment and Tree need to concordate. One of them is a file and one is a directory")
