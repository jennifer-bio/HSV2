import numpy as np
import pyvolve
import os
import re
from multiprocessing import Pool
from multiprocessing import Process
import tqdm
import time
import sys
from ete3 import Tree

rates = np.logspace(-3,-9,49)
rateLabels = ["1_0e3", "7_4e4", "5_6e4", "4_2e4","3_1e4", "2_3e4", "1_7e4", "1_3e4", "1_0e4", "7_4e5", "5_6e5", "4_2e5","3_1e5", "2_3e5", "1_7e5", "1_3e5","1_0e5", "7_4e6", "5_6e6", "4_2e6","3_1e6",\
 "2_3e6", "1_7e6", "1_3e6","1_0e6", "7_4e7", "5_6e7", "4_2e7","3_1e7", "2_3e7", "1_7e7", "1_3e7","1_0e7", "7_4e8", "5_6e8",\
  "4_2e8","3_1e8", "2_3e8", "1_7e8", "1_3e8", "1_0e8", "7_4e9", "5_6e9",\
  "4_2e9","3_1e9", "2_3e9", "1_7e9", "1_3e9", "1_0e9"]


iqtree2Path = "TODO"
THREADS = "AUTO"

relativeRates = [0.1612, 0.2472, 18.83, 162.3]
partitionLen = [60448, 42663, 2019, 290]


PARTIONSIZE = 104520

REPS = 200
AGENTS =  "TODO"

basePath =  "TODO"
dateFile = basePath+ "TODO"
regionFile = basePath+ "TODO"
treeFile = basePath+ "TODO"
seqFile = basePath+ "TODO"

#outDirPath = basePath+"05_pyvolve2/"
outDirPath = basePath + "TODO"
print("OUTDIR IS:", outDirPath)


print("SIMULATION METHOD: ", simMethod)

if not os.path.exists(outDirPath):
	os.mkdir(outDirPath)
if not os.path.exists(outDirPath+'EmpiricalClocks/'):
	os.mkdir(outDirPath+'EmpiricalClocks/')
if not os.path.exists(outDirPath+'SimulatedClocks/'):
	os.mkdir(outDirPath+'SimulatedClocks/')
if not os.path.exists(outDirPath+



	 'EmpiricalClocks_Mugration/'):
	os.mkdir(outDirPath+'EmpiricalClocks_Mugration/')
if not os.path.exists(outDirPath+'EmpiricalClocks_Mugration_isAfrica/'):
	os.mkdir(outDirPath+'EmpiricalClocks_Mugration_isAfrica/')

outDataFile = open(outDirPath+"simulatedRates.csv", "w")
outline = "Lable, SimulatedRate, AdjustedRate, Rep, EstimatedRate, rateR_2,  rootNodeDate, rootNodeNumDate \n"
outDataFile.write(outline)
outDataFile.close()

def evolver_seqGen(clockRateLavel, clockRate, rep):

	for cat in range(1,5):
		outFile = simPath+'/simulated_alignment.'+clockRateLavel+'.Rep'+str(rep) + "_cat" + str(cat) +".fasta"
		command_l = "seq-gen -of -l " + str(partitionLen[cat-1]) + " -m GTR -f 0.1607,0.3361,0.3411,0.1621 -r 0.8970,3.4346,0.3170,0.4374,3.2641,1.0 -s " + str(relativeRates[cat-1]*clockRate) + \
			" < " + cleanTreeFileName + " > " + outFile
		os.system(command_l)
		print(command_l)
	simSeq_d = {}

	for cat in range(1,5):
		outFile = simPath+'/simulated_alignment.'+clockRateLavel+'.Rep'+str(rep) + "_cat" + str(cat) +".fasta"
		simGenomes = open(outFile, "r")
		seq = ""
		for line in simGenomes:	
			if ">" in line:
				if seq != "":
					if cat == 1:
						simSeq_d[seqId] = seq
					else:
						simSeq_d[seqId] += seq
				seq = ""
				seqId = line.strip()
			if ">" not in line:
				seq += line.strip()
		if cat == 1:
			simSeq_d[seqId] = seq
		else:
			simSeq_d[seqId] += seq

		simGenomes.close()

	simulated_out = open(simPath+'/simulated_alignment.'+clockRateLavel+'.Rep'+str(rep)+".fasta", "w")
	for seqId in simSeq_d:
		outline = seqId + "\n" + simSeq_d[seqId] + "\n"
		simulated_out.write(outline)
	simulated_out.close()


def evolveAndEstimate(params_l):
	print(" --------------------------- WORKING ON", params_l, " --------------------------- ")
	r = params_l[0:3]
	rep = params_l[3]
	simulatedFasta = simPath+'/simulated_alignment.'+r[1]+'.Rep'+str(rep)+'.fasta'
	if os.path.exists(simulatedFasta) and os.stat(simulatedFasta).st_size != 0:
		print("using previous ", simPath+'/simulated_alignment.'+r[1]+'.Rep'+str(rep)+'.fasta')

	else:
		print (" --------------------------- EVOLVING SEQUENCES --------------------------- ")

		evolver_seqGen(r[1], r[2], str(rep))
		
		time.sleep(5)
		print("COMPETED EVOLVE SIM", r[1]+'.Rep'+str(rep))


	if not os.path.exists(simPath+'/simulated_alignment.'+r[1]+'.Rep'+str(rep)+'.treefile'): 
		print (" --------------------------- BUILD IQTREE --------------------------- ")
		command = iqtree2Path + ' --fast -m GTR+F+R4 -s '+simPath+'/simulated_alignment.'+r[1]+'.Rep'+str(rep)+'.fasta -te ' + treeFile + ' -nt ' + THREADS + \
			 ' --redo --prefix '+simPath+'/simulated_alignment.'+r[1]+'.Rep'+str(rep) +" > temp 2>%1"
		print(command)
		os.system(command)
		print("COMPETED IQTREE", r[1]+'.Rep'+str(rep))
		time.sleep(5)
	else:
		print("using previous", simPath+'/simulated_alignment.'+r[1]+'.Rep'+str(rep)+'.treefile')
	
	ttOutDir = simPath+'/simulated_alignment.'+r[1]+'_TreeTime_Rep'+str(rep)
	


	if not os.path.exists(ttOutDir+"/molecular_clock.txt"):
		print (" --------------------------- ESTIMATED SIMULATED TREETIME --------------------------- ")	
		os.system('treetime --aln '+simPath+'/simulated_alignment.'+r[1]+'.Rep'+str(rep)+'.fasta --tree '+simPath+'/simulated_alignment.'+r[1]+'.Rep'+str(rep)+'.treefile --dates '+dateFile+ \
				' --branch-length-mode input --clock-filter 0 --reroot min_dev --outdir ' + ttOutDir)
		
		print("COMPETED TREETIME", r[1]+'.Rep'+str(rep))
		time.sleep(30)
		
	else:
		print("useing previous", ttOutDir+"/molecular_clock.txt")


	#"Lable, SimulatedRate, AdjustedRate, Rep, EstimatedRate, rateR_2,  rootNodeDate, rootNodeNumDate"
	outline = str(r[1]) + ", " + str(r[0]) + ", " + str(r[2]) + ", " + str(rep)
	#	TODO - extracted needed info for log
	treetimeClockFile = open(ttOutDir+"/molecular_clock.txt")
	for line in treetimeClockFile:
		line_l = line.strip().split()
		if "Root-Tip-Regression" not in line:
			if len(line_l) == 2:
				#print("from molecular_clock:", line_l[1])
				outline += ", " + line_l[1]
	treetimeClockFile.close()


	with open(ttOutDir+"/dates.tsv", "r") as treetimeDatesFile: 
		for i in range(2):
			line2 = next(treetimeDatesFile)
			if "#node" not in line2:
				line2_l = line2.strip().split()
				outline += ", " + str(line2_l[1]) + ", " + str(line2_l[2])
	outline += "\n"
	outDataFile = open(outDirPath+"simulatedRates.csv", "a")
	outDataFile.write(outline)
	outDataFile.close()
	time.sleep(5)


	return 1



for r in zip(rates,rateLabels,rates):
	empPath = outDirPath+'EmpiricalClocks/'+r[1]
	simPath = outDirPath+'SimulatedClocks/'+r[1]
	mugPath = outDirPath+'EmpiricalClocks_Mugration/'+r[1]
	mugPath_Africa = outDirPath+'EmpiricalClocks_Mugration_isAfrica/'+r[1]
	mugPath_groupedAfrica = outDirPath+'EmpiricalClocks_Mugration_groupedAfrica/'+r[1]


	if not os.path.exists(empPath):
		os.mkdir(empPath)
	if not os.path.exists(simPath):
		os.mkdir(simPath)

	print('CURRENTLY RUNNING: '+r[1]+'\n')
	if not os.path.exists(empPath+'/timetree.nexus'):
		print('treetime --keep-root --aln '+seqFile+' --tree '+treeFile+' --dates '+dateFile+' --clock-rate '+str(r[2])+' --branch-length-mode input --clock-filter 0 --outdir '+empPath)
		os.system('treetime --keep-root --aln '+seqFile+' --tree '+treeFile+' --dates '+dateFile+' --clock-rate '+str(r[2])+' --branch-length-mode input --clock-filter 0 --outdir '+empPath)
		time.sleep(20)
	else:
		"Reusing previous Treetime run"
	for line in open(empPath+'/timetree.nexus'):
		if ' Tree tree1' in line: break

	cleanTreeFileName = empPath+"/cleanTree.nwk"
	cleanTree = ''
	skip = 0
	for AB in line[12:].rstrip():
		if AB == '[': skip = 1
		if skip == 0: cleanTree += AB
		if AB == ']': skip = 0
	cleanTree = re.sub("NODE_\d+","",cleanTree)
	# REMOVE NEGATIVE LENGTH BRANCHES
	cleanTree = re.sub("-0.00001","0.00000",cleanTree)

	t=Tree(cleanTree)
	t.resolve_polytomy()
	t.write(outfile=cleanTreeFileName, format=5)

	my_tree = pyvolve.read_tree(tree=cleanTree,scale_tree = r[2])

	os.system('mkdir '+simPath)

	cleanTreeFileName = empPath+"/cleanTree.nwk"

	if not os.path.exists(mugPath+'/annotated_tree.nexus'):
		osCommand = "treetime mugration --tree " + cleanTreeFileName + " --states " + regionFile + " --attribute Region --missing-data Unknown --confidence --outdir " + mugPath
		print(osCommand)
		os.system(osCommand)
	if not os.path.exists(mugPath_Africa+'/annotated_tree.nexus'):
		osCommand = "treetime mugration --tree " + cleanTreeFileName + " --states " + regionFile + " --attribute isAfrica --missing-data Unknown --confidence --outdir " + mugPath_Africa
		print(osCommand)
		os.system(osCommand)
	if not os.path.exists(mugPath_groupedAfrica+'/annotated_tree.nexus'):
		osCommand = "treetime mugration --tree " + cleanTreeFileName + " --states " + regionFile + " --attribute groupedAfrica --missing-data Unknown --confidence --outdir " + mugPath_groupedAfrica
		print(osCommand)
		os.system(osCommand)


	params_l = []

	for rep in range(1,REPS + 1,1):
		temp = list(r)
		temp += [str(rep)]
		params_l.append(temp)

	time.sleep(5)


	print("--- Starting Pool ---")
	pool = Pool(processes=AGENTS)
	pool.map(evolveAndEstimate, params_l)
