import sys
import string

recombination_file_name = "filtering/data/combinedCatOnlyBestWithPVals.txt"
sampleinfo_file_names = ["filtering/data/sampleInfo.txt"]
all_fasta_file_name = "filtering/fastas/extractedSeqs.fa" 
reference_file_name = "filtering/fastas/reference.fa"

print_all = (sys.argv[1] == '-a')

F = open(recombination_file_name, 'r')
lines = F.readlines()


labels = lines[0].split('\t')
labels[11] = labels[11].strip()

del lines[0]

parsimony_change = []

for i in range(len(lines)):
	lines[i] = lines[i].split('\t')
	
	parsimony_change.append(int(lines[i][10]) - int(lines[i][11]))



lines_sorted_by_parsimony_change = [x for _, x in sorted(zip(parsimony_change, lines),reverse=True)]


how_many_to_see = 0 if (len(sys.argv) == 1 or print_all) else int(sys.argv[1])



for i in range(how_many_to_see):

	print("Example", i)

	for j in range(len(labels)):

		print("\t", labels[j], lines_sorted_by_parsimony_change[i][j])


how_many_to_see = (len(lines_sorted_by_parsimony_change)) if print_all else how_many_to_see


make_fasta = (len(sys.argv) > 2) or print_all

if not make_fasta:
	sys.exit()


nodes_to_get = []
nodes_by_example = [[] for x in range(how_many_to_see)]

for i in range(how_many_to_see):
	if not (lines_sorted_by_parsimony_change[i][0] in nodes_to_get):
		nodes_to_get.append(lines_sorted_by_parsimony_change[i][0])
	if not (lines_sorted_by_parsimony_change[i][3] in nodes_to_get):
		nodes_to_get.append(lines_sorted_by_parsimony_change[i][3])
	if not (lines_sorted_by_parsimony_change[i][6] in nodes_to_get):
		nodes_to_get.append(lines_sorted_by_parsimony_change[i][6])
	
	nodes_by_example[i].append(lines_sorted_by_parsimony_change[i][0])
	nodes_by_example[i].append(lines_sorted_by_parsimony_change[i][3])
	nodes_by_example[i].append(lines_sorted_by_parsimony_change[i][6])



samples_to_get = []
samples_by_example = [[] for x in range(how_many_to_see)]


samplelines = []
for file_name in sampleinfo_file_names:
	samplefile = open(file_name, 'r')
	samplelines.extend(samplefile.readlines())
	samplefile.close()




del samplelines[0]

for i in range(len(samplelines)):
	samplelines[i] = samplelines[i].split('\t')

	samples = samplelines[i][1].split(',')

	if samplelines[i][0] in nodes_to_get:
		for k in range(len(samples)):
			if not samples[k] in samples_to_get:
				samples_to_get.append(samples[k])
		
		for j in range(how_many_to_see):
			if samplelines[i][0] in nodes_by_example[j]:
				for k in range(len(samples)):
					if not samples[k] in samples_by_example[j]:
						samples_by_example[j].append(samples[k])




all_samples = open(all_fasta_file_name, 'r')
sample_lines = all_samples.readlines()

reference = open(reference_file_name, 'r')
reference_lines = reference.readlines()

reference.close()




singleton_files = []
singleton_file_lines = [[] for x in range(how_many_to_see)]


for i in range(how_many_to_see):

	for j in range(len(reference_lines)):
		singleton_file_lines[i].append(reference_lines[j])
	singleton_file_lines[i][len(singleton_file_lines[i]) - 1] = singleton_file_lines[i][len(singleton_file_lines[i]) - 1] + '\n'






in_relevent_sample = False
sample_is_relevent_to_examples = []

for i in range(len(sample_lines)):
	possible_name = sample_lines[i]
	possible_name.strip()

	if possible_name[0] == '>':
		in_relevent_sample = False
		sample_is_relevent_to_examples = []

		possible_name = possible_name[1:-1]

		#Accounting for coords that jalview adds to some samplenames
		if possible_name.rfind('/') > possible_name.rfind('|'):
			possible_name = possible_name[:-8]
			

		

		if possible_name in samples_to_get:
			
			in_relevent_sample = True
			for j in range(how_many_to_see):
				if possible_name in samples_by_example[j]:
					sample_is_relevent_to_examples.append(j)
					
	

	if in_relevent_sample:
		
		for example in sample_is_relevent_to_examples:
			singleton_file_lines[example].append(sample_lines[i])





for i in range(how_many_to_see):
	a_file = open("filtering/fastas/OrderedRecombs/%d.fa" % i, 'w')
	
	for line in singleton_file_lines[i]:
		a_file.write(line)
	a_file.close()


