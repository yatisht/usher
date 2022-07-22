import sys
import string
import argparse
from collections import Counter

recombination_file_name = "filtering/data/combinedCatOnlyBestWithPVals.txt"
sampleinfo_file_name = "filtering/data/sampleInfo.txt"
relevent_sites_file_name = "filtering/data/allRelevantNodesInfSites.txt"

F = open(recombination_file_name, 'r')
lines = F.readlines()
F.close()

labels = lines[0].split('\t')
for i in range(len(labels)):
    labels[i] = labels[i].strip()

del lines[0]

parsimony_change = []

for i in range(len(lines)):
	lines[i] = lines[i].split('\t')
	
	parsimony_change.append(int(lines[i][10]) - int(lines[i][11]))

lines_sorted_by_parsimony_change = [x for _, x in sorted(zip(parsimony_change, lines) , reverse=True)]

index_to_check = int(sys.argv[1])


recomb_id   = int(lines_sorted_by_parsimony_change[index_to_check][0])
donor_id    = int(lines_sorted_by_parsimony_change[index_to_check][3])
acceptor_id = int(lines_sorted_by_parsimony_change[index_to_check][6])

#breakpoints
break1_1 = int(lines_sorted_by_parsimony_change[index_to_check][1].split(',')[0][1:])
break1_2 = int(lines_sorted_by_parsimony_change[index_to_check][1].split(',')[1][:-1])
break2_1 = int(lines_sorted_by_parsimony_change[index_to_check][2].split(',')[0][1:])
break2_2 = int(lines_sorted_by_parsimony_change[index_to_check][2].split(',')[1][:-1])


samplefile = open(sampleinfo_file_name, 'r')
samplelines = samplefile.readlines()
samplefile.close()

del samplelines[0]

recomb_sample = ""
donor_sample = ""
acceptor_sample = ""

recomb_samples = []
donor_samples = []
acceptor_samples = []

recomb_mutations = []
donor_mutations = []
acceptor_mutations = []


for i in range(len(samplelines)):
    samplelines[i] = samplelines[i].split('\t')

    samples = samplelines[i][1].split(',')


    if int(samplelines[i][0]) == recomb_id:
        recomb_sample = samples[0]
        recomb_samples = samples[:10]
        recomb_mutations = samplelines[i][2].split(',')

    if int(samplelines[i][0]) == donor_id:
        donor_sample = samples[0]
        donor_samples = samples[:10]
        donor_mutations = samplelines[i][2].split(',')

    if int(samplelines[i][0]) == acceptor_id:
        acceptor_sample = samples[0]
        acceptor_samples = samples[:10]
        acceptor_mutations = samplelines[i][2].split(',')


inf_file = open(relevent_sites_file_name, 'r')
inf_lines = inf_file.readlines()
inf_file.close()

mutagens = []

for line in inf_lines:
    line = line.split('\t')
    if int(line[0]) == recomb_id and int(line[1]) == donor_id and int(line[2]) == acceptor_id:
        mutagens = line[7].split(',')
        break


mutagens = [int(x) for x in mutagens]
mutagens.sort()

mutations = list(set().union(recomb_mutations, donor_mutations, acceptor_mutations))
mutations = [int(x) for x in mutations]
mutations.sort()

if len(mutagens) > 0:
    mutations = mutagens
    recomb_mutations = []
    donor_mutations = []
    acceptor_mutations = []
    
mutations_new_coords = []
for x in mutations:
    mutations_new_coords.append(int(x))

pre_recomb_mutations = []
pre_donor_mutations = []
pre_acceptor_mutations = []

trecomb_mutations = []
for x in recomb_mutations:
    pre_recomb_mutations.append(int(x))
    trecomb_mutations.append(int(x))
    
tdonor_mutations = []
for x in donor_mutations:
    pre_donor_mutations.append(int(x))
    tdonor_mutations.append(int(x))

tacceptor_mutations = []
for x in acceptor_mutations:
    pre_acceptor_mutations.append(int(x))
    tacceptor_mutations.append(int(x))


dna = open("filtering/fastas/AlignedRecombs/%d.fa" % index_to_check)
dna_lines = dna.readlines()

base_count = 0

#Matching mutations to reference
for i in range(1, len(dna_lines)):
    line = dna_lines[i].strip()
    
    if line[0] == '>':
        break
    
    for j in range(len(line)):
        if line[j] == '-':
            for k in range(len(mutations_new_coords)):
                if (mutations_new_coords[k] >= base_count):
                    mutations_new_coords[k] += 1
            for k in range(len(trecomb_mutations)):
                if (trecomb_mutations[k] >= base_count):
                    trecomb_mutations[k] += 1 
            for k in range(len(tdonor_mutations)):
                if (tdonor_mutations[k] >= base_count):
                    tdonor_mutations[k] += 1
            for k in range(len(tacceptor_mutations)):
                if (tacceptor_mutations[k] >= base_count):
                    tacceptor_mutations[k] += 1
        base_count += 1


trecomb_mutations.sort()
tdonor_mutations.sort()
tacceptor_mutations.sort()

recomb_mutations = trecomb_mutations
donor_mutations = tdonor_mutations
acceptor_mutations = tacceptor_mutations


base_count = 0

in_recomb = False
in_donor = False
in_acceptor = False


recomb_reads = [[] for x in mutations]
donor_reads = [[] for x in mutations]
acceptor_reads = [[] for x in mutations]

trecomb_reads   = [[[] for y in mutations] for x in recomb_samples]
tdonor_reads    = [[[] for y in mutations] for x in donor_samples]
tacceptor_reads = [[[] for y in mutations] for x in acceptor_samples]

recomb_index   = 0
donor_index    = 0
acceptor_index = 0



#grabbing reads, nearest 50 bp to mutation sites
for i in range(len(dna_lines)):
    possible_name = dna_lines[i]

    possible_name.strip()

    if in_recomb or in_donor or in_acceptor:
        base_list = dna_lines[i]
        base_list = base_list.strip()
        base_list = [char for char in base_list]
        
        

        for j in range(len(base_list)):
            base_count += 1
            for h in range(len(mutations)):
                if int(mutations_new_coords[h]) + 50 >= base_count and int(mutations_new_coords[h]) - 50 <= base_count:
                    if in_recomb:
                        #recomb_reads[h].append(base_list[j])
                        trecomb_reads[recomb_index][h].append(base_list[j])
                    if in_donor:
                        #donor_reads[h].append(base_list[j])
                        tdonor_reads[donor_index][h].append(base_list[j])
                    if in_acceptor:
                        #acceptor_reads[h].append(base_list[j])
                        tacceptor_reads[acceptor_index][h].append(base_list[j])
    
    if possible_name[0] == '>':
        base_count = 0
        in_recomb = False
        in_donor = False
        in_acceptor = False

        possible_name = possible_name[1:-1]

        #Accounting for coords that jalview adds to some samplenames
        if possible_name.rfind('/') > possible_name.rfind('|'):
            possible_name = possible_name[:-8]

        if possible_name in recomb_samples:
            recomb_index = recomb_samples.index(possible_name)
            in_recomb = True
        if possible_name in donor_samples:
            donor_index = donor_samples.index(possible_name)
            in_donor = True
        if possible_name in acceptor_samples:
            acceptor_index = acceptor_samples.index(possible_name)
            in_acceptor = True

        if possible_name == recomb_sample:
            in_recomb = True
        if possible_name == donor_sample:
            in_donor = True
        if possible_name == acceptor_sample:
            in_acceptor = True
    



#removing empty reads
recomb_sample_missing = [True for x in recomb_samples]
donor_sample_missing = [True for x in donor_samples]
acceptor_sample_missing = [True for x in acceptor_samples]


temp_recomb_reads   = []
temp_donor_reads    = []
temp_acceptor_reads = []

for i in range(len(trecomb_reads)):
    if len(trecomb_reads[i][0]) > 1:
        temp_recomb_reads.append(trecomb_reads[i])
        recomb_sample_missing[i] = False
trecomb_reads = temp_recomb_reads

for i in range(len(tdonor_reads)):
    if len(tdonor_reads[i][0]) > 1:
        temp_donor_reads.append(tdonor_reads[i])
        donor_sample_missing[i] = False
tdonor_reads = temp_donor_reads

for i in range(len(tacceptor_reads)):
    if len(tacceptor_reads[i][0]) > 1:
        temp_acceptor_reads.append(tacceptor_reads[i])
        acceptor_sample_missing[i] = False
tacceptor_reads = temp_acceptor_reads


for i in range(len(trecomb_reads)):
    for j in range(len(trecomb_reads[0])):
        trecomb_reads[i][j] = trecomb_reads[i][j][:101]

for i in range(len(tdonor_reads)):
    for j in range(len(tdonor_reads[0])):
        tdonor_reads[i][j] = tdonor_reads[i][j][:101]

for i in range(len(tacceptor_reads)):
    for j in range(len(tacceptor_reads[0])):
        tacceptor_reads[i][j] = tacceptor_reads[i][j][:101]


error_string = "%d\t%d\t%d\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%d\t" % (recomb_id, donor_id, acceptor_id, index_to_check)
error = False

if len(trecomb_reads) == 0:
    error_string += "Missing all recomb samples\t"
    print("Missing all recomb samples")
    error = True

if len(tdonor_reads) == 0:
    error_string += "Missing all donor samples\t"
    print("Missing all donor samples")
    error = True

if len(tacceptor_reads) == 0:
    error_string += "Missing all acceptor samples\t"
    print("Missing all acceptor samples")
    error = True

error_string += "\n"


if error:
    if(len(sys.argv) > 2):
        report_file = open("filtering/data/report.txt", 'a')
        report_file.write(error_string)
        report_file.close()

    sys.exit()



#Creating consensus sequence
for mutation in range(len(trecomb_reads[0])):

    for b in range(len(trecomb_reads[0][mutation])):
        

        base_over_samples = []

        for i in range(len(trecomb_reads)):
            base = trecomb_reads[i][mutation][b]
            base_over_samples.append(base)
        
        c = Counter(base_over_samples)
        common = c.most_common(2)
        base = common[0][0]
        
        
        if (base == 'n' or base == '-') and len(common) > 1:
            next_best_base = common[1][0]

            if next_best_base != 'n' and next_best_base != '-' and common[0][1] == common[1][1]:
                base = next_best_base

        recomb_reads[mutation].append(base)

for mutation in range(len(tdonor_reads[0])):

    for b in range(len(tdonor_reads[0][mutation])):

        base_over_samples = []

        for i in range(len(tdonor_reads)):
            base = tdonor_reads[i][mutation][b]
            base_over_samples.append(base)

        c = Counter(base_over_samples)
        common = c.most_common(2)
        base = common[0][0]
        
        
        if (base == 'n' or base == '-') and len(common) > 1:
            next_best_base = common[1][0]

            if next_best_base != 'n' and next_best_base != '-' and common[0][1] == common[1][1]:
                base = next_best_base
        
        donor_reads[mutation].append(base)


for mutation in range(len(tacceptor_reads[0])):

    for b in range(len(tacceptor_reads[0][mutation])):

        base_over_samples = []

        for i in range(len(tacceptor_reads)):
            base = tacceptor_reads[i][mutation][b]
            base_over_samples.append(base)

        c = Counter(base_over_samples)
        common = c.most_common(2)
        base = common[0][0]
        
        
        if (base == 'n' or base == '-') and len(common) > 1:
            next_best_base = common[1][0]

            if next_best_base != 'n' and next_best_base != '-' and common[0][1] == common[1][1]:
                base = next_best_base
        
        acceptor_reads[mutation].append(base)
            


#finding specific mutations
if len(mutagens) > 0:
    recomb_mutations = []
    donor_mutations = []
    acceptor_mutations = []

    for i in range(len(mutations)):
        if recomb_reads[i][50] != donor_reads[i][50] and recomb_reads[i][50] != acceptor_reads[i][50]:
            recomb_mutations.append(mutations[i])
    
        if donor_reads[i][50] != recomb_reads[i][50] and donor_reads[i][50] != acceptor_reads[i][50]:
            donor_mutations.append(mutations[i])
    
        if acceptor_reads[i][50] != recomb_reads[i][50] and acceptor_reads[i][50] != donor_reads[i][50]:
            acceptor_mutations.append(mutations[i])




nearest_weirdness = [50] * len(mutations)

#region 1 is before breakpoint 1
#region 2 is within breakpoint 1 interval
#region 3 is after breakpoint 1, before breakpoint 2
#region 4 is within breakpoint 1 interval
#region 5 is after breakpoint 2

acptr_matches_in_region_1 = 0
donor_matches_in_region_1 = 0
acptr_matches_in_region_2 = 0
donor_matches_in_region_2 = 0
acptr_matches_in_region_3 = 0
donor_matches_in_region_3 = 0
acptr_matches_in_region_4 = 0
donor_matches_in_region_4 = 0
acptr_matches_in_region_5 = 0
donor_matches_in_region_5 = 0


INDEL_in_wrong_place = False
too_many_mutations_near_indel = False
mutations_on_indels = 0



for i in range(len(mutations)):

    #finding nearest weirdness
    for j in range(51):
        not_all_dashes_p = recomb_reads[i][50+j] != '-' or donor_reads[i][50+j] != '-' or acceptor_reads[i][50+j] != '-'
        not_all_dashes_n = recomb_reads[i][50-j] != '-' or donor_reads[i][50-j] != '-' or acceptor_reads[i][50-j] != '-'

        if (recomb_reads[i][50+j] == '-' or recomb_reads[i][50+j] == 'n') and not_all_dashes_p:
            nearest_weirdness[i] = min(nearest_weirdness[i], j)
        if (recomb_reads[i][50-j] == '-' or recomb_reads[i][50-j] == 'n') and not_all_dashes_n:
            nearest_weirdness[i] = min(nearest_weirdness[i], j)
        if (donor_reads[i][50+j] == '-' or donor_reads[i][50+j] == 'n') and not_all_dashes_p:
            nearest_weirdness[i] = min(nearest_weirdness[i], j)
        if (donor_reads[i][50-j] == '-' or donor_reads[i][50-j] == 'n') and not_all_dashes_n:
            nearest_weirdness[i] = min(nearest_weirdness[i], j)
        if (acceptor_reads[i][50+j] == '-' or acceptor_reads[i][50+j] == 'n') and not_all_dashes_p:
            nearest_weirdness[i] = min(nearest_weirdness[i], j)
        if (acceptor_reads[i][50-j] == '-' or acceptor_reads[i][50-j] == 'n') and not_all_dashes_n:
            nearest_weirdness[i] = min(nearest_weirdness[i], j)

    weird = nearest_weirdness[i] < 5

    if nearest_weirdness[i] == 0:
        if recomb_reads[i][50] != '-' or donor_reads[i][50] != '-' or acceptor_reads[i][50] != '-':
            mutations_on_indels += 1

    print()
    where_we_are_string = ""

    in_good_region = 200 < mutations[i] and mutations[i] < 29800

    if(mutations[i] <= break1_1):
        where_we_are_string = "before breakpoint 1, should match acceptor"
        donor_matches_in_region_1 += (recomb_reads[i][50] == donor_reads[i][50]) and not weird
        acptr_matches_in_region_1 += (recomb_reads[i][50] == acceptor_reads[i][50]) and not weird


        if recomb_reads[i][50] != '-' and recomb_reads[i][50] != 'n' and in_good_region:
            INDEL = True
            for j in range(-5,6):
                if acceptor_reads[i][50 + j] != '-' and acceptor_reads[i][50 + j] != 'n':
                    INDEL = False
            if INDEL:
                INDEL_in_wrong_place = True
        
        if acceptor_reads[i][50] != '-' and acceptor_reads[i][50] != 'n' and in_good_region:
            INDEL = True
            for j in range(-5,6):
                if recomb_reads[i][50 + j] != '-' and recomb_reads[i][50 + j] != 'n':
                    INDEL = False
            if INDEL:
                INDEL_in_wrong_place = True

    if(break1_1 < mutations[i] and mutations[i] < break1_2):
        where_we_are_string = "within breakpoint 1 interval"
        donor_matches_in_region_2 += (recomb_reads[i][50] == donor_reads[i][50]) and not weird
        acptr_matches_in_region_2 += (recomb_reads[i][50] == acceptor_reads[i][50]) and not weird

    if(break1_2 <= mutations[i] and mutations[i] <= break2_1):
        where_we_are_string = "after breakpoint 1, before breakpoint 2, should match donor"
        donor_matches_in_region_3 += (recomb_reads[i][50] == donor_reads[i][50]) and not weird
        acptr_matches_in_region_3 += (recomb_reads[i][50] == acceptor_reads[i][50]) and not weird
        

        if recomb_reads[i][50] != '-' and recomb_reads[i][50] != 'n' and in_good_region:
            INDEL = True
            for j in range(-5,6):
                if donor_reads[i][50 + j] != '-' and donor_reads[i][50 + j] != 'n':
                    INDEL = False
            if INDEL:
                INDEL_in_wrong_place = True
        
        if donor_reads[i][50] != '-' and donor_reads[i][50] != 'n' and in_good_region:
            INDEL = True
            for j in range(-5,6):
                if recomb_reads[i][50 + j] != '-' and recomb_reads[i][50 + j] != 'n':
                    INDEL = False
            if INDEL:
                INDEL_in_wrong_place = True

    if(break2_1 < mutations[i] and mutations[i] < break2_2):
        where_we_are_string = "within breakpoint 2 interval"
        donor_matches_in_region_4 += (recomb_reads[i][50] == donor_reads[i][50]) and not weird
        acptr_matches_in_region_4 += (recomb_reads[i][50] == acceptor_reads[i][50]) and not weird

    if(break2_2 <= mutations[i]):
        where_we_are_string = "after breakpoint 2, should match acceptor"
        donor_matches_in_region_5 += (recomb_reads[i][50] == donor_reads[i][50]) and not weird
        acptr_matches_in_region_5 += (recomb_reads[i][50] == acceptor_reads[i][50]) and not weird

        if recomb_reads[i][50] != '-' and recomb_reads[i][50] != 'n' and in_good_region:
            INDEL = True
            for j in range(-5,6):
                if acceptor_reads[i][50 + j] != '-' and acceptor_reads[i][50 + j] != 'n':
                    INDEL = False
            if INDEL:
                INDEL_in_wrong_place = True
        
        if acceptor_reads[i][50] != '-' and acceptor_reads[i][50] != 'n' and in_good_region:
            INDEL = True
            for j in range(-5,6):
                if recomb_reads[i][50 + j] != '-' and recomb_reads[i][50 + j] != 'n':
                    INDEL = False
            if INDEL:
                INDEL_in_wrong_place = True

    

    print("mutation site", mutations[i], where_we_are_string)
    print("\t\t", "                                                  v")

    recomb_read = ""
    recomb_read = recomb_read.join(recomb_reads[i])
    print(recomb_id, "\trecomb\t", recomb_read)
    donor_read = ""
    donor_read = donor_read.join(donor_reads[i])
    print(donor_id, "\tdonor\t", donor_read)
    acceptor_read = ""
    acceptor_read = acceptor_read.join(acceptor_reads[i])
    print(acceptor_id, "\taccptr\t", acceptor_read)
    print()
    for j in range(len(trecomb_reads)):
        recomb_read = ""
        recomb_read = recomb_read.join(trecomb_reads[j][i])
        print(recomb_id, "\trecomb\t", recomb_read)
    for j in range(len(tdonor_reads)):
        donor_read = ""
        donor_read = donor_read.join(tdonor_reads[j][i])
        print(donor_id, "\tdonor\t", donor_read)
    for j in range(len(tacceptor_reads)):
        acceptor_read = ""
        acceptor_read = acceptor_read.join(tacceptor_reads[j][i])
        print(acceptor_id, "\taccptr\t", acceptor_read)


print()
print("Nearest weirdness to mutations:")
print("site\tNearest weirdness")

weird_mutations = 0

for i in range(len(mutations)):
    print(mutations[i], '\t', nearest_weirdness[i])
    if nearest_weirdness[i] < 2:
        weird_mutations += 1
    if nearest_weirdness[i] <= 5:
        
        bad = 0
        iffy = 0

        for j in range(5):
            dashes = 0
            ns = 0
            if recomb_reads[i][50+nearest_weirdness[i] + j] == '-':
                dashes += 1
            if donor_reads[i][50+nearest_weirdness[i] + j] == '-':
                dashes += 1
            if acceptor_reads[i][50+nearest_weirdness[i] + j] == '-':
                dashes += 1
            if recomb_reads[i][50+nearest_weirdness[i] + j] == 'n':
                ns += 1
            if donor_reads[i][50+nearest_weirdness[i] + j] == 'n':
                ns += 1
            if acceptor_reads[i][50+nearest_weirdness[i] + j] == 'n':
                ns += 1
            if 0 < dashes and dashes < 3:
                bad += 1
            if dashes == 3:
                iffy +=1
            if 0 < ns:
                bad += 1
            
            dashes = 0
            ns = 0
            if recomb_reads[i][50-nearest_weirdness[i] - j] == '-':
                dashes += 1
            if donor_reads[i][50-nearest_weirdness[i] - j] == '-':
                dashes += 1
            if acceptor_reads[i][50-nearest_weirdness[i] - j] == '-':
                dashes += 1
            if recomb_reads[i][50-nearest_weirdness[i] - j] == 'n':
                ns += 1
            if donor_reads[i][50-nearest_weirdness[i] - j] == 'n':
                ns += 1
            if acceptor_reads[i][50-nearest_weirdness[i] - j] == 'n':
                ns += 1
            if 0 < dashes and dashes < 3:
                bad += 1
            if dashes == 3:
                iffy +=1
            if 0 < ns:
                bad += 1
        
        if bad + iffy >= 5 and iffy < 5:
            too_many_mutations_near_indel = True
            




if weird_mutations >= 6:
    too_many_mutations_near_indel = True
if mutations_on_indels >= 1:
    too_many_mutations_near_indel = True


#matches are only with no weirdness 4 bp near
print()
print("Matches:")
print("region\tdonor\taccptr\tshould skew")
print("1\t", donor_matches_in_region_1,"\t", acptr_matches_in_region_1, "\tacceptor")
print("2\t", donor_matches_in_region_2,"\t", acptr_matches_in_region_2)
print("3\t", donor_matches_in_region_3,"\t", acptr_matches_in_region_3, "\tdonor")
print("4\t", donor_matches_in_region_4,"\t", acptr_matches_in_region_4)
print("5\t", donor_matches_in_region_5,"\t", acptr_matches_in_region_5, "\tacceptor")

print()
for j in range(len(recomb_samples)):
    print("recomb sample", recomb_samples[j] , "\tMISSING" if recomb_sample_missing[j] else "")
for j in range(len(donor_samples)):
    print("donor  sample", donor_samples[j] , "\tMISSING" if donor_sample_missing[j]  else "")
for j in range(len(acceptor_samples)):
    print("accptr sample", acceptor_samples[j], "\tMISSING" if acceptor_sample_missing[j]  else "")
print()

for j in range(len(labels)):
	print(labels[j], lines_sorted_by_parsimony_change[index_to_check][j])



pre_b1_skew                 = acptr_matches_in_region_1 - donor_matches_in_region_1
btw_b1_b2_skew              = donor_matches_in_region_3 - acptr_matches_in_region_3
post_b2_skew                = acptr_matches_in_region_5 - donor_matches_in_region_5

not_enough_acceptor_matches = pre_b1_skew <= 0 and post_b2_skew <= 0
not_enough_donor_matches    = btw_b1_b2_skew <= 0

unnecessary_breakpoint = (pre_b1_skew < 0 and post_b2_skew > 0) or (pre_b1_skew > 0 and post_b2_skew < 0)


#detecting mutation clumping ########################################
                                                                    #
most_mutations_in_20 = 0
site_of_most_mutations = 0
lower_index_of_most_mutations = 0
index_of_most_mutations = 0
most_clumps_in = ""

for i in range(len(recomb_mutations)):
    mutations_in_20 = 1
    temp_lower_index = 0
    for j in range(i):
        if recomb_mutations[j] >= recomb_mutations[i] - 20:
            if mutations_in_20 == 1:
                temp_lower_index = j
            mutations_in_20 += 1
    if mutations_in_20 >= most_mutations_in_20:
        most_mutations_in_20 = mutations_in_20
        site_of_most_mutations = recomb_mutations[i]
        lower_index_of_most_mutations = temp_lower_index
        index_of_most_mutations = i
        most_clumps_in = "recomb"

for i in range(len(donor_mutations)):
    mutations_in_20 = 1
    temp_lower_index = 0
    for j in range(i):
        if donor_mutations[j] >= donor_mutations[i] - 20:
            if mutations_in_20 == 1:
                temp_lower_index = j
            mutations_in_20 += 1
    if mutations_in_20 >= most_mutations_in_20:
        most_mutations_in_20 = mutations_in_20
        site_of_most_mutations = donor_mutations[i]
        lower_index_of_most_mutations = temp_lower_index
        index_of_most_mutations = i
        most_clumps_in = "donor"

for i in range(len(acceptor_mutations)):
    mutations_in_20 = 1
    temp_lower_index = 0
    for j in range(i):
        if acceptor_mutations[j] >= acceptor_mutations[i] - 20:
            if mutations_in_20 == 1:
                temp_lower_index = j
            mutations_in_20 += 1
    if mutations_in_20 >= most_mutations_in_20:
        most_mutations_in_20 = mutations_in_20
        site_of_most_mutations = acceptor_mutations[i]
        lower_index_of_most_mutations = temp_lower_index
        index_of_most_mutations = i
        most_clumps_in = "acceptor"
                                                                    #
#####################################################################



#Find out size of indel near clump, near meaning less than 1 base pair away####
                                                                              #
indel_near_clump_size = 0
indel_near_clump = False

clump_reads = [[],[],[],[],[],[]]

for i in range(len(mutations)):
    if most_clumps_in == "recomb":
        if mutations[i] == recomb_mutations[lower_index_of_most_mutations]:
            clump_reads[0] = recomb_reads[i]
            clump_reads[1] = donor_reads[i]
            clump_reads[2] = acceptor_reads[i]
        if mutations[i] == recomb_mutations[index_of_most_mutations]:
            clump_reads[3] = recomb_reads[i]
            clump_reads[4] = donor_reads[i]
            clump_reads[5] = acceptor_reads[i]
    if most_clumps_in == "donor":
        if mutations[i] == donor_mutations[lower_index_of_most_mutations]:
            clump_reads[0] = recomb_reads[i]
            clump_reads[1] = donor_reads[i]
            clump_reads[2] = acceptor_reads[i]
        if mutations[i] == donor_mutations[index_of_most_mutations]:
            clump_reads[3] = recomb_reads[i]
            clump_reads[4] = donor_reads[i]
            clump_reads[5] = acceptor_reads[i]
    if most_clumps_in == "acceptor":
        if mutations[i] == acceptor_mutations[lower_index_of_most_mutations]:
            clump_reads[0] = recomb_reads[i]
            clump_reads[1] = donor_reads[i]
            clump_reads[2] = acceptor_reads[i]
        if mutations[i] == acceptor_mutations[index_of_most_mutations]:
            clump_reads[3] = recomb_reads[i]
            clump_reads[4] = donor_reads[i]
            clump_reads[5] = acceptor_reads[i]




temp_indel_near_clump_size = 0
for j in range(50):

    if j > 5 and not indel_near_clump:
        break

    lbase = clump_reads[0][50 - j]
    rbase = clump_reads[3][50 + j]

        
    if lbase == '-' or lbase == 'n' or rbase == '-' or rbase == 'n':
        if (j <= 5):
            indel_near_clump = True
            
        temp_indel_near_clump_size += 1

    else:
        if (j <= 5):
            temp_indel_near_clump_size = 0
        else:
            break
    
indel_near_clump_size = max(temp_indel_near_clump_size,indel_near_clump_size)
temp_indel_near_clump_size = 0
for j in range(50):

    if j > 5 and not indel_near_clump:
        break

    lbase = clump_reads[1][50 - j]
    rbase = clump_reads[4][50 + j]


        
    if lbase == '-' or lbase == 'n' or rbase == '-' or rbase == 'n':
        if (j <= 5):
            indel_near_clump = True
            
        temp_indel_near_clump_size += 1

    else:
        if (j <= 5):
            temp_indel_near_clump_size = 0
        else:
            break

indel_near_clump_size = max(temp_indel_near_clump_size,indel_near_clump_size)
temp_indel_near_clump_size = 0
for j in range(50):

    if j > 5 and not indel_near_clump:
        break

    lbase = clump_reads[2][50 - j]
    rbase = clump_reads[5][50 + j]

    if lbase == '-' or lbase == 'n' or rbase == '-' or rbase == 'n':
        if (j <= 5):
            indel_near_clump = True
            
        temp_indel_near_clump_size += 1

    else:
        if (j <= 5):
            temp_indel_near_clump_size = 0
        else:
            break
    
indel_near_clump_size = max(temp_indel_near_clump_size,indel_near_clump_size)

                                                                              #
#indel_near_clump_size is populated############################################



suspicious_mutation_clump = most_mutations_in_20 >= 6 or (most_mutations_in_20 >= 3 and indel_near_clump_size >= 10)




most_informative_sites_in_20 = 0
informative_clump_site = 0

for i in range(len(mutations_new_coords)):
    informative_sites_in_20 = 1

    for j in range(i):
        if mutations_new_coords[j] >= mutations_new_coords[i] - 20:
            informative_sites_in_20 += 1
    
    if most_informative_sites_in_20 <= informative_sites_in_20:
        most_informative_sites_in_20 = informative_sites_in_20
        informative_clump_site = mutations_new_coords[i]

informative_sites_clump = most_informative_sites_in_20 >= 6




should_check = False

if not_enough_acceptor_matches:
    should_check = True
if not_enough_donor_matches:
    should_check = True
if unnecessary_breakpoint:
    should_check = True

if too_many_mutations_near_indel:
    should_check = True

if suspicious_mutation_clump:
    should_check = True

if informative_sites_clump:
    should_check = True

if INDEL_in_wrong_place:
    should_check = True



print()
print("Report:")
print("recomb_node_id")
print("donor_node_id")
print("acceptor_node_id")
print("should check") #3
print("Not enough acceptor matches")
print("Not enough donor matches")
print("Unnecessary breakpoint")
print("pre_b1_skew") #
print("btw_b1_b2_skew") #
print("post_b2_skew") #
print("Too many mutations near INDEL")
print("Mutations near INDEL") #
print("Mutations on INDEL")
print("Suspicious mutation clump")
print("Mutations clumped") #
print("Mutation clump site") #
print("Clump in") #
print("Size of INDEL near clump") #

print("Informative sites clump")
print("Informative sites clumped")
print("site of informative sites clump")

print("INDEL in the wrong region") #
print(labels[1])
print(labels[2])
print(labels[9])
print(labels[10])
print(labels[11])
print(labels[12])
print(labels[13]) #
#print(labels[14])
#print(labels[15])
print("Parsimony change rank") #
print()

to_report = [recomb_id, donor_id, acceptor_id, should_check,
not_enough_acceptor_matches, not_enough_donor_matches, unnecessary_breakpoint,
pre_b1_skew, btw_b1_b2_skew, post_b2_skew,
too_many_mutations_near_indel, weird_mutations, mutations_on_indels,
suspicious_mutation_clump, most_mutations_in_20, site_of_most_mutations, most_clumps_in, indel_near_clump_size,
informative_sites_clump, most_informative_sites_in_20, informative_clump_site,
INDEL_in_wrong_place,
lines_sorted_by_parsimony_change[index_to_check][1],
lines_sorted_by_parsimony_change[index_to_check][2],
lines_sorted_by_parsimony_change[index_to_check][9],
lines_sorted_by_parsimony_change[index_to_check][10],
lines_sorted_by_parsimony_change[index_to_check][11],
lines_sorted_by_parsimony_change[index_to_check][12],
lines_sorted_by_parsimony_change[index_to_check][13],
#lines_sorted_by_parsimony_change[index_to_check][14],
#lines_sorted_by_parsimony_change[index_to_check][15],
index_to_check, '\n']

report_str = '\t'.join([str(x) for x in to_report])

print(report_str)


output_to_report_file = len(sys.argv) > 2
if not output_to_report_file:
    sys.exit()

report_file = open("filtering/data/report.txt", 'a')
report_file.write(report_str)
report_file.close()
