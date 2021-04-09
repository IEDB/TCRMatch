import csv
import sys

output_map = {}
with open("../data/IEDB_data.tsv", "r") as inf:
    #Skip header
    inf.readline()
    for line in inf:
        line_l = line.rstrip().split("\t")
        seq = line_l[0]
        # We don't use untrimmed sequence currently but may be useful later
        orig_seq = line_l[1]
        receptor_group = line_l[2]
        epitopes = line_l[3]
        try:
            org = line_l[4]
            antigen = line_l[5]
        except:
            org = " "
            antigen = " "
        # This simplifies receptor group output by creating lists that are later joined
        if seq in output_map:
            output_map[seq].append((receptor_group, epitopes, org, antigen))
        else:
            output_map[seq] = [(receptor_group, epitopes, org, antigen)]


input_file = sys.argv[1]
tcrmatch_output = []
with open(input_file, "r") as infile:
    for line in infile:
        tcrmatch_output.append((line.rstrip().split(" ")))

with open(sys.argv[2], "w") as outf:
    outf.write("input_sequence\tmatch_sequence\tscore\tepitopes\treceptor_group\tantigen\tsource_organism\n")
    for match in tcrmatch_output:
        cur_grp = ','.join([x[0] for x in output_map[match[1]]])
        cur_epi = ','.join([x[1] for x in output_map[match[1]]])
        cur_org = ','.join([x[2] for x in output_map[match[1]]])
        cur_anti = ','.join([x[3] for x in output_map[match[1]]])
        outf.write(match[0] + "\t" + match[1] + "\t" + match[2]  + "\t" + cur_epi + "\t" + cur_grp + "\t" + cur_anti + "\t" + cur_org + "\n")