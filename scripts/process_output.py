#!/usr/bin/env python3

import argparse
import csv
import sys

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', '--results_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)

    args = parser.parse_args()

    results_file = args.results_file
    output_file = args.output_file
    
    return results_file, output_file

def get_output_map():
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

    return output_map

results_file, output_file = parse_arguments()
output_map = get_output_map()

tcrmatch_output = []
with open(results_file, "r") as infile:
    for line in infile:
        tcrmatch_output.append((line.rstrip().split(" ")))

with open(output_file, "w") as outf:
    outf.write("input_sequence\tmatch_sequence\tscore\tepitopes\treceptor_group\tantigen\tsource_organism\n")
    for match in tcrmatch_output:
        cur_grp = ','.join([x[0] for x in output_map[match[1]]])
        cur_epi = ','.join([x[1] for x in output_map[match[1]]])
        cur_org = ','.join([x[2] for x in output_map[match[1]]])
        cur_anti = ','.join([x[3] for x in output_map[match[1]]])
        outf.write(match[0] + "\t" + match[1] + "\t" + match[2]  + "\t" + cur_epi + "\t" + cur_grp + "\t" + cur_anti + "\t" + cur_org + "\n")
