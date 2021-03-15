import argparse
import sys
import os
import math
from tcrmatch_c import tcrmatch
import multiprocessing as mp
import numpy as np
import urllib.request

package_dir = os.path.dirname(os.path.abspath(__file__))
iedb_file = os.path.join(package_dir, 'data/IEDB_data.tsv')

# Parse and IEDB data and prep for additional output
# This is added here to make multiprocessing easier
output_map = {}
iedb_seqs = set()
with open(iedb_file, "r") as inf:
    #Skip header
    inf.readline()
    for line in inf:
        line_l = line.rstrip().split("\t")
        # The least elegant way to unpack this line
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
        iedb_seqs.add(seq)
        # This simplifies receptor group output by creating lists that are later joined
        if seq in output_map:
            output_map[seq].append((receptor_group, epitopes, org, antigen))
        else:
            output_map[seq] = [(receptor_group, epitopes, org, antigen)]

def parse_input(infile, fformat="text"):
    input_seqs = []

    if fformat == "airr":
        with open(args.i, "r") as inf:
            inf.readline() #Skip header
            for line in inf:
                items = line.rstrip().split("\t")
                if items[3] == "T":
                    input_seqs.append(items[48])
    elif fformat == "text":
        with open(args.i, "r") as inf:
            for line in inf:
                input_seqs.append(line.rstrip())
    else:
        raise ValueError("Invalid file format provided. Options are \"text\" or \"airr\"")

    return input_seqs

def run_tcrmatch(input_seqs, iedb_seqs=iedb_seqs):
    res = []
    for in_seq in input_seqs:
        for iedb_seq in iedb_seqs:
            res.append(tcrmatch(in_seq, iedb_seq))

    return res

tasks = ["match", "update"]
if len(sys.argv) < 2 or (sys.argv[1] not in tasks):
    usage = """
  ________________  __  ___      __       __  
 /_  __/ ____/ __ \/  |/  /___ _/ /______/ /_ 
  / / / /   / /_/ / /|_/ / __ `/ __/ ___/ __ \\
 / / / /___/ _, _/ /  / / /_/ / /_/ /__/ / / /
/_/  \____/_/ |_/_/  /_/\__,_/\__/\___/_/ /_/  (TCRMatch)
                       
	Usage: python -m TCRMatch <task> [options]

	Available tasks are:
		match\n\t\tupdate

	For help:
	> python -m TCRMatch <task> -h
	"""
    print(usage)
    sys.exit()

elif sys.argv[1] == "match":
    prsr = argparse.ArgumentParser(
        prog='match',
        description='Find matching CDR3beta sequences based on TCRMatch')
    prsr.add_argument(
        '-p',
        help="number of threads (default is 1)",
        type=int,
        metavar='num_threads',
        required=False)
    prsr.add_argument(
        '-i',
        help="input file containing a list of TCR CDR3beta sequences",
        type=str,
        metavar='infile_name',
        required=True)
    prsr.add_argument(
        '-o',
        help="output file name and location, ex: data/my_outfile.csv",
        type=str,
        metavar='outfile_name',
        required=True)
    prsr.add_argument(
        '-f',
        help="input file format (options are airr or text -- default is text)",
        type=str,
        metavar='input_format',
        default='text',
        required=False)
    prsr.add_argument(
        '-t',
        help="minimum threshold to constitute a match (default is .97)",
        type=float,
        metavar='threshold',
        default=.97,
        required=False)

    args = prsr.parse_args(sys.argv[2:])
    input_seqs = parse_input(args.i, args.f)

    if args.p:
        pool = mp.Pool(processes = args.p)
        res = pool.map(run_tcrmatch, np.array_split(input_seqs, args.p))
        pool.close()
        pool.join()
        with open(args.o, "w") as outf:
            outf.write("input_sequence\tmatch_sequence\tscore\tepitopes\treceptor_group\tantigen\tsource_organism\n")
            for chunk in res:
                for line in chunk:
                    if line[2] >= args.t:
                        cur_grp = ','.join([x[0] for x in output_map[line[1]]])
                        cur_epi = ','.join([x[1] for x in output_map[line[1]]])
                        cur_org = ','.join([x[2] for x in output_map[line[1]]])
                        cur_anti = ','.join([x[3] for x in output_map[line[1]]])
                        outf.write(line[0] + "\t" + line[1] + "\t" + "{:.2f}".format(line[2])  + "\t" + cur_epi + "\t" + cur_grp + "\t" + cur_anti + "\t" + cur_org + "\n")
        
    else:
        res = run_tcrmatch(input_seqs, iedb_seqs)
        with open(args.o, "w") as outf:
            outf.write("input_sequence\tmatch_sequence\tscore\tepitopes\treceptor_group\tantigen\tsource_organism\n")
            for line in res:
                if line[2] >= args.t:
                    cur_grp = ','.join([x[0] for x in output_map[line[1]]])
                    cur_epi = ','.join([x[1] for x in output_map[line[1]]])
                    cur_org = ','.join([x[2] for x in output_map[line[1]]])
                    cur_anti = ','.join([x[3] for x in output_map[line[1]]])
                    outf.write(line[0] + "\t" + line[1] + "\t" + "{:.2f}".format(line[2])  + "\t" + cur_epi + "\t" + cur_grp + "\t" + cur_anti + "\t" + cur_org + "\n")

elif sys.argv[1] == "update":
    print("Updating the IEDB data file")
    curr_iedb = urllib.request.urlopen('https://downloads.iedb.org/misc/TCRMatch/IEDB_data.tsv').read()
    print("...")
    with open(iedb_file, "wb") as fh:
        fh.write(curr_iedb)
    print("Done!")
