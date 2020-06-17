import argparse
import sys
import os
import math
from tcrmatch_c import tcrmatch_c
import multiprocessing as mp
import functools

package_dir = os.path.dirname(os.path.abspath(__file__))
iedb_file = os.path.join(package_dir, 'data/iedb_tcr.tsv')

def parse_input(infile, airr=False, multiproc=False):
    input_seqs = []
    if not airr and multiproc:
        with open(args.i, "r") as inf:
            for line in inf:
                input_seqs.append([line.rstrip().encode()])
    if not airr and not multiproc:
        with open(args.i, "r") as inf:
            for line in inf:
                input_seqs.append(line.rstrip().encode())
    if airr and multiproc:
        with open(args.i, "r") as inf:
            inf.readline() #Skip header
            for line in inf:
                items = line.rstrip().split("\t")
                if items[3] == "T":
                    input_seqs.append([items[48].encode()])
    if airr and not multiproc:
        with open(args.i, "r") as inf:
            inf.readline() #Skip header
            for line in inf:
                items = line.rstrip().split("\t")
                if items[3] == "T":
                    input_seqs.append(items[48].encode())
    return input_seqs

tasks = ["match"]
if len(sys.argv) < 2 or (sys.argv[1] not in tasks):
    usage = """
  ________________  __  ___      __       __  
 /_  __/ ____/ __ \/  |/  /___ _/ /______/ /_ 
  / / / /   / /_/ / /|_/ / __ `/ __/ ___/ __ \\
 / / / /___/ _, _/ /  / / /_/ / /_/ /__/ / / /
/_/  \____/_/ |_/_/  /_/\__,_/\__/\___/_/ /_/  (TCRMatch)
                       
	Usage: python -m TCRMatch <task> [options]

	Available tasks are:
		{:s}

	For help:
	> python -m TCRMatch <task> -h
	""".format(*tasks)
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
    args = prsr.parse_args(sys.argv[2:])



    # Parse and encode IEDB strings as binary strings
    iedb_seqs = []
    with open(iedb_file, "r") as inf:
        for line in inf:
            iedb_seqs.append(line.rstrip().encode())

    if args.f == 'text' and args.p:
        input_seqs = parse_input(args.i, False, True)
    elif args.f == 'airr' and args.p:
        input_seqs = parse_input(args.i, True, True)
    elif args.f == 'text' and not args.p:
        input_seqs = parse_input(args.i, False, False)
    elif args.f == 'airr' and not args.p:
        input_seqs = parse_input(args.i, True, False)
    else:
        input_seqs = parse_input(args.i, False, False)
        


    if args.p:
        pool = mp.Pool(processes = args.p)
        res = pool.map(functools.partial(tcrmatch_c, iedb_seqs), input_seqs)
        pool.close()
        pool.join()
        with open(args.o, "w") as outf:
            outf.write("input_sequence\tmatch_sequence\tscore\n")
            for chunk in res:
                for line in chunk:
                    outf.write(line[0].decode() + "\t" + line[1].decode() + "\t" + "{:.2f}".format(line[2]) + "\n")
    else:
        res = tcrmatch_c(input_seqs, iedb_seqs)
        with open(args.o, "w") as outf:
            outf.write("input_sequence\tmatch_sequence\tscore\n")
            for line in res:
                outf.write(line[0].decode() + "\t" + line[1].decode() + "\t" + "{:.2f}".format(line[2]) + "\n")

