import argparse
import sys
from TCRMatch.mait_match_cython import tcrmatch

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
        prog='match', description='Find matching CDR3beta sequences based on TCRMatch')
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
    args = prsr.parse_args(sys.argv[2:])
    
    if args.p:
        tcrmatch(args.i, args.o, args.p)
    else:
        tcrmatch(args.i, args.o, 1)
