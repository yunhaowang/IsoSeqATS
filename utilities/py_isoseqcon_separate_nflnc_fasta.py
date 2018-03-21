#!/usr/bin/env python
import sys,re,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	separate_nflnc_fa(args.input,args.sense,args.unknown)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def separate_nflnc_fa(input_fa,output_sense,output_unknown):
	flag = 0
	for line in input_fa:
		if line.startswith(">"):
			if "strand=NA;" in line:
				flag = 0
				print >>output_unknown, line,
			else:
				flag = 1
				print >>output_sense, line,
		else:
			if flag == 0:
				print >>output_unknown, line,
			else:
				print >>output_sense, line,

	input_fa.close()
	output_sense.close()
	output_unknown.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Function: separate the nflnc fasta file")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: nflnc fasta file (i.e. 'nflnc.fasta' under the  output folder of PacBio SMRTanalysis2.3 'pbtranscript.py classify' with the parameter '--detect_chimera_nfl -d '")
	parser.add_argument('-s','--sense',type=argparse.FileType('w'),required=True,help="Output: nflnc fasta with sense strand")
	parser.add_argument('-u','--unknown',type=argparse.FileType('w'),required=True,help="Output: nflnc fasta with unknown strand")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
