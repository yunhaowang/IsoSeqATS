#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	extract_last_exon(args.input,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def extract_last_exon(input_gpd,output_txt):
	for line in input_gpd:
		gene,iso,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")[:11]
		if strand == "+":
			print >>output_txt, "\t".join([gene,iso,chr,strand,exon_start.split(",")[-2],exon_end.split(",")[-2]])
		elif strand == "-":
			print >>output_txt, "\t".join([gene,iso,chr,strand,exon_start.split(",")[0],exon_end.split(",")[0]])
		else:
			pass
	input_gpd.close()
	output_txt.close()

def do_inputs():
	output_txt_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. start site of last exon (+)
6. end site of last exon (+)'''

	parser = argparse.ArgumentParser(description="Function: extract last exon position from gpd file",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: txt file including last exon information")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
