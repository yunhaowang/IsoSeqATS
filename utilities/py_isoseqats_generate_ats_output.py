#!/usr/bin/env python
import sys,time,re,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	if args.cage_bed:
		dic_chr_strand_peak = parse_cage_bed(args.cage_bed)
	else:
		dic_chr_strand_peak = {}
	dic_lr_tss = parse_long_read_tss(args.lr_gpd)
	generate_ats_output(args.input,args.output,dic_chr_strand_peak,dic_lr_tss,args.dis_group)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def parse_cage_bed(input_cage_bed): #
	dic_chr_strand_peak = {}
	for line in input_cage_bed:
		chr,start,end,peak_name,score,strand = line.strip().split("\t")[:6]
		chr_strand = chr + "&" + strand
		if chr_strand not in dic_chr_strand_peak.keys():
			dic_chr_strand_peak[chr_strand] = {}
			dic_chr_strand_peak[chr_strand]["start"] = []
			dic_chr_strand_peak[chr_strand]["end"] = []
			dic_chr_strand_peak[chr_strand]["start"].append(int(start)+1)
			dic_chr_strand_peak[chr_strand]["end"].append(int(end))
		else:
			dic_chr_strand_peak[chr_strand]["start"].append(int(start)+1)
			dic_chr_strand_peak[chr_strand]["end"].append(int(end))
	input_cage_bed.close()
	return dic_chr_strand_peak

def parse_long_read_tss(input_lr_gpd):
	dic_lr_tss = {}
	for line in input_lr_gpd:
		read_id,chr,strand,tss,tts,mapq,sf_flag = line.strip().split("\t")[1:8]
		if read_id.split("_")[-1].startswith("F1"): # check 5 end completeness
			if strand == "+":
				if sf_flag.split("_")[0] == "0": # check soft-clip
					dic_lr_tss[read_id] = (int(tss) + 1)
			elif strand == "-":
				if sf_flag.split("_")[1] == "0": # check soft-clip
					dic_lr_tss[read_id] = tts
			else:
				pass
	input_lr_gpd.close()
	return dic_lr_tss

def generate_ats_output(input_con_gpd,output_ats_gpd,dic_chr_strand_peak,dic_lr_tss,gather_dis):
	for line in input_con_gpd:
		chr,strand = line.strip().split("\t")[2:4]
		read_set = line.strip().split("\t")[-1].split(",")
		dic_lr_tss_set = {}
		for read in read_set:
			if read.split("_")[-1].startswith("F0"): continue# remove the read with truncated 5'end
			if read in dic_lr_tss.keys():
				lr_tss = int(dic_lr_tss[read])
				if lr_tss not in dic_lr_tss_set.keys():
					dic_lr_tss_set[lr_tss] = 1
				else:
					dic_lr_tss_set[lr_tss] += 1
		if dic_lr_tss_set == {}: # No valid TSS
			print >>output_ats_gpd, "\t".join(["\t".join(line.strip().split("\t")[:-1]),"NA","NA"])
		else: # Have valid TSS
			pos_list = dic_lr_tss_set.keys()
			pos_list.sort()

			dic_tss_lr_cage = {}
			chr_strand = chr+"&"+strand
			if chr_strand in dic_chr_strand_peak.keys(): # cage
				for tss in pos_list:
					dic_tss_lr_cage[tss] = str(tss) + "_" + str(dic_lr_tss_set[tss]) + "_0"
					for i in range(0,len(dic_chr_strand_peak[chr_strand]["start"])):
						if tss >= dic_chr_strand_peak[chr_strand]["start"][i] and tss <= dic_chr_strand_peak[chr_strand]["end"][i]:
							dic_tss_lr_cage[tss] = str(tss) + "_" + str(dic_lr_tss_set[tss]) + "_1"
							break
			else:
				for tss in pos_list:
					dic_tss_lr_cage[tss] = str(tss) + "_" + str(dic_lr_tss_set[tss]) + "_0"

			if len(pos_list) == 1: # one tss
				idv_tss_set = dic_tss_lr_cage[pos_list[0]]
				grp_tss_set = idv_tss_set
				print >>output_ats_gpd, "\t".join(["\t".join(line.strip().split("\t")[:-1]),idv_tss_set,grp_tss_set])
			else:
				# individual tss
				idv_tss_list = []
				for idv_tss in pos_list:
					idv_tss_list.append(dic_tss_lr_cage[idv_tss])
				idv_tss_set = ",".join(idv_tss_list)

				# gather tss
				group_idx_set = []
				group_idx = [pos_list[0]]
				for i in range(1,len(pos_list)):
					if (pos_list[i] - gather_dis) > group_idx[-1]:
						group_idx_set.append(group_idx)
						group_idx = [pos_list[i]]
						if pos_list[i] == pos_list[-1]:
							group_idx_set.append(group_idx)
					else:
						group_idx.append(pos_list[i])
						if pos_list[i] == pos_list[-1]:
							group_idx_set.append(group_idx)
				
				# grouped tss
				grp_tss_list = []
				for gp in group_idx_set:
					sum_pos_frq_lr = 0 # long read data
					sum_pos_frq_sr = 0 # cage data
					dic_pos_frq_tss = {}
					for tss in gp:
						pos_frq = int(dic_tss_lr_cage[tss].split("_")[1]) + int(dic_tss_lr_cage[tss].split("_")[2])
						sum_pos_frq_lr += int(dic_tss_lr_cage[tss].split("_")[1])
						sum_pos_frq_sr += int(dic_tss_lr_cage[tss].split("_")[2])
						if pos_frq not in dic_pos_frq_tss.keys():
							dic_pos_frq_tss[pos_frq] = []
							dic_pos_frq_tss[pos_frq].append(tss)
						else:
							dic_pos_frq_tss[pos_frq].append(tss)
					max_pos_frq_tss = dic_pos_frq_tss[max(dic_pos_frq_tss.keys())][0]
					grp_tss_list.append(str(max_pos_frq_tss)+"_"+str(sum_pos_frq_lr)+"_"+str(sum_pos_frq_sr))
				grp_tss_set = ",".join(grp_tss_list)
				print >>output_ats_gpd, "\t".join(["\t".join(line.strip().split("\t")[:-1]),idv_tss_set,grp_tss_set])

	input_con_gpd.close()
	output_ats_gpd.close()

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS (+)
6. TTS (+)
7. number of support full-length long reads
8. number of support total long reads
9. exon count
10. exon start set
11. exon end set
12. For novel isoform, derived genic locus
13. For novel isoform, overlap percentage with derived genic locus
14. For novel singleton isoform, if it is located at the last exon of any known isoform. If yes, isoform ID otherwise '-'
15. For novel singleton isoform, the overlap percentage with the the last exon
16. For novel multi-exon isoform, number of splice sites are detected by anno and/or short reads; and the total number of splice sites
17. For novel multi-exon isoform, if the multi-exon isoform is the subset (based on splice junction combination) of known multi-exon isoform, isoform ID if yes otherwise '-'
18. For novel isoform, maximal length of polyA track in defined region
19. For novel isoform, maximal percentage of nucleotide A in defined region
20. Individual TSS information: TSS postion; number of supporting long reads; number of supporting CAGE position
21. Grouped TSS information: TSS postion; number of supporting long reads; number of supporting CAGE position'''

	parser = argparse.ArgumentParser(description="Function: generate ATSS-specific output file (gpd format)",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: constructed isoform (GPD file) generated by 'py_isoseqats_generate_output.py'")
	parser.add_argument('-l','--lr_gpd',type=argparse.FileType('r'),required=True,help="Input: polished long read (GPD file) generated by 'py_isoseqcon_polish.py'")
	parser.add_argument('-c','--cage_bed',type=argparse.FileType('r'),help="Optional input: CAGE data (BED format)")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: constructed isoform with TSS information (GPD file)")
	parser.add_argument('-d','--dis_group',type=int,default=5,help="If the distance (bp) between two TSSs is <= the set, group them together")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
