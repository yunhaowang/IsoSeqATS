#!/usr/bin/env python
import sys,time,re,argparse,itertools,collections

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	dic_chr_seq = parse_genome(args.genome)
	dic_chr_strand_iso_info = parse_last_exon(args.last_exon)
	generate_output(args.input,args.output,args.sgt_k_fl,args.sgt_n_fl,args.sgt_k_lr,args.sgt_n_lr,args.mlt_k_fl,args.mlt_n_fl,args.mlt_k_lr,args.mlt_n_lr,args.mlt_n_ss,args.mlt_n_sub,dic_chr_seq,args.up_len,args.down_len,args.max_a_pct,args.max_a_len,dic_chr_strand_iso_info,args.max_over_pct)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def parse_genome(genome_fa):
	dic_chr_seq = {}
	chr = ""
	seq_list = [""]
	for line in genome_fa:
		if line.startswith(">"):
			dic_chr_seq[chr] = "".join(seq_list)
			seq_list = []
			chr = line.strip().split()[0][1:]
		else:
			seq_list.append(line.strip())
	dic_chr_seq[chr] = "".join(seq_list) # add last chromosome
	del dic_chr_seq[""] # delete null key/value pair
	return dic_chr_seq
	genome_fa.close()

def determine_polya_track(dic_chr_seq,up_len,down_len,chr,strand,tss,tts):
	if strand == "+":
		def_region = dic_chr_seq[chr][(int(tts)-up_len):(int(tts)+down_len+1)]
		def_region = def_region.upper()
		# length of polyA track
		nuc_A_frq_list = []
		groups = itertools.groupby(def_region)
		for nuc_len in [(label,sum(1 for _ in group)) for label, group in groups]:
			if nuc_len[0] == "A":
				nuc_A_frq_list.append(nuc_len[1])
		if nuc_A_frq_list == []:
			nuc_A_frq_max = 0
		else:
			nuc_A_frq_max = max(nuc_A_frq_list)
		# nucleotide A percentage
		nuc_sum = collections.Counter(def_region)
		nuc_A_sum = nuc_sum["A"]
		nuc_A_pct = float(nuc_A_sum)/len(def_region)
	else:
		def_region = dic_chr_seq[chr][(int(tss)-down_len):(int(tss)+up_len+1)]
		def_region = def_region.upper()
		# length of polyT track
		nuc_A_frq_list = []
		groups = itertools.groupby(def_region)
		for nuc_len in [(label,sum(1 for _ in group)) for label, group in groups]:
			if nuc_len[0] == "T":
				nuc_A_frq_list.append(nuc_len[1])
		if nuc_A_frq_list == []:
			nuc_A_frq_max = 0
		else:
			nuc_A_frq_max = max(nuc_A_frq_list)
		# nucleotide A percentage
		nuc_sum = collections.Counter(def_region)
		nuc_A_sum = nuc_sum["T"]
		nuc_A_pct = float(nuc_A_sum)/len(def_region)

	return nuc_A_frq_max,nuc_A_pct

def parse_last_exon(last_exon_txt):
	dic_chr_strand_iso_info = {}
	for line in last_exon_txt:
		gene,iso,chr,strand,start,end = line.strip().split("\t")[:6]
		chr_strand = chr + "&" + strand
		if chr_strand not in dic_chr_strand_iso_info.keys():
			dic_chr_strand_iso_info[chr_strand] = {}
			dic_chr_strand_iso_info[chr_strand][iso] = start + "&" + end
			dic_chr_strand_iso_info[chr_strand]["iso_list"] = []
			dic_chr_strand_iso_info[chr_strand]["iso_list"].append(iso)
		else:
			dic_chr_strand_iso_info[chr_strand][iso] = start + "&" + end
			dic_chr_strand_iso_info[chr_strand]["iso_list"].append(iso)
	last_exon_txt.close()
	return dic_chr_strand_iso_info

def determine_last_exon(dic_chr_strand_iso_info,chr,strand,tss,tts):
	chr_strand = chr + "&" + strand
	if chr_strand in dic_chr_strand_iso_info.keys():
		dic_over_pct = {}
		for iso in dic_chr_strand_iso_info[chr_strand]["iso_list"]:
			if int(tts) <= int(dic_chr_strand_iso_info[chr_strand][iso].split("&")[0]): # break the loop when the TTS of novel singleton isoform is < the start site of annotated last exon
				break
			else:
				if int(tss) < int(dic_chr_strand_iso_info[chr_strand][iso].split("&")[1]) and int(tts) > int(dic_chr_strand_iso_info[chr_strand][iso].split("&")[0]): # with overlap
					over_pct = float(min(int(tts),int(dic_chr_strand_iso_info[chr_strand][iso].split("&")[1]))-max(int(tss),int(dic_chr_strand_iso_info[chr_strand][iso].split("&")[0])))/(int(tts)-int(tss))
					if over_pct not in dic_over_pct.keys():
						dic_over_pct[over_pct] = []
						dic_over_pct[over_pct].append(iso)
					else:
						dic_over_pct[over_pct].append(iso)
		if dic_over_pct != {}: # choose the maximal overlap
			iso_over_list = []
			max_pct = max(dic_over_pct.keys())
			for iso in dic_over_pct[max_pct]:
				iso_over_list.append(re.sub("^refiso_","",iso))
			iso_over = "|".join(iso_over_list)
		else:
			max_pct = 0.0
			iso_over = "-"
	else:
		max_pct = 0.0
		iso_over = "-"
	return max_pct,iso_over

def generate_output(input_gpd_list,output_gpd,sgt_k_fl,sgt_n_fl,sgt_k_lr,sgt_n_lr,mlt_k_fl,mlt_n_fl,mlt_k_lr,mlt_n_lr,mlt_n_ss,mlt_n_sub,dic_chr_seq,up_len,down_len,max_a_pct,max_a_len,dic_chr_strand_iso_info,max_over_pct):
	novel_sgt_loci_idx = 0
	novel_sgt_iso_idx = 0
	novel_mlt_loci_idx = 0
	novel_mlt_iso_idx = 0

	for input_gpd in input_gpd_list:
		for line in input_gpd:
			gene_id_set,iso_id_set,chr,strand,tss,tts,read_set,read_number,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")[:11]
			# gene id
			gene_id_list = []
			for gene in gene_id_set.split(","):
				gene_id_list.append(re.sub("^refgene_","",gene))
			gene_id = "|".join(gene_id_list)
			# isoform id
			iso_id_list = []
			for iso in iso_id_set.split(","):
				iso_id_list.append(re.sub("^refiso_","",iso))
			iso_id = "|".join(iso_id_list)
			# full-length long read
			fl_c = 0
			for lr in read_set.split(","):
				if lr.endswith("_F1P1T1"):
					fl_c += 1

			if int(exon_number) == 1: # singleton isoform
				if iso_id.startswith("novel_sgt_iso_"): # novel
					if int(read_number) >= sgt_n_lr and fl_c >= sgt_n_fl: # check read count
						nuc_A_frq_max,nuc_A_pct = determine_polya_track(dic_chr_seq,up_len,down_len,chr,strand,tss,tts)
						if nuc_A_frq_max > max_a_len or nuc_A_pct > max_a_pct: continue # check polyA track info
						last_exon_max_pct,last_exon_iso_over = determine_last_exon(dic_chr_strand_iso_info,chr,strand,tss,tts)
						if last_exon_max_pct > max_over_pct: continue # check the overlap with the last exon of known isoform
						derived_genic_loci,derived_genic_loci_over_pct = [gene_id,line.rstrip("\n").split("\t")[11]] # get derived genic locus
						novel_sgt_iso_idx += 1
						iso_id = "novel_sgt_iso_" + str(novel_sgt_iso_idx) # re-index novel singleton isoform
						if gene_id.startswith("novel_sgt_loci_"):
							novel_sgt_loci_idx += 1
							gene_id = "novel_sgt_loci_" + str(novel_sgt_loci_idx) # re-index novel singleton isoform loci
							derived_genic_loci = "-"
						print >>output_gpd, "\t".join([gene_id,iso_id,chr,strand,tss,tts,str(fl_c),read_number,exon_number,exon_start,exon_end,derived_genic_loci,derived_genic_loci_over_pct,last_exon_iso_over,str(last_exon_max_pct),"NA","NA",str(nuc_A_frq_max),str(nuc_A_pct),read_set])
				else: # known
					if int(read_number) >= sgt_k_lr and fl_c >= sgt_k_fl:
						print >>output_gpd, "\t".join([gene_id,iso_id,chr,strand,tss,tts,str(fl_c),read_number,exon_number,exon_start,exon_end,"NA","NA","NA","NA","NA","NA","NA","NA",read_set])
			else: # multi-exon isoform
				flag_ss,subset_iso_set,pct = line.rstrip("\n").split("\t")[11:14]
				# subset isoform id
				subset_iso_list = []
				for sub_iso in subset_iso_set.split(","):
					subset_iso_list.append(re.sub("^refiso_","",sub_iso))
				subset_iso = "|".join(subset_iso_list)

				if iso_id.startswith("novel_mlt_iso_"): # novel
					if int(read_number) >= mlt_n_lr and fl_c >= mlt_n_fl and (float(flag_ss.split("/")[0])/float(flag_ss.split("/")[1])) >= float(mlt_n_ss):
						if mlt_n_sub == "yes" and subset_iso != "-": continue # subset info
						nuc_A_frq_max,nuc_A_pct = determine_polya_track(dic_chr_seq,up_len,down_len,chr,strand,tss,tts)
						if nuc_A_frq_max > max_a_len or nuc_A_pct > max_a_pct: continue # polyA track info
						derived_genic_loci,derived_genic_loci_over_pct = [gene_id,pct] # derived genic locus
						novel_mlt_iso_idx += 1
						iso_id = "novel_mlt_iso_" + str(novel_mlt_iso_idx) # re-index novel multi-exon isoform
						if gene_id.startswith("novel_mlt_loci_"):
							novel_mlt_loci_idx += 1
							gene_id = "novel_mlt_loci_" + str(novel_mlt_loci_idx) # re-index novel multi-exon isoform loci
						print >>output_gpd, "\t".join([gene_id,iso_id,chr,strand,tss,tts,str(fl_c),read_number,exon_number,exon_start,exon_end,derived_genic_loci,derived_genic_loci_over_pct,"NA","NA",flag_ss,subset_iso,str(nuc_A_frq_max),str(nuc_A_pct),read_set])
				else: # known
					if int(read_number) >= mlt_k_lr and fl_c >= mlt_k_fl:
						print >>output_gpd, "\t".join([gene_id,iso_id,chr,strand,tss,tts,str(fl_c),read_number,exon_number,exon_start,exon_end,"NA","NA","NA","NA","NA","NA","NA","NA",read_set])
		input_gpd.close()
	output_gpd.close()

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
20. Long read set'''

	parser = argparse.ArgumentParser(description="Function: generate final output file (gpd format)",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,nargs="+",help="Input: gpd file generated by 'py_isoseqcon_construct_sgt.py' and/or 'py_isoseqcon_construct_mlt.py'")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: constructed isoform (gpd file)")

	group1 = parser.add_argument_group('Singleton isoform options')
	group1.add_argument('--sgt_k_fl',type=int,default=0,help="Number of full-length long reads for known singleton isoform")
	group1.add_argument('--sgt_n_fl',type=int,default=0,help="Number of full-length long reads for novel singleton isoform")
	group1.add_argument('--sgt_k_lr',type=int,default=1,help="Number of long reads for known singleton isoform")
	group1.add_argument('--sgt_n_lr',type=int,default=1,help="Number of long reads for novel singleton isoform")

	group2 = parser.add_argument_group('Multi-exon isoform options')
	group2.add_argument('--mlt_k_fl',type=int,default=0,help="Number of full-length long reads for known multi-exon isoform")
	group2.add_argument('--mlt_n_fl',type=int,default=0,help="Number of full-length long reads for novel multi-exon isoform")
	group2.add_argument('--mlt_k_lr',type=int,default=1,help="Number of long reads for known multi-exon isoform")
	group2.add_argument('--mlt_n_lr',type=int,default=1,help="Number of long reads for novel multi-exon isoform")
	group2.add_argument('--mlt_n_ss',type=float,default=0.5,help="For novel multi-exon isoform, percentage of splice sites are annotated by known annotation library. [0.0~1.0]")
	group2.add_argument('--mlt_n_sub',type=str,default="yes",choices=["yes","no"],help="For novel multi-exon isoform, exclude it if it is the subset of known isoforms")

	group3 = parser.add_argument_group('PolyA track detection (for novel isoform) options')
	group3.add_argument('-g','--genome',type=argparse.FileType('r'),required=True,help="Reference genome file (fasta format) to identify those fake novel isoforms that are casued by polyA track located intermediate region of transcripts")
	group3.add_argument('--up_len',type=int,default=5,help="Upstream length of the transcription end site for defined region (bp)")
	group3.add_argument('--down_len',type=int,default=20,help="Downstream length of the transcription end site for defined region (bp)")
	group3.add_argument('--max_a_pct',type=float,default=0.5,help="Maximal percentage of nucleotide A in defined region surrounding the transcription end site [0.0~1.0]")
	group3.add_argument('--max_a_len',type=int,default=8,help="Maximal length of polyA track in defined region surrounding the transcription end site (bp)")

	group4 = parser.add_argument_group('Last exon residual detection (for novel singleton isoform) options')
	group4.add_argument('--last_exon',type=argparse.FileType('r'),required=True,help="Last exon file (txt format), should be sorted by 'sort -k3,3 -k4,4 -k5,5n -k6,6n'")
	group4.add_argument('--max_over_pct',type=float,default=0.5,help="Maximal percentage of the overlap between novel singlton isoform and annotated last exon. [0.0~1.0]")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
