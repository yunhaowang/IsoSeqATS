#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	dic_chr_spliceSite_anno,dic_junction_set,dic_chr_iso_info = extract_info_from_annotation(args.anno)
	if args.short_reads:
		dic_chr_spliceSite_sr = extract_spliceSite_from_short_reads(args.short_reads)
		dic_chr_spliceSite = {}
		chr_set = set(dic_chr_spliceSite_anno.keys()+dic_chr_spliceSite_sr.keys())
		for chr in chr_set:
			if chr in dic_chr_spliceSite_anno.keys() and chr in dic_chr_spliceSite_sr.keys():
				dic_chr_spliceSite[chr] = dic_chr_spliceSite_anno[chr] + dic_chr_spliceSite_sr[chr]
			elif chr in dic_chr_spliceSite_anno.keys() and chr not in dic_chr_spliceSite_sr.keys():
				dic_chr_spliceSite[chr] = dic_chr_spliceSite_anno[chr]
			else:
				dic_chr_spliceSite[chr] = dic_chr_spliceSite_sr[chr]
	else:
		dic_chr_spliceSite = dic_chr_spliceSite_anno

	construction(dic_chr_spliceSite,dic_junction_set,dic_chr_iso_info,args.input,args.output,args.lr_known,args.lr_novel)

	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

#=== extract splice site from short reads ===
def extract_spliceSite_from_short_reads(sr_gpd):
	dic_chr_spliceSite = {}
	for line in sr_gpd:
		gene_id,isoform_id,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")[:11]
		if int(exon_number) > 1:
			if chr not in dic_chr_spliceSite.keys():
				dic_chr_spliceSite[chr] = exon_start.split(",")[1:-1] + exon_end.split(",")[:-2]
			else:
				dic_chr_spliceSite[chr] += (exon_start.split(",")[1:-1] + exon_end.split(",")[:-2])
				dic_chr_spliceSite[chr] = list(set(dic_chr_spliceSite[chr]))
	return dic_chr_spliceSite 
	sr_gpd.close()

#=== extract splice site, junction set, and gene region informations from known gene annotation library ===
def extract_info_from_annotation(anno_gpd):
	dic_chr_spliceSite = {}
	dic_junction_set = {}
	dic_chr_iso_info = {}
	for line in anno_gpd:
		gene,iso,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")[:11]
		if int(exon_number) > 1:

			chr_strand = chr + "_" + strand # junction set
			if chr_strand not in dic_junction_set.keys():
				dic_junction_set[chr_strand] = {}
				dic_junction_set[chr_strand]["-".join(exon_start.split(",")[1:-1])+"_"+"-".join(exon_end.split(",")[:-2])] = []
				dic_junction_set[chr_strand]["-".join(exon_start.split(",")[1:-1])+"_"+"-".join(exon_end.split(",")[:-2])].append(iso)

			else:
				if ("-".join(exon_start.split(",")[1:-1])+"_"+"-".join(exon_end.split(",")[:-2])) not in dic_junction_set[chr_strand].keys():
					dic_junction_set[chr_strand]["-".join(exon_start.split(",")[1:-1])+"_"+"-".join(exon_end.split(",")[:-2])] = []
					dic_junction_set[chr_strand]["-".join(exon_start.split(",")[1:-1])+"_"+"-".join(exon_end.split(",")[:-2])].append(iso)
				else:
					dic_junction_set[chr_strand]["-".join(exon_start.split(",")[1:-1])+"_"+"-".join(exon_end.split(",")[:-2])].append(iso)

			if chr not in dic_chr_spliceSite.keys(): # splice site
				dic_chr_spliceSite[chr] = exon_start.split(",")[1:-1] + exon_end.split(",")[:-2]
			else:
				dic_chr_spliceSite[chr] += (exon_start.split(",")[1:-1] + exon_end.split(",")[:-2])
				dic_chr_spliceSite[chr] = list(set(dic_chr_spliceSite[chr]))

		if chr not in dic_chr_iso_info.keys(): # gene region
			dic_chr_iso_info[chr] = {}
			dic_chr_iso_info[chr][iso] = gene + "&" + tss + "&" + tts
			dic_chr_iso_info[chr]["iso_list"] = []
			dic_chr_iso_info[chr]["iso_list"].append(iso)
		else:
			dic_chr_iso_info[chr][iso] = gene + "&" + tss + "&" + tts
			dic_chr_iso_info[chr]["iso_list"].append(iso)

	return dic_chr_spliceSite,dic_junction_set,dic_chr_iso_info
	anno_gpd.close()

#=== test splice sites of long reads are annotated by annotation library and/or detected by short reads ===
def test_splice_sites(lr_mlt_gpd,dic_chr_spliceSite):
	lr_id,lr_name,chr,strand,tss,tts,mapq,sf,exon_number,exon_start,exon_end = lr_mlt_gpd.split("\t")[:11]
	count_yes = 0
	if chr in dic_chr_spliceSite.keys():
		for ss in (exon_start.split(",")[1:-1]+exon_end.split(",")[:-2]):
			if ss in dic_chr_spliceSite[chr]:
				count_yes += 1
	flag_ss = str(count_yes) + "/" + str(2+(int(exon_number)-2)*2)
	return flag_ss

#==== construction ===
def construction(dic_chr_spliceSite,dic_junction_set,dic_chr_iso_info,concat_gpd,constructed_gpd,lr_count_known,lr_count_novel): # construct isoform
	novel_mlt_iso_index = 0
	novel_mlt_loci_index = 0
	for line in concat_gpd:
		gene_set,iso_set,chr,strand,tss,tts,gene_number,iso_number,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")[:11]
		if int(exon_number) == 1:
			continue
		if "refgene_" in gene_set: # known isoform
			refgene_set = set()
			refiso_set = set()
			lr_set = set()
			for gene_id in gene_set.split(","):
				if gene_id.startswith("refgene_"):
					refgene_set.add(gene_id)
				else:
					lr_set.add(gene_id)
			if len(lr_set) >= lr_count_known: # number of supported long reads
				for iso_id in iso_set.split(","):
					if iso_id.startswith("refiso_"):
						refiso_set.add(iso_id)
				refgene_set = list(refgene_set)
				refiso_set = list(refiso_set)
				refgene_set.sort()
				refiso_set.sort()
				flag_ss = str(2+(int(exon_number)-2)*2) + "/" + str(2+(int(exon_number)-2)*2)
				print >>constructed_gpd, "\t".join([",".join(refgene_set),",".join(refiso_set),chr,strand,tss,tts,",".join(lr_set),str(len(lr_set)),exon_number,exon_start,exon_end,flag_ss,",".join(refiso_set),"1.0"])
		else: # novel isoform
			if int(gene_number) < lr_count_novel: # number of supported long reads 
				continue
			novel_mlt_iso_index += 1
			novel_mlt_iso_name = "novel_mlt_iso_" + str(novel_mlt_iso_index)
			lr_gpd = line.rstrip("\n")
			flag_ss = test_splice_sites(lr_gpd,dic_chr_spliceSite)
			p5 = "-".join(exon_start.split(",")[1:-1])
			p3 = "-".join(exon_end.split(",")[:-2])
			chr_strand = chr+"_"+strand
			if chr_strand in dic_junction_set.keys():
				for p5_p3 in dic_junction_set[chr_strand].keys():
					if p5 in p5_p3.split("_")[0] and p3 in p5_p3.split("_")[1]: # if the junction combination set is a subset of any known isoform
						iso_subset = ",".join(dic_junction_set[chr_strand][p5_p3])
						break
					else:
						iso_subset = ""
			else:
				iso_subset = ""
				
			if iso_subset != "":
				novel_mlt_gene_name = dic_chr_iso_info[chr][iso_subset.split(",")[0]].split("&")[0]
				print >>constructed_gpd, "\t".join([novel_mlt_gene_name,novel_mlt_iso_name,chr,strand,tss,tts,gene_set,gene_number,exon_number,exon_start,exon_end,flag_ss,iso_subset,"1.0"])

			else:
				if chr in dic_chr_iso_info.keys():
					dic_over_pct = {}
					for iso in dic_chr_iso_info[chr]["iso_list"]:
						if int(tts) <= int(dic_chr_iso_info[chr][iso].split("&")[1]): # break the loop when the TTS of reads is < the TSS of annotation
							break
						else:
							if int(tss) < int(dic_chr_iso_info[chr][iso].split("&")[2]) and int(tts) > int(dic_chr_iso_info[chr][iso].split("&")[1]): # with overlap
								over_pct = float(min(int(tts),int(dic_chr_iso_info[chr][iso].split("&")[2]))-max(int(tss),int(dic_chr_iso_info[chr][iso].split("&")[1])))/(int(tts)-int(tss))
								if over_pct not in dic_over_pct.keys():
									dic_over_pct[over_pct] = []
									dic_over_pct[over_pct].append(iso)
								else:
									dic_over_pct[over_pct].append(iso)
					if dic_over_pct != {}: # choose the maximumal overlap
						gene_over_set = set()
						max_pct = max(dic_over_pct.keys())
						for iso in dic_over_pct[max_pct]:
							gene_over_set.add(dic_chr_iso_info[chr][iso].split("&")[0])
						novel_mlt_gene_name = ",".join(list(gene_over_set))
					else:
						max_pct = "0.0"
						novel_mlt_loci_index += 1
						novel_mlt_gene_name = "novel_mlt_loci_" + str(novel_mlt_loci_index)
				else:
					max_pct = "0.0"
					novel_mlt_loci_index += 1
					novel_mlt_gene_name = "novel_mlt_loci_" + str(novel_mlt_loci_index)

				print >>constructed_gpd, "\t".join([novel_mlt_gene_name,novel_mlt_iso_name,chr,strand,tss,tts,gene_set,gene_number,exon_number,exon_start,exon_end,flag_ss,"-",str(max_pct)])

	concat_gpd.close()
	constructed_gpd.close()

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS
6. TTS
7. long read set
8. long read count
9. exon count
10. exon start set
11. exon end set
12. proportion of splice sites are detected by anno and/or short reads
13. if the isoform is the subset of known isoform, isoform ID if yes otherwise '-'
14. overlap percentage of constructed isoform over known gene/isoform '''

	parser = argparse.ArgumentParser(description="Function: construct multi-exon isoform",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-a','--anno',type=argparse.FileType('r'),required=True,help="Input: annotation library (gpd file), should be sorted by 'sort -k3,3 -k5,5n -k6,6n'")
	parser.add_argument('-s','--short_reads',type=argparse.FileType('r'),help="Optional input: short reads for confirming splice sites (gpd file)")
	parser.add_argument('--lr_known',type=int,default=1,help="number of long reads for supporting a known isoform")
	parser.add_argument('--lr_novel',type=int,default=1,help="number of long reads for supporting a novel isoform")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: concatenated singleton isoform (gpd file)")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: constructed isoform (gpd file)")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
