#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	construction(args.anno,args.input,args.output,args.lr_known,args.lr_novel)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

#=== extract gene region from known gene annotation library ===
def extract_gene_region_from_annotation(anno_gpd):
	dic_chr_iso_info = {}
	for line in anno_gpd:
		gene,iso,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")[:11]
		if chr not in dic_chr_iso_info.keys():
			dic_chr_iso_info[chr] = {}
			dic_chr_iso_info[chr][iso] = gene + "&" + tss + "&" + tts
			dic_chr_iso_info[chr]["iso_list"] = []
			dic_chr_iso_info[chr]["iso_list"].append(iso)
		else:
			dic_chr_iso_info[chr][iso] = gene + "&" + tss + "&" + tts
			dic_chr_iso_info[chr]["iso_list"].append(iso)
	anno_gpd.close()
	return dic_chr_iso_info

def construction(anno_gpd,concat_gpd,constructed_gpd,lr_count_known,lr_count_novel): # construct isoform
	dic_chr_iso_info = extract_gene_region_from_annotation(anno_gpd)
	novel_sgt_iso_index = 0
	novel_sgt_loci_index = 0
	for line in concat_gpd:
		gene_set,iso_set,chr,strand,tss,tts,gene_number,iso_number,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")[:11]
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
				print >>constructed_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t1.0" % (",".join(refgene_set),",".join(refiso_set),chr,strand,tss,tts,",".join(lr_set),str(len(lr_set)),exon_number,exon_start,exon_end)
		else: # novel isoform
			if int(gene_number) >= lr_count_novel: # number of supported long reads
				novel_sgt_iso_index += 1
				novel_sgt_iso_name = "novel_sgt_iso_" + str(novel_sgt_iso_index)
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
					if dic_over_pct != {}: # choose the maximal overlap
						gene_over_set = set()
						max_pct = max(dic_over_pct.keys())
						for iso in dic_over_pct[max_pct]:
							gene_over_set.add(dic_chr_iso_info[chr][iso].split("&")[0])
						novel_sgt_gene_name = ",".join(list(gene_over_set))
					else:
						max_pct = "0.0"
						novel_sgt_loci_index += 1
						novel_sgt_gene_name = "novel_sgt_loci_" + str(novel_sgt_loci_index)
				else:
					max_pct = "0.0"
					novel_sgt_loci_index += 1
					novel_sgt_gene_name = "novel_sgt_loci_" + str(novel_sgt_loci_index)
				print >>constructed_gpd, "\t".join([novel_sgt_gene_name,novel_sgt_iso_name,chr,strand,tss,tts,gene_set,gene_number,exon_number,exon_start,exon_end,str(max_pct)])
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
12. overlap percentage of constructed isoform over known gene/isoform '''

	parser = argparse.ArgumentParser(description="Function: construct singleton isoform",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-a','--anno',type=argparse.FileType('r'),required=True,help="Input: annotation library (gpd file), should be sorted by 'sort -k3,3 -k5,5n -k6,6n'")
	parser.add_argument('--lr_known',type=int,default=1,help="number of long reads for supporting a known isoform")
	parser.add_argument('--lr_novel',type=int,default=1,help="number of long reads for supporting a novel isoform")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: concatenated singleton isoform (gpd file)")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: constructed isoform (gpd file)")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
