#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	convert_gpd2gtf(args.input,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def convert_gpd2gtf(input_gpd,output_gtf):
	for line in input_gpd:
		gene_id,isoform_id,chr,strand,tss,tts,full_length_LR_count,LR_count,exon_number,exon_start,exon_end,derived_genic_loci,derived_genic_loci_overlap_pct,sgt_last_exon_iso,sgt_last_exon_overlap_pct,mlt_splice_site,mlt_subset_iso,polya_len,polya_pct = line.strip().split("\t")[:19]
		grouped_tss,tss_count = line.strip().split("\t")[-2:]
		info_transcript = 'gene_id "' + gene_id + '"; transcript_id "' + isoform_id + '"; full_length_LR_count "' + full_length_LR_count + '"; LR_count "' + LR_count + '"; derived_genic_loci "' + derived_genic_loci + '"; derived_genic_loci_overlap_pct "' + derived_genic_loci_overlap_pct + '"; sgt_last_exon_iso "' + sgt_last_exon_iso + '"; sgt_last_exon_overlap_pct "' + sgt_last_exon_overlap_pct + '"; mlt_splice_site "' + mlt_splice_site + '"; mlt_subset_iso "' + mlt_subset_iso + '"; polya_len "' + polya_len + '"; polya_pct "' + polya_pct + '"; grouped_tss "' + grouped_tss + '"; tss_count ' + tss_count + '";'
		print >>output_gtf, chr + '\t.\ttranscript\t' + str(int(tss)+1) + '\t' + str(tts) + '\t.\t' + strand + '\t.\t' + info_transcript
		if strand == "+":
			for i in range(0,int(exon_number)):
				info = info_transcript + ' exon_number "' + str(i+1) + '";'
				print >>output_gtf, chr + "\t.\texon\t" + str(int(exon_start.split(",")[i])+1) + "\t" + exon_end.split(",")[i] + "\t.\t" + strand + "\t.\t" + info
		else:
			for i in range(0,int(exon_number)):
				info = info_transcript + '"; exon_number "' + str(int(exon_number)-i) + '";'
				print >>output_gtf, chr + '\t.\texon\t' + str(int(exon_start.split(",")[i])+1) + '\t' + exon_end.split(",")[i] + '\t.\t' + strand + '\t.\t' + info

	input_gpd.close()
	output_gtf.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Convert gpd to gtf",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gtf file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
