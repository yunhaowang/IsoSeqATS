#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	concatenate(args.input,args.output,args.pct)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def concatenate(anno_gpd,concatenate_anno_gpd,overlap_pct):
	dic_sgt = {}
	chr_strand_set = set()
	for line in anno_gpd:
		gene_id,isoform_id,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")[:11]
		if int(exon_number) == 1: # singleton isoform/reads
			chr_strand = chr +"$" + strand
			if chr_strand not in chr_strand_set: 
				if chr_strand_set != set():
					print >>concatenate_anno_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s,\t%s," % (",".join(dic_sgt["gene_id"]),",".join(dic_sgt["isoform_id"]),dic_sgt["chr"],dic_sgt["strand"],dic_sgt["tss"],dic_sgt["tts"],str(len(dic_sgt["gene_id"])),str(len(dic_sgt["isoform_id"])),"1",dic_sgt["tss"],dic_sgt["tts"])
				chr_strand_set.add(chr_strand)
				dic_sgt["gene_id"] = set()
				dic_sgt["isoform_id"] = set()
				dic_sgt["chr"] = chr
				dic_sgt["strand"] = strand
				dic_sgt["tss"] = tss
				dic_sgt["tts"] = tts
				dic_sgt["gene_id"].add(gene_id)
				dic_sgt["isoform_id"].add(isoform_id)
			else:
				if int(tss) < int(dic_sgt["tts"]) and ((float(min(int(dic_sgt["tts"]),int(tts)))-int(tss))/(int(tts)-int(tss)) >= overlap_pct or float((min(int(dic_sgt["tts"]),int(tts)))-int(tss))/(int(dic_sgt["tts"])-int(dic_sgt["tss"])) >= overlap_pct): # overlap percentage
					dic_sgt["tts"] = str(max(int(tts),int(dic_sgt["tts"])))
					dic_sgt["gene_id"].add(gene_id)
					dic_sgt["isoform_id"].add(isoform_id)
				else:
					print >>concatenate_anno_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s,\t%s," % (",".join(dic_sgt["gene_id"]),",".join(dic_sgt["isoform_id"]),dic_sgt["chr"],dic_sgt["strand"],dic_sgt["tss"],dic_sgt["tts"],str(len(dic_sgt["gene_id"])),str(len(dic_sgt["isoform_id"])),"1",dic_sgt["tss"],dic_sgt["tts"])
					dic_sgt["gene_id"] = set()
					dic_sgt["isoform_id"] = set()
					dic_sgt["chr"] = chr
					dic_sgt["strand"] = strand
					dic_sgt["tss"] = tss
					dic_sgt["tts"] = tts
					dic_sgt["gene_id"].add(gene_id)
					dic_sgt["isoform_id"].add(isoform_id)
	if chr_strand_set != set():
		print >>concatenate_anno_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s,\t%s," % (",".join(dic_sgt["gene_id"]),",".join(dic_sgt["isoform_id"]),dic_sgt["chr"],dic_sgt["strand"],dic_sgt["tss"],dic_sgt["tts"],str(len(dic_sgt["gene_id"])),str(len(dic_sgt["isoform_id"])),"1",dic_sgt["tss"],dic_sgt["tts"])
	anno_gpd.close()
	concatenate_anno_gpd.close()

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS
6. TTS
7. gene count
8. isoform count
9. exon count
10. exon start set
11. exon end set'''
	parser = argparse.ArgumentParser(description="Function: concatenate singleton isoforms and/or reads if they have the same strand and the overlap at setted percentage. Note: gpd file must be sorted by chromosome,strand,TSS and TTS using the command 'sort -k3,3 -k4,4 -k5,5n -k6,6n'",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: concatenated gpd file")
	parser.add_argument('-p','--pct',type=float,default=0.8,help="Percentage of overlap between two reads [0.0~1.0]")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
