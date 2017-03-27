import sys, os, argparse 

def create_position_list(fin, outdir):
	fout=open(outdir+'/somatic.positions.txt','w')
	for line in open(fin):
		if not line.startswith("#"):
			ls=line.strip().split("\t")
			pos=(ls[0],ls[1],str(int(ls[1])))
			fout.write('\t'.join(pos)+'\n')

def main():
	parser = argparse.ArgumentParser(description='NeoEpitome-Mut2Antigen (v1.0)')
    parser.add_argument('vcf_input', help='input annotated somatic mutation VCF file.')
    parser.add_argument('-o','--outdir', help='')
    parser.add_argument('-r','--reference',help='Reference genome used for short reads mapping and GATK realigments.')
    parser.add_argument('--tumor-dna-bam', type=argparse.FileType('r'), help='bam-readcount output file for tumor DNA BAM and snvs')
    parser.add_argument('--tumor-rna-bam', type=argparse.FileType('r'), help='bam-readcount output file for tumor RNA BAM and snvs')
	args = parser.parse_args()

	fin=args.vcf_input
	outdir=args.outdir

	create_position_list(fin,outdir)
	cmd1='bam-readcount -f '+args.reference+' -w 1 -l '+outdir+'/somatic.positions.txt '+args.tumor_dna_bam+'> '+outdir+'/DNA.readcount'
	os.system(cmd1)
	cmd2='bam-readcount -f '+args.reference+' -w 1 -l '+outdir+'/somatic.positions.txt '+args.tumor_rna_bam+'> '+outdir+'/RNA.readcount'
	os.system(cmd2)

if __name__ == '__main__':
	main()
