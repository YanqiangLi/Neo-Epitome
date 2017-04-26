import argparse,os,sys,time,logging,datetime

def STAR_alignment(readsFilesCaseRNA, starGenomeDir,outDir):
	gz=readsFilesCaseRNA.split(',')[0].endswith('gz')
	readsFiles_split=' '.join(readsFilesCaseRNA.split(','))
	if gz:
		cmd1='STAR --genomeDir '+starGenomeDir+' --twopassMode Basic --readFilesIn '+readsFiles_split+' --runThreadN 6 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --readFilesCommand zcat --outFileNamePrefix '+outDir+'/'
	else:
		cmd1='STAR --genomeDir '+starGenomeDir+' --twopassMode Basic --readFilesIn '+readsFiles_split+' --runThreadN 6 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outFileNamePrefix '+outDir+'/'	
	logging.debug('[RNA-seq] Running command 1: '+cmd1+'\n')
	os.system(cmd1)
	
def cufflinks(gtf, sampleID):
	cmd2='cufflinks --multi-read-correct -p 8 --GTF '+gtf+' -o '+sampleID+'/cufflinks '+sampleID+'/Aligned.sortedByCoord.out.bam'
	logging.debug('[RNA-seq] Running command 2: '+cmd2+'\n')
	os.system(cmd2)


parser = argparse.ArgumentParser(description='NeoEpitome-Fastq2Mut (v1.0)')
parser.add_argument('--starGenomeDir',help='Reference genome used for short reads mapping and GATK realigments.')
parser.add_argument('--gtf',help='')
parser.add_argument('-p','--sampleID', default='NeoEpitomeOut', help='Sample ID will be used for output folder name and reads group name.')
parser.add_argument('readsFilesCaseRNA',help='Tumor sample paired-end fastq files seperated by ",". ')

args = parser.parse_args()

sampleID=args.sampleID.rstrip('/')
os.system('mkdir -p '+args.sampleID)
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=sampleID+'/NeoEpitome-Fastq2Exp.log'+ str(datetime.datetime.now())+'.txt' ,
                    filemode='w')

logging.debug('[RNA-seq] # Start STAR 2pass alignment.')
print '[RNA-seq] # Start STAR 2pass alignment.'
if os.path.exists(sampleID+'/Aligned.sortedByCoord.out.bam')==False:
	STAR_alignment(args.readsFilesCaseRNA, args.starGenomeDir, sampleID)
	if os.path.exists(sampleID+'/Aligned.sortedByCoord.out.bam')==False:
		sys.exit('[RNA-seq] # An Error Occured. cufflinks Incomplete. Exit!')
else:
	logging.debug('[RNA-seq] # Skipped STAR 2pass alignment.')
	print '[RNA-seq] # Skipped STAR 2pass alignment.'

logging.debug('[RNA-seq] # Start cufflinks.')
print '[RNA-seq] # Start cufflinks.'
if os.path.exists(sampleID+'/cufflinks/genes.fpkm_tracking')==False:
	cufflinks(args.gtf,sampleID)
	if os.path.exists(sampleID+'/cufflinks/genes.fpkm_tracking')==False:
		sys.exit('[RNA-seq] # An Error Occured. cufflinks Incomplete. Exit!')
else:
	logging.debug('[RNA-seq] # Skipped cufflinks.')
	print '[RNA-seq] # Skipped cufflinks.'
logging.debug('[RNA-seq] # Completed.')
print '[RNA-seq] # Completed.'
