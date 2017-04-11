import argparse,os,sys,time,logging,datetime
import multiprocessing as mp

startTime = time.time()
JAVA7='/usr/bin/java'
JAVA8='/u/local/apps/java/jdk1.8.0_111/bin/java'


def SNP_calling_MuSE(outPath, sampleID, reference, known_snps):
	logging.debug('[DNA-seq] # Start MuSE.')
	print '[DNA-seq] # Start MuSE.'
	if os.path.exists(outPath+sampleID+'/'+sampleID+'.muse.snp.vcf')==False:

		cmd9='MuSE call -f '+reference+' '+outPath+sampleID+'.case/case.brsq.idrealn.mkdup.merged.sorted.output.bam '+outPath+sampleID+'.ctrl/ctrl.brsq.idrealn.mkdup.merged.sorted.output.bam '+'-O '+outPath+sampleID+'/'+sampleID+'.intermediate'
		logging.debug('[DNA-seq] Running command 9: '+cmd9+'\n')
		os.system(cmd9)
		cmd10='MuSE sump -E -I '+outPath+sampleID+'/'+sampleID+'.intermediate.MuSE.txt'+' -O '+outPath+sampleID+'/'+sampleID+'.muse.snp.vcf -D '+known_snps
		logging.debug('[DNA-seq] Running command 10: '+cmd10+'\n')
		os.system(cmd10)

	if os.path.exists(outPath+sampleID+'/'+sampleID+'.muse.snp.vcf')==False:
		sys.exit('[DNA-seq] # An Error Occured. MuSE Incomplete. Exit!')
	logging.debug('[DNA-seq] # MuSE SNV calling completed.')
	print '[DNA-seq] # MuSE SNV calling completed.'

def SNP_calling_MuTect_command(outPath, sampleID, reference, known_snps, jarPath, n):

	cmd9=JAVA7+' -jar '+jarPath+'/mutect.jar --analysis_type MuTect --intervals '+outPath+sampleID+'/'+sampleID+'.intervals.'+str(n)+'.bed --reference_sequence '+reference+' --dbsnp '+known_snps+' --input_file:normal '+outPath+sampleID+'.ctrl/ctrl.brsq.idrealn.mkdup.merged.sorted.output.bam --input_file:tumor '+outPath+sampleID+'.case/case.brsq.idrealn.mkdup.merged.sorted.output.bam --out '+outPath+sampleID+'/'+sampleID+'.'+str(n)+'.mutect.snp.txt'
	logging.debug('[DNA-seq] Running command 9: '+cmd9+'\n')
	os.system(cmd9)

def SNP_calling_MuTect(outPath, sampleID, reference, known_snps, jarPath):
	logging.debug('[DNA-seq] # Start MuTect.')
	print '[DNA-seq] # Start MuTect.'

	if os.path.exists(outPath+sampleID+'/'+sampleID+'.mutect.snp.txt')==False:
		base=0
		m=1
		fout=open(outPath+sampleID+'/'+sampleID+'.intervals.'+str(m)+'.bed','w')
		for l in open(reference+'.fai'):
			ls=l.strip().split('\t')
			base+=int(ls[1])
			if base<=300000000:
				fout.write('{}\t{}\t{}\n'.format(ls[0],0,int(ls[1])-1))
			else:
				base=0
				fout.close()
				m+=1
				fout=open(outPath+sampleID+'/'+sampleID+'.intervals.'+str(m)+'.bed','w')
				fout.write('{}\t{}\t{}\n'.format(ls[0],0,int(ls[1])-1))
		fout.close()

		mutect_calling=[]
		for n in xrange(1,m+1):
			mutect_calling.append(mp.Process(target=SNP_calling_MuTect_command,args=(outPath, sampleID, reference, known_snps, jarPath, n)))
			mutect_calling[n-1].start()
		for n in xrange(1,m+1):
			mutect_calling[n-1].join()
		
		os.system('cat '+outPath+sampleID+'/'+sampleID+'.*.mutect.snp.txt >'+outPath+sampleID+'/'+sampleID+'.mutect.snp.txt')

	if os.path.exists(outPath+sampleID+'/'+sampleID+'.mutect.snp.txt')==False:
		sys.exit('[DNA-seq] # An Error Occured. MuTect Incomplete. Exit!')
	logging.debug('[DNA-seq] # MuTect SNV calling completed.')
	print '[DNA-seq] # MuTect SNV calling completed.'
	
def SNP_calling_VarScan2(outPath, sampleID, reference, jarPath):
	logging.debug('[DNA-seq] # Start VarScan2.')
	print '[DNA-seq] # Start VarScan2.'
	if os.path.exists(outPath+sampleID+'/'+sampleID+'.varscan.snp.vcf')==False:

		cmd9='samtools mpileup -f '+reference+' -q 1 -B '+outPath+sampleID+'.ctrl/ctrl.brsq.idrealn.mkdup.merged.sorted.output.bam '+outPath+sampleID+'.case/case.brsq.idrealn.mkdup.merged.sorted.output.bam > '+outPath+sampleID+'/intermediate.pileup'
		logging.debug('[DNA-seq] Running command 9: '+cmd9+'\n')
		os.system(cmd9)

		cmd10=JAVA8+' -jar '+jarPath+'/VarScan.jar somatic '+outPath+sampleID+'/intermediate.pileup '+outPath+sampleID+'/'+sampleID+'.varscan --mpileup 1 --min-coverage 8 --min-coverage-normal 8 --min-coverage-tumor 6 --min-var-freq 0.10 --min-freq-for-hom 0.75 --normal-purity 1.0 --tumor-purity 1.00 --p-value 0.99 --somatic-p-value 0.05 --strand-filter 0 --output-vcf'
		logging.debug('[DNA-seq] Running command 10: '+cmd10+'\n')
		os.system(cmd10)

		cmd11=JAVA8+' -jar '+jarPath+'/VarScan.jar processSomatic '+outPath+sampleID+'/'+sampleID+'.varscan.snp.vcf --min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07'
		logging.debug('[DNA-seq] Running command 11: '+cmd11+'\n')
		os.system(cmd11)

		cmd12=JAVA8+' -jar '+jarPath+'/VarScan.jar processSomatic '+outPath+sampleID+'/'+sampleID+'.varscan.indel.vcf --min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07'
		logging.debug('[DNA-seq] Running command 12: '+cmd12+'\n')
		os.system(cmd12)

	if os.path.exists(outPath+sampleID+'/'+sampleID+'.varscan.snp.vcf')==False:
		sys.exit('[DNA-seq] # An Error Occured. VarScan Incomplete. Exit!')
	os.system('rm '+outPath+sampleID+'/intermediate.pileup')
	logging.debug('[DNA-seq] # VarScan SNV calling completed.')
	print '[DNA-seq] # VarScan SNV calling completed.'

def BRSQ(reference,binDir,prefixOut,known_snps):
	
	cmd7=JAVA8+' -jar '+binDir+'/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 6 -R '+reference+' -I '+prefixOut+'.idrealn.mkdup.merged.sorted.output.bam -knownSites '+known_snps+' -o '+prefixOut+'.bqsr.grp'
	logging.debug('[DNA-seq] Running command 7: '+cmd7+'\n')
	os.system(cmd7)

	cmd8=JAVA8+' -jar '+binDir+'/GenomeAnalysisTK.jar -T PrintReads -nct 6 -R '+reference+' -I '+prefixOut+'.idrealn.mkdup.merged.sorted.output.bam --BQSR '+prefixOut+'.bqsr.grp -o '+prefixOut+'.brsq.idrealn.mkdup.merged.sorted.output.bam'
 	logging.debug('[DNA-seq] Running command 8: '+cmd8+'\n')
 	os.system(cmd8)

def DNA_mapping(readsFiles,binDir,prefixOut,sampleID,reference):

	readsFiles_split=' '.join(readsFiles.split(','))
	folderPath=prefixOut.rsplit('/',1)[0]
	os.system('mkdir -p '+folderPath)

	cmd1='bwa mem -v 1 -t 8 -T 0 -R \'@RG\\tID:'+sampleID+'\\tSM:'+sampleID+'\\tPL:ILLUMINA\\tPU:lane1\\tLB:'+sampleID+'\' '+reference+' '+readsFiles_split+' |samtools view -Shb -o '+prefixOut+'.raw.output.bam - '
	logging.debug('[DNA-seq] Running command 1: '+cmd1+'\n')
	os.system(cmd1)

	cmd2=JAVA8+' -jar '+binDir+'/picard.jar SortSam QUIET=true CREATE_INDEX=true INPUT='+prefixOut+'.raw.output.bam OUTPUT='+prefixOut+'.sorted.output.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT'
	logging.debug('[DNA-seq] Running command 2: '+cmd2+'\n')
	os.system(cmd2)

	
	cmd3=JAVA8+' -jar '+binDir+'/picard.jar MergeSamFiles QUIET=true ASSUME_SORTED=false CREATE_INDEX=true INPUT='+prefixOut+'.sorted.output.bam MERGE_SEQUENCE_DICTIONARIES=false OUTPUT='+prefixOut+'.merged.sorted.output.bam SORT_ORDER=coordinate USE_THREADING=true VALIDATION_STRINGENCY=STRICT'
	logging.debug('[DNA-seq] Running command 3: '+cmd3+'\n')
	os.system(cmd3)

	cmd4=JAVA8+' -jar '+binDir+'/picard.jar MarkDuplicates QUIET=true CREATE_INDEX=true VALIDATION_STRINGENCY=STRICT INPUT='+prefixOut+'.merged.sorted.output.bam O='+prefixOut+'.mkdup.merged.sorted.output.bam M='+prefixOut+'.mkdup.merged.sorted.output.txt'
	logging.debug('[DNA-seq] Running command 4: '+cmd4+'\n')
	os.system(cmd4)

parser = argparse.ArgumentParser(description='NeoEpitome-Fastq2Mut (v1.0)')
parser.add_argument('-r','--reference',help='Reference genome used for short reads mapping and GATK realigments.')
parser.add_argument('-i','--known-indels',help='known indel for GATK use. Sorted VCF format. Detail see GATK.')
parser.add_argument('-s','--known-snps',help='dbSNP for GATK use. Sorted VCF format. Detail see GATK.')
parser.add_argument('readsFilesCase',help='Tumor sample paired-end fastq files seperated by ",". ')
parser.add_argument('readsFilesCtrl', help='Tumor sample paired-end fastq files seperated by ",".')
parser.add_argument('-d','--binDir', help='Directory for java applications indluding GATK, picard.')
parser.add_argument('-p','--sampleID', default='NeoEpitomeOut', help='Sample ID will be used for output folder name and reads group name.')
parser.add_argument('--snv-calling-method',default='VarScan2', help='SNV calling method can be set as VarScan2, MuTect, MuSE, or concensus of the three. Default is VarScan2.')
args = parser.parse_args()
if args.sampleID.rstrip('/').find('/')!=-1:
	[outPath,sampleID]=args.sampleID.rstrip('/').rsplit('/',1)
	outPath+='/'
else:
	sampleID=args.sampleID.rstrip('/')
	outPath=''


jarPath='/'+args.binDir.strip('/')

os.system('mkdir -p '+outPath+sampleID)

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=outPath+sampleID+'/NeoEpitome-Fastq2Mut.log'+ str(datetime.datetime.now())+'.txt' ,
                    filemode='w')
#Step 1
logging.debug('[DNA-seq] # Start DNA mapping.')
print '[DNA-seq] # Start DNA mapping.'
if os.path.exists(outPath+sampleID+'.case/case.mkdup.merged.sorted.output.bam')==False and os.path.exists(outPath+sampleID+'.ctrl/ctrl.mkdup.merged.sorted.output.bam')==False:
	
	mappings=[]
	mappings.append(mp.Process(target=DNA_mapping,args=(args.readsFilesCase,jarPath,outPath+sampleID+'.case/case',sampleID, args.reference)))
	mappings[0].start()
	mappings.append(mp.Process(target=DNA_mapping,args=(args.readsFilesCtrl,jarPath,outPath+sampleID+'.ctrl/ctrl',sampleID, args.reference)))
	mappings[1].start()
	mappings[0].join()
	mappings[1].join()

	if os.path.exists(outPath+sampleID+'.case/case.mkdup.merged.sorted.output.bam')==False and os.path.exists(outPath+sampleID+'.ctrl/ctrl.mkdup.merged.sorted.output.bam')==False:
		sys.exit('[DNA-seq] # An Error Occured. DNA mapping Incomplete. Exit!')
	os.system('rm '+outPath+sampleID+'.case/case.merged.sorted.output.bam')
	os.system('rm '+outPath+sampleID+'.case/case.sorted.output.bam')
	os.system('rm '+outPath+sampleID+'.ctrl/ctrl.merged.sorted.output.bam')
	os.system('rm '+outPath+sampleID+'.ctrl/ctrl.sorted.output.bam')

else:
	logging.debug('[DNA-seq] # Skipped DNA mapping.')
	print '[DNA-seq] # Skipped DNA mapping.'


#Step 2
logging.debug('[DNA-seq] # GATK bam refining - Indel.')
print '[DNA-seq] # GATK bam refining - Indel.'
os.system('mkdir -p '+outPath+sampleID)
if os.path.exists(outPath+sampleID+'.case/case.idrealn.mkdup.merged.sorted.output.bam')==False and os.path.exists(outPath+sampleID+'.ctrl/ctrl.idrealn.mkdup.merged.sorted.output.bam')==False:

	cmd5=JAVA8+' -jar '+jarPath+'/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 8 -R '+args.reference+' -known '+args.known_indels+' -I '+outPath+sampleID+'.case/case.mkdup.merged.sorted.output.bam -I '+outPath+sampleID+'.ctrl/ctrl.mkdup.merged.sorted.output.bam -o '+outPath+sampleID+'/realign_target.intervals'
	logging.debug('[DNA-seq] Running command 5: '+cmd5+'\n')
	os.system(cmd5)

	os.system('rm -f '+outPath+sampleID+'/output.map')
	guide_file=open(outPath+sampleID+'/output.map','w')
	guide_file.write('case.mkdup.merged.sorted.output.bam\t'+outPath+sampleID+'.case/case.idrealn.mkdup.merged.sorted.output.bam\nctrl.mkdup.merged.sorted.output.bam\t'+outPath+sampleID+'.ctrl/ctrl.idrealn.mkdup.merged.sorted.output.bam')
	guide_file.close()

	cmd6=JAVA8+' -jar '+jarPath+'/GenomeAnalysisTK.jar -T IndelRealigner -R '+args.reference+' -known '+args.known_indels+' -targetIntervals '+outPath+sampleID+'/realign_target.intervals --noOriginalAlignmentTags -I '+outPath+sampleID+'.case/case.mkdup.merged.sorted.output.bam -I '+outPath+sampleID+'.ctrl/ctrl.mkdup.merged.sorted.output.bam -nWayOut '+outPath+sampleID+'/output.map'
	logging.debug('[DNA-seq] Running command 6: '+cmd6+'\n')
	os.system(cmd6)

	os.system('mv '+outPath+sampleID+'/output.map '+outPath+sampleID+'/output.map_USED')

	if os.path.exists(outPath+sampleID+'.case/case.idrealn.mkdup.merged.sorted.output.bam')==False and os.path.exists(outPath+sampleID+'.ctrl/ctrl.idrealn.mkdup.merged.sorted.output.bam')==False:
		sys.exit('[DNA-seq] # An Error Occured. GATK bam refining - Indel Incomplete. Exit!')
	#os.system('rm '+outPath+sampleID+'.case/case.mkdup.merged.sorted.output.bam')
	#os.system('rm '+outPath+sampleID+'.ctrl/ctrl.mkdup.merged.sorted.output.bam')
else:
	logging.debug('[DNA-seq] # Skipped GATK bam refining - Indel.')
	print '[DNA-seq] # Skipped GATK bam refining - Indel.'

#Step 3
logging.debug('[DNA-seq] # GATK bam refining - SNP.')
print '[DNA-seq] # GATK bam refining - SNP.'

if os.path.exists(outPath+sampleID+'.case/case.brsq.idrealn.mkdup.merged.sorted.output.bam')==False and os.path.exists(outPath+sampleID+'.ctrl/ctrl.brsq.idrealn.mkdup.merged.sorted.output.bam')==False:
	
	recalibratings=[]
	recalibratings.append(mp.Process(target=BRSQ,args=(args.reference,jarPath,outPath+sampleID+'.case/case',args.known_snps)))
	recalibratings[0].start()
	recalibratings.append(mp.Process(target=BRSQ,args=(args.reference,jarPath,outPath+sampleID+'.ctrl/ctrl',args.known_snps)))
	recalibratings[1].start()
	recalibratings[0].join()
	recalibratings[1].join()

	if os.path.exists(outPath+sampleID+'.case/case.brsq.idrealn.mkdup.merged.sorted.output.bam')==False and os.path.exists(outPath+sampleID+'.ctrl/ctrl.brsq.idrealn.mkdup.merged.sorted.output.bam')==False:
		sys.exit('[DNA-seq] # An Error Occured. GATK bam refining - SNP Incomplete. Exit!')
	#os.system('rm '+outPath+sampleID+'.case/case.idrealn.mkdup.merged.sorted.output.bam')
	#os.system('rm '+outPath+sampleID+'.ctrl/ctrl.idrealn.mkdup.merged.sorted.output.bam')

else:
	logging.debug('[DNA-seq] # Skipped GATK bam refining - SNP.')
	print '[DNA-seq] # Skipped GATK bam refining - SNP.'

#Step 4
if args.snv_calling_method=='VarScan2':
	SNP_calling_VarScan2(outPath, sampleID, args.reference, jarPath)
elif args.snv_calling_method=='MuTect':
	SNP_calling_MuTect(outPath, sampleID, args.reference, args.known_snps, jarPath)
elif args.snv_calling_method=='MuSE':
	SNP_calling_MuSE(outPath, sampleID, args.reference, args.known_snps)
elif args.snv_calling_method=='consensus':
	callings=[]
	callings.append(mp.Process(target=SNP_calling_VarScan2,args=(outPath, sampleID, args.reference, jarPath)))
	callings[0].start()
	callings.append(mp.Process(target=SNP_calling_MuSE,args=(outPath, sampleID, args.reference, args.known_snps)))
	callings[1].start()
	callings.append(mp.Process(target=SNP_calling_MuTect,args=(outPath, sampleID, args.reference, args.known_snps, jarPath)))
	callings[2].start()
	callings[0].join()
	callings[1].join()
	callings[2].join()
	mut_dict={}
	for l in open(outPath+sampleID+'/'+sampleID+'.varscan.snp.Somatic.hc.vcf'):
		if l.startswith('#'):
			continue
		ls=l.strip().split('\t')
		mut_dict[ls[0]+'_'+ls[1]+'_'+ls[3]+'_'+ls[4]]=set(['VarScan2'])
	for l in open(outPath+sampleID+'/'+sampleID+'.varscan.indel.Somatic.hc.vcf'):
		if l.startswith('#'):
			continue
		ls=l.strip().split('\t')
		if mut_dict.has_key(ls[0]+'_'+ls[1]+'_'+ls[3]+'_'+ls[4]):
			mut_dict[ls[0]+'_'+ls[1]+'_'+ls[3]+'_'+ls[4]].add(['VarScan2'])
		else:
			mut_dict[ls[0]+'_'+ls[1]+'_'+ls[3]+'_'+ls[4]]=set(['VarScan2'])
	for l in open(outPath+sampleID+'/'+sampleID+'.muse.snp.vcf'):
		if l.startswith('#'):
			continue
		ls=l.strip().split('\t')
		if mut_dict.has_key(ls[0]+'_'+ls[1]+'_'+ls[3]+'_'+ls[4]):
			mut_dict[ls[0]+'_'+ls[1]+'_'+ls[3]+'_'+ls[4]].add('MuSE')
		else:
			mut_dict[ls[0]+'_'+ls[1]+'_'+ls[3]+'_'+ls[4]]=set(['MuSE'])
	for l in open(outPath+sampleID+'/'+sampleID+'.mutect.snp.txt'):
		if l.startswith('#'):
			continue
		ls=l.strip().split('\t')
		if ls[-1]=='KEEP':
			if mut_dict.has_key(ls[0]+'_'+ls[1]+'_'+ls[3]+'_'+ls[4]):
				mut_dict[ls[0]+'_'+ls[1]+'_'+ls[3]+'_'+ls[4]].add('MuTect')
			else:
				mut_dict[ls[0]+'_'+ls[1]+'_'+ls[3]+'_'+ls[4]]=set(['MuTect'])

	fout=open(outPath+sampleID+'/'+sampleID+'.consensus.snv.vcf','w')
	fout.write('##fileformat=VCFv4.0\n##source=Formated\n##phasing=none\n')
	fout.write(('{}\t'*9+'{}\n').format('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sampleID))
	for key in mut_dict.keys():
		if len(mut_dict[key])>=2:
			ks=key.split('_')
			fout.write(('{}\t'*9+'{}\n').format(ks[0],ks[1],'.',ks[2],ks[3],len(mut_dict[key]),'.','.',';'.join(mut_dict[key]),'.'))
	fout.close()
else:
	sys.exit('Unrecognized method of SNV calling. Exit!')

logging.debug('[DNA-seq] # DNA-seq based SNV calling completed.\n')
print '[DNA-seq] # DNA-seq based SNV calling completed.'



