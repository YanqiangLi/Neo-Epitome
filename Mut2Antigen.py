import sys, argparse, subprocess, csv, re, numpy, datetime, os, uuid
import multiprocessing as mp

#python Mut2Antigen.py --cov-site-dna example_data/snvs.bam_readcount --exp-gene-rna example_data/genes.fpkm_tracking example_data/input.vcf
#python Mut2Antigen.py --step 50 -a HLA-A*01:01,HLA-A*02:01,HLA-B*57:01,HLA-B*15:01,HLA-C*06:02,HLA-C*03:04 --exp-gene-rna neoantigen_2899/genes.fpkm_tracking --exp-isoform-rna neoantigen_2899/isoforms.fpkm_tracking --cov-site-dna neoantigen_2899/DNA.readcount --cov-site-rna neoantigen_2899/RNA.readcount neoantigen_2899/neoantigen_2899.snp.vep.vcf
# USER need to set VEP path and VEP cache path for annotation step to be functional.

now = datetime.datetime.now()
ID = str(uuid.uuid4()).split('-')[0]

#VEP settings
PLUGIN_PATH='~/nobackup-yxing/.vep/Plugins/'
CACHE_PATH='~/nobackup-yxing/.vep/'
VEP_PATH='/u/home/p/panyang/nobackup-yxing/ensembl-vep-release-88'


def getTerminalSize():
	env = os.environ
	def ioctl_GWINSZ(fd):
		try:
			import fcntl, termios, struct, os
			cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
		'1234'))
		except:
			return
		return cr
	cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
	if not cr:
		try:
			fd = os.open(os.ctermid(), os.O_RDONLY)
			cr = ioctl_GWINSZ(fd)
			os.close(fd)
		except:
			pass
	if not cr:
		cr = (env.get('LINES', 25), env.get('COLUMNS', 80))
	return int(cr[1]), int(cr[0])
		
def progressBar(progress,out = sys.stdout):
	status = ""
	if isinstance(progress, int):
		progress = float(progress)
	if not isinstance(progress, float):
		progress = 0
		status = "error: progress var must be float"
	if progress < 0:
		progress = 0
		status = "Halt..."
	if progress >= 1:
		progress = 1
		status = "Done...\r\n"
	barLength = getTerminalSize()[0] - 25 - len(status)
	if barLength<10: 
		barLength = 10
	block = int(round(barLength*progress))
	text = "\rPercent: [{0}] {1:.2f}% {2}".format( "="*block+">" + "-"*(barLength-block-1), progress*100, status)
	out.write(text)
	out.flush()

def fileLen(fin):
	with open(fin) as f:
		for i, l in enumerate(f):
			pass
	return i

def loadSampleMHC(f_s_in):
	if f_s_in.find(',')!=-1:
		hla_allele_list=f_s_in.split(',')
	else:
		n=0
		hla_allele_list=[]
		for l in open(f_s_in):
			if n==0:
				n+=1
				continue
			ls=l.strip().split('\t')
			hla_allele_list.append('HLA-'+ls[1].strip("'"))
			hla_allele_list.append('HLA-'+ls[3].strip("'"))
	return hla_allele_list

def vcfAnnotation(vcf_input, outdir, assembly):
	cmd1=VEP_PATH+'/vep --species homo_sapiens --assembly '+assembly+' --offline --dir_cache '+CACHE_PATH+' --dir_plugins '+PLUGIN_PATH+' --input_file '+vcf_input+' --vcf --symbol --terms SO --plugin Downstream --plugin Wildtype --uniprot --canonical -hgvs --output_file '+outdir+'/'+vcf_input.split('/')[-1].split('.vcf')[0]+'.vep.vcf'
	os.system(cmd1)
	if os.path.exists(outdir+'/'+vcf_input.split('/')[-1].split('.vcf')[0]+'.vep.vcf')==False:
		sys.exit('[VEP annotation] # An Error Has Occured. VEP Annotation Incomplete. Exit!')
	return outdir+'/'+vcf_input.split('/')[-1].split('.vcf')[0]+'.vep.vcf'

def loadMHCType():
	mhc_type_dict={x.strip():'MHCI' for x in open('data/MHCI.txt')}
	mhc_type_dict.update({x.strip():'MHCII' for x in open('data/MHCII.txt')})
	return mhc_type_dict

def mhcPredType(mhc_type_dict, hla_allele):
	if hla_allele in mhc_type_dict:
		if mhc_type_dict[hla_allele]=='MHCI':
			return 'mhci'
		else:
			return 'mhcii'
	else:
		sys.exit("# Unsupported HLA type: "+hla_allele+". Exit! ")

def parsePred(stdout):
	ref_dict={}
	ref_out=csv.DictReader(stdout,delimiter='\t')
	for r in ref_out:
		ic50=[k for k in r.keys() if k.find('ic50')!=-1]
		p = re.compile('\d+(\.\d+)?')
		ic50_value=[float(r[k]) for k in ic50 if p.match(r[k]) != None]
		ref_dict[r['allele']+'_'+r['seq_num']+'_'+r['start']+'_'+r['length']]=[numpy.median(ic50_value),r['peptide']]
	return ref_dict
		
def seqPred(prot_seq_list,hla_allele_list,epitope_len_list, mhc_type_dict):
	prot_seq= '%0A'.join('%3Epredict'+str(i)+'%0A'+seq for i,seq in enumerate(prot_seq_list))
	stdout_total=[]
	for hla_allele in hla_allele_list:
		cmd = 'curl --data method=recommended&sequence_text='+prot_seq+'&allele='+','.join([hla_allele]*len(epitope_len_list))+'&length='+','.join(epitope_len_list)+' http://tools-cluster-interface.iedb.org/tools_api/'+mhcPredType(mhc_type_dict, hla_allele)+'/'
		args = cmd.split()
		process = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = process.communicate()
		stdout_total+=[stdout]
	return [parsePred(item.rstrip().splitlines()) for item in stdout_total]

def localIEDBCommand(iedb_path, hla_allele, epitope_len, outdir):
	cmd = 'python '+iedb_path+'/predict_binding.py IEDB_recommended '+hla_allele+' '+epitope_len+' '+outdir+'/tmp.fasta > '+outdir+'/tmp_'+hla_allele.replace('*','_').replace(':','_')+'_'+epitope_len+'_iedb.txt'
	os.system(cmd)

def seqPredLocal(prot_seq_list,hla_allele_list,epitope_len_list, mhc_type_dict, iedb_path, outdir):
	with open(outdir+'/tmp.fasta', 'w') as fw:
		for seq in prot_seq_list:
			if seq.find('X')!=-1 or len(seq)<8:
				continue
			fw.writelines('>seq\n'+seq+'\n')
	predicting=[]
	n=0
	file_list=[]
	for hla_allele in hla_allele_list:
		for epitope_len in epitope_len_list:
			os.system('rm -f '+outdir+'/tmp_'+hla_allele.replace('*','_').replace(':','_')+'_'+epitope_len+'_iedb.txt')
			file_list.append(outdir+'/tmp_'+hla_allele.replace('*','_').replace(':','_')+'_'+epitope_len+'_iedb.txt')
			predicting.append(mp.Process(target=localIEDBCommand,args=(iedb_path, hla_allele, epitope_len, outdir)))
			predicting[n].start()
			n+=1
	for i in range(0,n):
		predicting[i].join()
	parsed_dict=[parsePred(open(tmp_file)) for tmp_file in file_list]
	os.system('rm -f '+outdir+'/tmp*')
	return parsed_dict

def mergePeps2Database(fastafile, reference, outdir): 
	##need to change
	os.system("cat "+reference +' '+ fastafile + ">" + outdir + "/merge_" + fastafile.split('/')[-1])
	return os.path.basename(outdir + "/merge_" + fastafile.split('/')[-1])

def findCSQNum(info):
	csq_num=-1
	for t,info_field in enumerate(info):
		if info_field.startswith('CSQ='):
			csq_num=t
	if csq_num!=-1:
		return csq_num
	else:
		sys.exit('# VCF file annotation error. Exit!')

def frameshiftPred(ls, annotation, iedb_path, mhc_type_dict, hla_allele_list, epitope_len_list, ic50_cut_off, fout, fout_seq):
	if ls[0].startswith('chr'):
		dna_pos_fs=ls[0]+':'+ls[1]
	else:
		dna_pos_fs='chr'+ls[0]+':'+ls[1]
	dna_var_fs=ls[3]+'/'+ls[4]
	gene_name_fs=annotation['SYMBOL']
	gene_ac_fs=annotation['Gene']
	trnscrpt_ac_fs=annotation['Feature']
	prot_ac_fs=annotation['HGVSp'].split(':')[0]
	prot_var_fs=annotation['HGVSp'].split(':')[1].split('.')[1]
	prot_pos_fs=annotation['Protein_position'].split('-')[0]
	prot_seq_ref_fs=annotation['WildtypeProtein'].strip()
	prot_seq_mut_fs=annotation['DownstreamProtein'].strip()
	if iedb_path:
		pred_result_mut=seqPredLocal([prot_seq_mut_fs], hla_allele_list, epitope_len_list, mhc_type_dict, iedb_path, fout.rsplit('/',1)[0])
	else:
		pred_result_mut=seqPred([prot_seq_mut_fs], hla_allele_list, epitope_len_list, mhc_type_dict)
	for i,pred_allele_result_mut in enumerate(pred_result_mut):
		for k in pred_allele_result_mut:
			foc='-'
			ks=k.split('_')
			start_pos_seq=int(ks[2])
			allele=ks[0]
			epitope_len=ks[3]
			fout.write(('{}\t'*15+'{}\n').format(dna_pos_fs,dna_var_fs,gene_name_fs,gene_ac_fs,trnscrpt_ac_fs,prot_ac_fs,prot_pos_fs,prot_var_fs,start_pos_seq,epitope_len,allele,'-',pred_allele_result_mut[k][0],foc,pred_allele_result_mut[k][1].lower(),pred_allele_result_mut[k][1]))
			if pred_allele_result_mut[k][0]>int(ic50_cut_off):
				continue
			fout_seq.write('>{}\n{}\n'.format(dna_pos_fs+'_'+dna_var_fs+'_'+gene_name_fs+'_'+prot_ac_fs+'_'+str(prot_pos_fs)+'_'+str(start_pos_seq)+'_'+prot_var_fs+'_'+str(epitope_len)+'_'+allele,prot_seq_mut_fs)) #FASTA with full length mutant protein sequences, let the library search algorithm to remove duplicated fragments.

def mutationPipeline(fin, hla_allele_list, epitope_len_list, step, ic50_cut_off, outdir, iedb_path, fs):
	# Defining varaibles 
	fin_len=fileLen(fin)
	peptide_len = max(map(int,epitope_len_list))
	fout=open(outdir+'/'+fin.split('/')[-1]+'.step1.out','w')
	fout.write(('{}\t'*15+'{}\n').format('GenomicLocation','GenomicMut','GeneName','GeneAC','TranscriptAC','ProteinAC','ProteinPos','ProteinVar','PeptideStartPos','EpitopeLength','HLAAllele','RefPredScore','MutPredScore','FoC','MutSeq','EpitopeSeq'))
	fout_seq=open(outdir+'/'+fin.split('/')[-1]+'.step1.out.fa','w')
	seq_num,dna_pos,dna_var,gene_name,gene_ac,trnscrpt_ac,prot_ac,prot_pos,prot_var,prot_seq_ref,prot_seq_mut,start_pos,mut_seq,ref_seq=([] for i in range(14))
	mhc_type_dict=loadMHCType()
	# Parse Input file
	csq_num=-1
	info_header=[]
	for file_index,l in enumerate(open(fin)): 
		progressBar(file_index/(0.0+fin_len))
		if l.startswith('#'):
			if l.startswith('##INFO=<ID=CSQ'):
				csq_header=l.strip().strip('>').split('Format: ')[1].strip('"').split('|')
			continue
		if len(info_header):
			sys.exit('# VCF file annotation error. Exit!')
		ls=l.strip().split('\t')
		info=ls[7].split(';')
		if csq_num==-1:
			csq_num=findCSQNum(info)
		for csq_value_str in info[csq_num].split(','):
			csq_value=csq_value_str.split('|')
			annotation=dict(zip(csq_header,csq_value))

			consequence= annotation['Consequence'].split('&')
			if 'missense_variant' in consequence or 'inframe_insertion' in consequence or 'inframe_deletion' in consequence:
				if ls[0].startswith('chr'):
					dna_pos+=[ls[0]+':'+ls[1]]
				else:
					dna_pos+=['chr'+ls[0]+':'+ls[1]]
				dna_var+=[ls[3]+'/'+ls[4]]
				gene_name+=[annotation['SYMBOL']]
				gene_ac+=[annotation['Gene']]
				trnscrpt_ac+=[annotation['Feature']]
				prot_ac+=[annotation['HGVSp'].split(':')[0]]
				prot_var+=[annotation['Amino_acids']]
				prot_vars=prot_var[-1].split('/')
				[ref_aa,mut_aa]=prot_vars
				if prot_vars[1]=='-': #handle pure insertion
					mut_aa=''
				if prot_vars[0]=='-': #handle pure deletion
					ref_aa=''
				prot_pos+=[annotation['Protein_position'].split('-')[0]]
				prot_seq_ref+=[annotation['WildtypeProtein'].strip()]
				prot_seq_mut+=[annotation['WildtypeProtein'].strip()[:int(prot_pos[-1])-1]+mut_aa.lower()+annotation['WildtypeProtein'].strip()[int(prot_pos[-1])-1+len(ref_aa):]]
				start_pos+=[max(0,int(prot_pos[-1])-1-(peptide_len-1))]
				end_pos=min(len(prot_seq_mut[-1]),int(prot_pos[-1])-1+(peptide_len-1)+len(mut_aa)-len(ref_aa)+1)
				mut_seq+=[prot_seq_mut[-1][start_pos[-1]:end_pos]]
				ref_seq+=[prot_seq_ref[-1][start_pos[-1]:end_pos]]
			elif 'frameshift_variant' in consequence and fs==True:
				frameshiftPred(ls, annotation, iedb_path, mhc_type_dict, hla_allele_list, epitope_len_list, ic50_cut_off, fout, fout_seq)
				continue

			if len(dna_pos)==step or file_index==fin_len-1:
				if len(dna_pos)==0:
					continue
				if iedb_path:
					pred_result_mut=seqPredLocal(mut_seq, hla_allele_list, epitope_len_list, mhc_type_dict, iedb_path, outdir)
					pred_result_ref=seqPredLocal(ref_seq, hla_allele_list, epitope_len_list, mhc_type_dict, iedb_path, outdir)
				else:
					pred_result_mut=seqPred(mut_seq, hla_allele_list, epitope_len_list, mhc_type_dict)
					pred_result_ref=seqPred(ref_seq, hla_allele_list, epitope_len_list, mhc_type_dict)
				for i,pred_allele_result_mut in enumerate(pred_result_mut):
					for k in pred_allele_result_mut:
						foc=pred_allele_result_mut[k][0]/pred_result_ref[i][k][0]
						ks=k.split('_')
						seq_index=int(ks[1])-1
						start_pos_seq=int(start_pos[seq_index])+int(ks[2])
						#filtering predictions beyond the epitope length
						allele=ks[0]
						epitope_len=ks[3]
						
						if start_pos_seq <=(prot_pos[seq_index]-epitope_len) or start_pos_seq>prot_pos[seq_index]:
							continu
						fout.write(('{}\t'*15+'{}\n').format(dna_pos[seq_index],dna_var[seq_index],gene_name[seq_index],gene_ac[seq_index],trnscrpt_ac[seq_index],prot_ac[seq_index],prot_pos[seq_index],prot_var[seq_index],start_pos_seq,epitope_len,allele,pred_result_ref[i][k][0],pred_allele_result_mut[k][0],foc,mut_seq[seq_index],pred_allele_result_mut[k][1]))
						if pred_allele_result_mut[k][0]>int(ic50_cut_off):
							continue
						fout_seq.write('>{}\n{}\n'.format(dna_pos[seq_index]+'_'+dna_var[seq_index]+'_'+gene_name[seq_index]+'_'+prot_ac[seq_index]+'_'+str(prot_pos[seq_index])+'_'+str(start_pos_seq)+'_'+prot_var[seq_index]+'_'+str(epitope_len)+'_'+allele+'_'+pred_allele_result_mut[k][1],prot_seq_mut[seq_index])) #FASTA with full length mutant protein sequences, let the library search algorithm to remove duplicated fragments.
				dna_pos,dna_var,prot_ac,gene_name,gene_ac,trnscrpt_ac,prot_pos,prot_var,prot_seq_ref,prot_seq_mut,start_pos,mut_seq,ref_seq=([] for i in range(13))

def appendExpandCov(args):
	msg=False
	if (args.exp_gene_rna is not None
		or args.exp_isoform_rna is not None
		or args.cov_site_dna is not None
		or args.cov_site_rna is not None):
		msg=True
	return msg

def parseCoverageFile(bam_readcount_file): 
	cov={}
	for l in bam_readcount_file:
		ls=l.strip().split('\t')
		(chrom, pos, ref, depth, base_info)= (ls[0], ls[1], ls[2], ls[3], ls[4:])
		if chrom.startswith('chr'):
			genome_loc=chrom+':'+pos
		else:
			genome_loc='chr'+chrom+':'+pos
		if genome_loc not in cov:
			cov[genome_loc]={}
		for base in base_info[1:]:
			base_metric=base.split(':')
			if int(depth)==0:
				cov[genome_loc][base_metric[0]]='\t'.join(['NA']*4)
				continue
			strand_dir=str(int(base_metric[5])/(0.000001+int(base_metric[6])))
			if int(base_metric[6])==0:
				strand_dir=str(float('inf'))
			mut_ratio=str(int(base_metric[1])/(0.0+int(depth)))
			cov[genome_loc][base_metric[0]]=depth+'\t'+mut_ratio+'\t'+base_metric[2]+'\t'+strand_dir
	return cov

def parseQuantFile(fin_Quant_gene_case):
	exp={}
	for l in csv.DictReader(fin_Quant_gene_case,dialect='excel-tab'):
		exp[l['tracking_id'].split('.')[0]] = l['FPKM']
	return exp

def appendSeqInfo(args):
	fin=args.vcf_input
	outdir=args.outdir.strip('/')
	cov_dna={}
	print '# Loading expression/coverage information.'
	if args.cov_site_dna:
		cov_dna=parseCoverageFile(args.cov_site_dna)
	cov_rna={}
	if args.cov_site_rna:
		cov_rna=parseCoverageFile(args.cov_site_rna)
	exp_gene={}
	if args.exp_gene_rna:
		exp_gene=parseQuantFile(args.exp_gene_rna)
	exp_isoform={}
	if args.exp_isoform_rna:
		exp_isoform=parseQuantFile(args.exp_isoform_rna)
	fout=open(outdir+'/'+fin.split('/')[-1]+'.step1.out.step2.out','w')
	print '# Matching information to predicted neo-epitopes.'
	for n,l in enumerate(open(outdir+'/'+fin.split('/')[-1]+'.step1.out')):
		if n==0:
			fout.write('{}\t{}\t{}\t{}\t{}\n'.format(l.strip(),'ReadCov_DNA\tMutReadRatio_DNA\tMutAvgMapQual_DNA\tMutStrDir_DNA','ReadCov_RNA\tMutReadRatio_RNA\tMutAvgMapQual_RNA\ttMutStrDir_RNA','geneFPKM','transcriptFPKM'))
			continue
		ls=l.strip().split('\t')
		try:
			cov_dna_str=cov_dna[ls[0]][ls[1].split('/')[1]]
		except:
			cov_dna_str='\t'.join(['NA']*4)
		try:
			cov_rna_str=cov_rna[ls[0]][ls[1].split('/')[1]]
		except:
			cov_rna_str='\t'.join(['NA']*4)

		exp_gene_str='NA'
		if ls[3].split('.')[0] in exp_gene:
			exp_gene_str=exp_gene[ls[3].split('.')[0]]
		
		exp_isoform_str='NA'
		if ls[4].split('.')[0] in exp_isoform:
			exp_isoform_str=exp_isoform[ls[4].split('.')[0]]
			
		fout.write('{}\t{}\t{}\t{}\t{}\n'.format(l.strip(),cov_dna_str,cov_rna_str,exp_gene_str,exp_isoform_str))

def main():
	#Define parameters
	parser = argparse.ArgumentParser(description='NeoEpitome-Mut2Antigen (v1.0)')
	parser.add_argument('vcf_input', help='input annotated somatic mutation VCF file.')
	parser.add_argument('-j', '--junction-input', help='input of somatic junctions file.')
	parser.add_argument('-e', '--epitope-len-list', default='9,10,11', help='epitope length for prediction. Default is 9,10.')
	parser.add_argument('-a', '--hla-allele-list', default='HLA-A*01:01,HLA-B*07:02', help='a list of HLA types. Default is HLA-A*01:01,HLA-B*01:01.')
	parser.add_argument('-o', '--outdir', default= 'Result.'+ID, help='The output directory.')
	parser.add_argument('--vcf-annotation', action='store_true', help='Specify local IEDB location if it is installed.')
	parser.add_argument('--assembly',default='GRCh37', help='Specify the annotation version: GRCh37 or GRCh38. Default is GRCh37.')
	parser.add_argument('--iedb-local', default=False, help='Specify local IEDB location if it is installed.')
	parser.add_argument('--step', default=100, help='Number of entries per time sending to prediction. Default is 50.')
	parser.add_argument('--ic50-cut-off', default=500, help='Cut-off based on median value of concensus predicted IC50 values. Default is 1000.')
	parser.add_argument('--frame-shift',action='store_true', help='Parse and predict frame-shift mutation outcomes.')
	parser.add_argument('--protein-ms', type=argparse.FileType('r'), help='mzML format is recommended. Currently only support library search. Quantatitive pipeline is under development.')
	parser.add_argument('--protein-reference', default= 'data/.fasta', help='fasta format reference protein sequences. Details see MSGFPlus manu.')
	parser.add_argument('--protein-mod', default= 'data/mod.txt', help='protein modification file needed for labeled MS data library search. Details see MSGFPlus manu.')
	parser.add_argument('--exp-gene-rna', type=argparse.FileType('r'), help='genes.fpkm_tracking file from Cufflinks')
	parser.add_argument('--exp-isoform-rna', type=argparse.FileType('r'), help='isoforms.fpkm_tracking file from Cufflinks')
	parser.add_argument('--cov-site-dna', type=argparse.FileType('r'), help='bam-readcount output file for tumor DNA BAM and snvs')
	parser.add_argument('--cov-site-rna', type=argparse.FileType('r'), help='bam-readcount output file for tumor RNA BAM and snvs')
	args = parser.parse_args()

	fin=args.vcf_input
	outdir=args.outdir.strip('/')
	os.system('mkdir -p '+outdir)
	if args.vcf_annotation:
		if fin.split('/')[-1].endswith('.vep.vcf')!=True:
			print str(datetime.datetime.now()),'# Run VEP annotation.'
			fin=vcfAnnotation(fin, outdir, args.assembly)
			print fin
		else:
			sys.exit('# The input VCF file seems already annotated. Please rename the file without ".vep" or remove argument "--vcf-annotation"')
	epitope_len_list=args.epitope_len_list.split(',')
	if min(epitope_len_list)<8:
		sys.exit("# The request epitope length is too small. Exit.")
	iedb_path=args.iedb_local
	if args.iedb_local!=False:
		iedb_path=args.iedb_local.rstrip('/')
	hla_allele_list=loadSampleMHC(args.hla_allele_list)
	if hla_allele_list==[]:
		sys.exit("# No HLA Alleles. Exit.")
	print str(datetime.datetime.now()),'# Searching Neoepitopes for allele types',','.join(hla_allele_list), 'with ', ','.join(epitope_len_list),'long...'
	
	# Step 1. Predicting the Neo-epitope by HLA binding affinity. Generating TSV and FASTA files. 
	if os.path.exists(outdir+'/'+fin.split('/')[-1]+'.step1.out')==False:
		mutationPipeline(fin, hla_allele_list, epitope_len_list, int(args.step), args.ic50_cut_off, outdir, iedb_path, args.frame_shift)
		print str(datetime.datetime.now()),'# Finished searching and predicting. The result is saved as: '+outdir+'/'+fin.split('/')[-1]+'.step1.out \n'
	else:
		print str(datetime.datetime.now()),'# Skipped searching and predicting. The file alraedy exists: '+outdir+'/'+fin.split('/')[-1]+'.step1.out \n'

	# Step 2 (Optional). Confirming with coverage and epxression level information. Appending to TSV.
	if appendExpandCov(args):
		print str(datetime.datetime.now()),'# Confirming with coverage and/or expression information from DNA/RNA sequencing.'
		appendSeqInfo(args)
		print str(datetime.datetime.now()),'# Finished appending information. The result is saved as :'+outdir+'/'+fin.split('/')[-1]+'.step1.out.step2.out \n'
  
	# Step 3 (Optional). Searching proteome level evidence by MS data library search. Appending to TSV.
	if args.protein_ms:
		print str(datetime.datetime.now()),'# Searching for protein evidence of predicted neo-epitopes.'
		custom_database = mergePeps2Database(outdir+'/'+fin.split('/')[-1]+'.step1.out.fa', args.protein_reference, outdir)
		#print str(now),'# Running MSGFPlus for proteome data searching.. The result is saved as: '+outdir+'/'+fin.split('/')[-1]+'.step1.out '
		#os.system('java -Xmx3500M -jar MSGFPlus.jar -s '+args.protein_ms_input+' -d '+custom_database+' -t 20ppm -protocol 2 -ti -1,2 -ntt 2 -tda 1 -inst 3 -m 1 -mod '+args.protein_mod)
		## os.system('java -Xmx3500M -jar MSGFPlus.jar -s test.mzXML -d IPI_human_3.79.fasta -t 0.5Da,2.5Da -ntt 2 -protocol 2 -tda 1 -o testMSGFPlus.mzid -mod Mods.txt')
		print str(datetime.datetime.now()),'# Finished searching and filtering.\n'
	print str(datetime.datetime.now()),'# Neo-Epitome Detection is Completed.\n'
if __name__ == '__main__':
	main()
