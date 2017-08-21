import sys, argparse, subprocess, csv, re, numpy, datetime, os, uuid, logging
import multiprocessing as mp
from Bio.Blast import NCBIXML

now = datetime.datetime.now()
ID = str(uuid.uuid4()).split('-')[0]


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

def getProtSeq(id, seq_list):
	print id
	seq_left,seq_right=seq_list
	prot_seq_left = translate_dna(seq_left)
	prot_seq_left_blast = {}
	query_name = id.strip('>').split('_')
	if len(prot_seq_left)==0:
		return []
	else:
		fout=open('prot.tmp.fa','w')
		for p_seq in prot_seq_left:
			fout.write('>'+str(p_seq[1])+'_'+str(len(p_seq[0]))+'\n')
			fout.write(p_seq[0]+'\n')
		fout.close()
		cmd1='blastp -query prot.tmp.fa -db ~/nobackup-yxing/references.annotations/Homo_sapiens.GRCh38.pep.all.fa -matrix PAM30 -out prot.tmp.blast.xml -outfmt 5'
		logging.debug('[blastp] blast for '+';'.join([ps[0] for ps in prot_seq_left]))
		print '[blastp] blast for '+';'.join([ps[0] for ps in prot_seq_left])
		os.system(cmd1)
		for blast_result in NCBIXML.parse(open('prot.tmp.blast.xml')): 
			for desc in  blast_result.descriptions:
				hit_name = desc.title.split()
				hit_name_pos = hit_name[3].split(':')
				if hit_name_pos[2] == query_name[0].strip('chr'):
					if abs(int(query_name[1])-int(hit_name_pos[3]))<=10000000:
						print query_name
						print hit_name
						query_info=str(blast_result.query).split('_')
						prot_seq_left_blast[query_info[0]]=query_info[1]
						break
		prot_seqs=[]
		for seq,orf in prot_seq_left:
			if str(orf) in prot_seq_left_blast:
				prot_full_junc=translate_dna(seq_left[orf:]+seq_right,1, False)
				header_line=id+'_orf:'+str(orf)+'_junction:'+prot_seq_left_blast[str(orf)]+'\n'
				prot_seqs.append(header_line+prot_full_junc[0][0])
	return prot_seqs

def getProtSeq_fast(id, seq_list):
	seq_left,seq_right=seq_list
	prot_seq_left = translate_dna(seq_left)
	if len(prot_seq_left)==0:
		return []
	else:
		prot_seqs=[]
		for seq,orf in prot_seq_left:
			prot_full_junc=translate_dna(seq_left[orf:]+seq_right,1, False)
			header_line=id+'_orf:'+str(orf)+'_junction:'+str(len(prot_seq_left[0][0]))+'\n'
			prot_seqs.append(header_line+prot_full_junc[0][0])
	return prot_seqs

def translate_dna(seq, orf_option=3, rm_early_stop=True):
	codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

	prot_seqs = []
	for start_site in range(orf_option):
		cds = seq[start_site:]
		prot_seq = ''
		for n in range(0, len(cds), 3):
			current_triple = cds[n:n+3]

			if len(current_triple)!=3:
				break
			if current_triple.upper() in codontable:
				letter = codontable[current_triple.upper()]
				if letter=='-':
					prot_seq += letter
					break
				if current_triple == current_triple.upper():
					prot_seq += letter
				else:
					prot_seq += letter.lower()
			else:
				prot_seq += '*'
				print 'unknown codon',current_triple.upper() 
				break
		if prot_seq.find('_')!=-1:
			if rm_early_stop:
				continue #totally remove stop codon containing seq
			else:
				prot_seq = prot_seq.split('_')[0] #trim aa after stop codon
		prot_seqs.append((prot_seq,start_site))

	return prot_seqs

def getDNASeq(start_j, seq):
	DNA_seq_left = seq[:start_j+1]
	DNA_seq_right = seq[start_j+1:]
	return (DNA_seq_left,DNA_seq_right)

def rev_complement(seq):
	rev_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N', '|':'|'}
	rev_seq = seq[::-1]
	rev_comp_seq = map(lambda x: rev_dict[x], rev_seq)
	return ''.join(rev_comp_seq)

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
	logging.debug('[iedb] '+cmd)
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
	return parsed_dict


def main():
	#Define parameters
	parser = argparse.ArgumentParser(description='NeoEpitome-Seq2Antigen (v1.0)')
	parser.add_argument('fasta_input', help='input alternative splicing or fusion sequences in Fasta format.')
	parser.add_argument('-e', '--epitope-len-list', default='9,10,11', help='epitope length for prediction. Default is 9,10,11.')
	parser.add_argument('-a', '--hla-allele-list', default='HLA-A*01:01,HLA-B*08:01,HLA-C*07:01', help='a list of HLA types. Default is HLA-A*01:01,HLA-B*08:01,HLA-C*07:01.')
	parser.add_argument('-o', '--outdir', default= 'Result.'+ID, help='The output directory.')
	parser.add_argument('--iedb-local', default='~/nobackup-yxing/local/bin/mhc_i/src/', help='Specify local IEDB location if it is installed.')
	parser.add_argument('--ic50-cut-off', default=500, help='Cut-off based on median value of concensus predicted IC50 values. Default is 1000.')
	parser.add_argument('--protein-ms', type=argparse.FileType('r'), help='mzML format is recommended. Currently only support library search. Quantatitive pipeline is under development.')
	parser.add_argument('--protein-reference', default= 'data/.fasta', help='fasta format reference protein sequences. Details see MSGFPlus manu.')
	parser.add_argument('--protein-mod', default= 'data/mod.txt', help='protein modification file needed for labeled MS data library search. Details see MSGFPlus manu.')
	args = parser.parse_args()


	fin=args.fasta_input
	outdir=args.outdir.rstrip('/')
	os.system('mkdir -p '+outdir)
	iedb_path=args.iedb_local
	if args.iedb_local!=False:
		iedb_path=args.iedb_local.rstrip('/')
	hla_allele_list=loadSampleMHC(args.hla_allele_list)
	if hla_allele_list==[]:
		sys.exit("# No HLA Alleles. Exit.")

	logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=outdir+'/NeoEpitome-GenomeSeq2Antigen.log'+ str(datetime.datetime.now())+'.txt' ,
                    filemode='w')
	#Translation
	header=''
	if os.path.exists(outdir)==False:
		os.system('mkdir '+outdir)
	prot_fout=open(outdir+'/filtered_peptide.'+fin.split('/')[-1],'w')
	for row_num,l in enumerate(open(fin)):
		if row_num%2==0:
			header=l.strip()
		else:
			ls=l.strip()
			logging.debug('[# Query: '+str(row_num/2+1) +' Info: '+header)
			print '# Query:',str(row_num/2+1), 'Info:',header
			header=header.strip().split('_')
			DNA_seq_left,DNA_seq_right=getDNASeq(int(header[5].split(':')[1]), ls) 
			protseqs_plus=getProtSeq_fast('_'.join(header)+'_translation:forward',(DNA_seq_left,DNA_seq_right))

			DNA_seq_left,DNA_seq_right=getDNASeq(len(ls)-int(header[6].split(':')[1])-1, rev_complement(ls)) 
			protseqs_minus=getProtSeq_fast('_'.join(header)+'_translation:reverse',(DNA_seq_left,DNA_seq_right))

			logging.debug('[writting] Confirmed peptide sequence(s): '+str(len(protseqs_plus+protseqs_minus)))
			print '[writting] Confirmed peptide sequence(s): ', len(protseqs_plus+protseqs_minus)
			peptide_final=protseqs_plus+protseqs_minus
			for entry in peptide_final:
				prot_fout.write(entry+'\n')
	prot_fout.close()

	# Prediction
	prot_seq_list=[]
	prot_seq_header_list=[]
	mhc_type_dict = loadMHCType()

	pred_fout=open(outdir+'/pred.filtered_peptide.'+fin.split('/')[-1]+'.txt','w')
	pred_fout_seq=open(outdir+'/pred.filtered_peptide.'+fin.split('/')[-1]+'.fasta','w')
	for n,r in enumerate(open(outdir+'/filtered_peptide.'+fin.split('/')[-1])):
		if n%2==0:
			header=r.strip()
		else:
			rs=r.strip()
			hs=header.split('_')
			testable_pep=rs[max(0,int(hs[-1].split(':')[1])-11):min(int(hs[-1].split(':')[1])+11,len(rs))]
			prot_seq_list.append(testable_pep)
			prot_seq_header_list.append(header)
	if prot_seq_list==[''] or prot_seq_list==[]:
		logging.debug('No valid peptide.') 
		sys.exit('No valid peptide.')

	pred_result_mut = seqPredLocal(prot_seq_list,hla_allele_list,args.epitope_len_list.split(','), mhc_type_dict, iedb_path, outdir)
	for i,pred_allele_result_mut in enumerate(pred_result_mut):
		for k in pred_allele_result_mut:
			ks=k.split('_')
			seq_index=int(ks[1])-1
			allele=ks[0]
			epitope_len=ks[3]
			fusion_pep = prot_seq_list[seq_index]
			fusion_info = prot_seq_header_list[seq_index].split('_')
			pred_fout.write(('{}\t'*9+'{}\n').format('_'.join(fusion_info[0:3]),fusion_info[-3], fusion_info[-2],fusion_info[-1], int(ks[2]) - int(fusion_info[-1].split(':')[1]),epitope_len, allele, fusion_pep, pred_allele_result_mut[k][0],pred_allele_result_mut[k][1]))
			if pred_allele_result_mut[k][0]>int(args.ic50_cut_off):
				continue
			pred_fout_seq.write('>{}\n{}\n'.format(prot_seq_header_list[seq_index]+'_'+epitope_len+'_'+allele+'_'+ks[2],fusion_pep))#FASTA with full length mutant protein sequences, let the library search algorithm to remove duplicated fragments.

if __name__ == '__main__':
	main()

