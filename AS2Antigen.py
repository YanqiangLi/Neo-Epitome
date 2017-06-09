import sys, os
from pyfaidx import Fasta

def trim_prot(prot_seq):
	pass

def prot_seq(id, seq_list, proteome_seq, id_mapper):
	seq_left,seq_right=seq_list
	prot_seq_left = translate_dna(seq_left)
	#print prot_seq_left
	filtered_orf=[]
	for (ps,orf) in prot_seq_left:
		for prot_id in id_mapper[id]:
			#print prot_id, str(proteome_seq[prot_id])
			if str(proteome_seq[prot_id]).find(ps.upper())!=-1:
				filtered_orf.append((ps,orf))
				break
	prot_seqs=[]
	for fs in filtered_orf:
		prot_full_junc=translate_dna(seq_left[fs[1]:]+seq_right,1)
		prot_seqs.append(prot_full_junc[0][0])
	return prot_seqs

def translate_dna(seq, orf_option=3):
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
			if current_triple.upper() in codontable:
				letter = codontable[current_triple.upper()]
				if current_triple == current_triple.upper():
					prot_seq += letter
				else:
					prot_seq += letter.lower()
			else:
				#prot_seq += '-'
				continue
		prot_seqs.append((prot_seq,start_site))

	return prot_seqs

def rev_complement(seq):
	rev_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N', '|':'|'}
	rev_seq = seq[::-1]
	rev_comp_seq = map(lambda x: rev_dict[x], rev_seq)
	return ''.join(rev_comp_seq)

def read_genome(genome_file):
	print 'Loading genome file...'
	genome = Fasta(genome_file, as_raw=True)
	print 'Loading genome file done.'
	return genome

def get_seqs(flank, sj_event, genomic_seq):##need to change
	sjs=sj_event.split('_')
	left_bound=max(int(sjs[2])-flank,int(sjs[1]))
	right_bound=min(int(sjs[3])+flank,int(sjs[4]))
	seq_left=str(genomic_seq[sjs[0]][left_bound:int(sjs[2])])
	seq_right=str(genomic_seq[sjs[0]][int(sjs[3]):right_bound])
	if sjs[5]=='-':
		seq_right=rev_complement(seq_left)
		seq_left=rev_complement(seq_right)
	return [seq_left,seq_right]

def read_proteome(proteome_file):
	print 'Loading proteome file...'
	proteome = Fasta(proteome_file, as_raw=True)
	print 'Loading proteome file done.'
	return proteome

def gene2prot_dict():
	id_mapper={}
	for l in open('P2G.mapper.txt'):
		ls=l.strip().split('\t')
		if ls[1] in id_mapper:
			id_mapper[ls[1]].append(ls[0])
		else:
			id_mapper[ls[1]]=[ls[0]]
	return id_mapper

genome_file = sys.argv[2]
proteome_file = sys.argv[3]
genomic_seq = read_genome(genome_file)
proteome_seq = read_proteome(proteome_file)
id_mapper = gene2prot_dict()
flank=36
fin=sys.argv[1]
n=0
for l in open(fin):
	ls=l.strip().split()
	sj_left=ls[3]+'_'+str(int(ls[7]))+'_'+str(int(ls[8]))+'_'+str(int(ls[5]))+'_'+str(int(ls[6]))+'_'+ls[4]+'_'+ls[1].strip('"').split('_')[0]
	sj_right=ls[3]+'_'+str(int(ls[5]))+'_'+str(int(ls[6]))+'_'+str(int(ls[9]))+'_'+str(int(ls[10]))+'_'+ls[4]+'_'+ls[1].strip('"').split('_')[0]

## just need to set the flanking region, and use genomic_seq for DNA seq (note the strandness)
	DNA_seqs_left = get_seqs(flank, sj_left, genomic_seq)
	DNA_seqs_right = get_seqs(flank, sj_right, genomic_seq)
	prot_seqs_left = prot_seq(sj_left.split('_')[6], DNA_seqs_left, proteome_seq, id_mapper)
	prot_seqs_right = prot_seq(sj_right.split('_')[6], DNA_seqs_right, proteome_seq, id_mapper)
	print 
	print prot_seqs_left
	print prot_seqs_right
	n+=1
	if n==100:
		break
	# prot_seqs_right = translate_dna(DNA_seqs_right)
	# prot_seqs_left = translate_dna(DNA_seqs_left)


	# tr_prot_seq_left = trim_prot(prot_seqs_left, trim_RK)
	# tr_prot_seq_right = trim_prot(prot_seqs_right, trim_RK)
	
	# left_trim = prot_seq.find(tr_prot_seq_left)
	
	# if len(tr_prot_seq) < 5 or tr_prot_seq.upper() == tr_prot_seq:
	# 	np = np + 1
	# 	continue
	# out.write(id_line + '_startTrim:'+str(left_trim)+'_ORF:' + str(j) + '\n' + tr_prot_seq + '\n')
