import sys, os, csv, argparse, logging, datetime


def seq2HLA(readsFilesCaseRNA,runname,bindir):
	readsFiles_split=readsFilesCaseRNA.split(',')
	os.system('mkdir -p '+runname)
	cmd = 'python '+bindir+'seq2HLA.py -1 '+readsFiles_split[0]+' -2 '+readsFiles_split[1]+' -r '+runname
	logging.debug('[seq2hla] Running command 1:'+cmd)
	os.system (cmd)

	HLA_type=[]
	for n,l in enumerate(open(runname+'-ClassI.HLAgenotype4digits')):
		if n==0:
			continue
		ls=l.strip().split('\t')
		#print ls
		if ls[2]!='NA':
			if float(l[2])<=0.05:
				HLA_type.append('HLA-'+l[1].strip('\''))
		if ls[4]!='NA':
			if float(ls[4])<=0.05:
				HLA_type.append('HLA-'+l[3].strip('\''))
		continue
	for n,l in enumerate(open(runname+'-ClassII.HLAgenotype4digits')):
		if n==0:
			continue
		ls=l.strip().split('\t')
		if ls[2]!='NA':
			if float(l[2])<=0.05:
				HLA_type.append('HLA-'+l[1].strip('\''))
		if ls[4]!='NA':
			if float(ls[4])<=0.05:
				HLA_type.append('HLA-'+l[3].strip('\''))
		continue
	if len(HLA_type)==0:
		print '# No HLA type predicted. Exit.'
		exit()
	HLA_type_str=','.join(list(set(HLA_type)))
	return HLA_type_str	

def main():
	parser = argparse.ArgumentParser(description='NeoEpitome-seq2HLA (v1.0)')
	parser.add_argument('-p','--sampleID', default='NeoEpitomeOut', help='Sample ID will be used for output folder name and reads group name.')
	parser.add_argument('readsFilesCaseRNA',help='Tumor sample paired-end fastq files seperated by ",". ')

	args = parser.parse_args()
	os.system('mkdir -p '+args.sampleID)
	sampleID = args.sampleID.strip('/')
	runname = sampleID+'/hla_types'
	bindir = '/u/home/p/panyang/NeoEpitope/bin/seq2hla/'
	logging.basicConfig(level=logging.DEBUG,
	                    format='%(asctime)s %(message)s',
	                    filename=sampleID+'/NeoEpitome-seq2HLA.log'+ str(datetime.datetime.now())+'.txt' ,
	                    filemode='w')

	logging.debug('[seq2hla] # Start seq2HLA.')
	print '[seq2hla] # Start seq2HLA.'

	hla=seq2HLA(args.readsFilesCaseRNA, runname, bindir)

	print '[seq2hla] # Completed. HLA types: '+hla
	logging.debug('[seq2hla] # Completed. HLA types: '+hla)

if __name__ == '__main__':
	main()