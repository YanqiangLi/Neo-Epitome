import sys, os, csv, argparse, logging, datetime


def seq2HLA(readsFilesCaseRNA,runname,bindir):
	if os.path.exists(runname+'/hla_types-ClassI.HLAgenotype4digits')==False:
		readsFiles_split=readsFilesCaseRNA.split(',')
		os.system('mkdir -p '+runname)
		cmd = 'python '+bindir+'seq2HLA.py -1 '+readsFiles_split[0]+' -2 '+readsFiles_split[1]+' -r '+runname+'/hla_types'
		logging.debug('[seq2hla] Running command 1:'+cmd)
		os.system (cmd)
		if os.path.exists(runname+'/hla_types-ClassI.HLAgenotype4digits')==False:
			sys.exit('[seq2hla] # An Error Has Occured. seq2hla Incomplete. Exit!')
	else:
		logging.debug('[seq2hla] # Skipped seq2HLA.')
		print '[seq2hla] # Skipped seq2HLA.'

	HLA_type=[]
	for n,l in enumerate(open(runname+'/hla_types-ClassI.HLAgenotype4digits')):
		if n==0:
			continue
		ls=l.strip().split('\t')
		#print ls
		if ls[2]!='NA':
			if float(ls[2])<=0.05:
				HLA_type.append('HLA-'+ls[1].strip('\''))
		if ls[4]!='NA':
			if float(ls[4])<=0.05:
				HLA_type.append('HLA-'+ls[3].strip('\''))
		continue
	# for n,l in enumerate(open(runname+'-ClassII.HLAgenotype4digits')):
	# 	if n==0:
	# 		continue
	# 	ls=l.strip().split('\t')
	# 	if ls[2]!='NA':
	# 		if float(l[2])<=0.05:
	# 			HLA_type.append('HLA-'+l[1].strip('\''))
	# 	if ls[4]!='NA':
	# 		if float(ls[4])<=0.05:
	# 			HLA_type.append('HLA-'+l[3].strip('\''))
	# 	continue

	if len(HLA_type)==0:
		sys.exit('# [seq2hla] No HLA type predicted. Exit.')

	HLA_type_str=','.join(list(set(HLA_type)))
	return HLA_type_str	

def main():
	parser = argparse.ArgumentParser(description='NeoEpitome-seq2HLA (v1.0)')
	parser.add_argument('-p','--sampleID', default='NeoEpitomeOut', help='Sample ID will be used for output folder name and reads group name.')
	parser.add_argument('readsFilesCaseRNA',help='Tumor sample paired-end fastq files seperated by ",". ')

	args = parser.parse_args()
	os.system('mkdir -p '+args.sampleID)
	sampleID = args.sampleID.rstrip('/')
	runname = sampleID+'/hla_types'
	bindir = '../../Neo-Epitome/bin/seq2hla/'
	logging.basicConfig(level=logging.DEBUG,
	                    format='%(asctime)s %(message)s',
	                    filename=sampleID+'/NeoEpitome-seq2HLA.log'+ str(datetime.datetime.now())+'.txt' ,
	                    filemode='w')

	logging.debug('[seq2hla] # Start HLA typing.')
	print '[seq2hla] # Start HLA typing.'

	hla=seq2HLA(args.readsFilesCaseRNA, runname, bindir)

	print '[seq2hla] # Completed. HLA types: '+hla
	logging.debug('[seq2hla] # Completed. HLA types: '+hla)

if __name__ == '__main__':
	main()
