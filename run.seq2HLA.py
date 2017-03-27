import sys, os, csv
#python run.seq2HLA.py fastq1 fastq2
readfile1=sys.argv[1]
readfile2=sys.argv[2]
#runname='pred_hla'
runname=sys.argv[3].strip('/')+'/result'
os.system('mkdir -p '+runname)
cmd = 'python bin/seq2hla/seq2HLA.py -1 '+readfile1+' -2 '+readfile2+' -r '+runname
print cmd
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
print HLA_type_str	
