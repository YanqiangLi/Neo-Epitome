import sys, os, argparse, subprocess, logging, datetime
import multiprocessing as mp

def create_position_list(fin, outdir):
    fout=open(outdir+'/somatic.positions.txt','w')
    for line in open(fin):
        if not line.startswith("#"):
            ls=line.strip().split("\t")
            pos=(ls[0],ls[1],str(int(ls[1])))
            fout.write('\t'.join(pos)+'\n')
    logging.debug('[Bam2Cov] Created somatic mutation list.')
    print '[Bam2Cov] Created somatic mutation list.'
def calculating_cov(reference, outdir, tumor_bam, data_type):
    if os.path.exists(tumor_bam.rsplit('.',1)[0]+'.bai')==False and os.path.exists(tumor_bam+'.bai')==False:
        logging.debug('[Bam2Cov] Making index for '+tumor_bam+'\n')
        print '[Bam2Cov] Making index for '+tumor_bam
        os.system('samtools index '+tumor_bam)
    cmd1='bam-readcount -f '+str(reference)+' -w 1 -l '+outdir+'/somatic.positions.txt '+str(tumor_bam)+'> '+outdir+'/'+data_type+'.readcount'
    os.system(cmd1)
    logging.debug('[Bam2Cov] Running command: '+cmd1+'\n')
    print '[Bam2Cov] Running command: '+cmd1

def main():
    parser = argparse.ArgumentParser(description='NeoEpitome-BAM2Cov (v1.0)')
    parser.add_argument('vcf_input', help='input somatic mutation VCF file.')
    parser.add_argument('-p','--sampleID', default='NeoEpitomeOut', help='Sample ID will be used for output folder name and reads group name.')
    parser.add_argument('-r','--reference',help='Reference genome used for short reads mapping and GATK realigments.')
    parser.add_argument('--tumor-dna-bam', help='BAM file for tumor DNA seq.')
    parser.add_argument('--tumor-rna-bam', help='BAM file for tumor RNA seq.')
    args = parser.parse_args()

    fin=args.vcf_input
    sampleID=args.sampleID.rstrip('/')#outdir

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(message)s',
                        filename=sampleID+'/NeoEpitome-Bam2Cov.log'+ str(datetime.datetime.now())+'.txt' ,
                        filemode='w')

    create_position_list(fin,sampleID)

    mappings=[]
    mappings.append(mp.Process(target=calculating_cov,args=(args.reference, sampleID, args.tumor_dna_bam, 'DNA')))
    mappings[0].start()
    mappings.append(mp.Process(target=calculating_cov,args=(args.reference, sampleID, args.tumor_rna_bam, 'RNA')))
    mappings[1].start()
    mappings[0].join()
    mappings[1].join()

    logging.debug('[Bam2Cov] # Completed.')
    print '[Bam2Cov] # Completed.'

if __name__ == '__main__':
    main()
