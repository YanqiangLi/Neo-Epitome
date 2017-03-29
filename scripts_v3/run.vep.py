import sys,os

input_vcf=sys.argv[1]
outdir=sys.argv[2]
pluginpath='~/nobackup-yxing/.vep/Plugins/'
veppath='~/nobackup-yxing/variant_effect_predictor'

cmd1='perl '+veppath+'/variant_effect_predictor.pl --input_file '+input_vcf+' --format vcf --output_file '+outdir+'/'+input_vcf.split('/')[-1].split('.vcf')[0]+'.vep.vcf --vcf --symbol --terms SO --plugin Downstream --plugin Wildtype --dir_plugins '+pluginpath+' --cache --uniprot --canonical -hgvs'
print 'run ',cmd1
os.system(cmd1)