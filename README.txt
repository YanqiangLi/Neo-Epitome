##Dependency:
1. Python 2.7 (biopython, numpy)
2. R
3. BWA,Samtools,
3. IEDB (optional)

##VCF mode:
usage: Mut2Antigen.py [-h] [-j JUNCTION_INPUT] [-e EPITOPE_LEN_LIST]
                      [-a HLA_ALLELE_LIST] [-o OUTDIR] [--vcf-annotation]
                      [--assembly ASSEMBLY] [--iedb-local IEDB_LOCAL]
                      [--step STEP] [--ic50-cut-off IC50_CUT_OFF]
                      [--frame-shift] [--protein-ms PROTEIN_MS]
                      [--protein-reference PROTEIN_REFERENCE]
                      [--protein-mod PROTEIN_MOD]
                      [--exp-gene-rna EXP_GENE_RNA]
                      [--exp-isoform-rna EXP_ISOFORM_RNA]
                      [--cov-site-dna COV_SITE_DNA]
                      [--cov-site-rna COV_SITE_RNA]
                      vcf_input

NeoEpitome-Mut2Antigen (v1.0)

positional arguments:
  vcf_input             input annotated somatic mutation VCF file.

optional arguments:
  -h, --help            show this help message and exit
  -j JUNCTION_INPUT, --junction-input JUNCTION_INPUT
                        input of somatic junctions file.
  -e EPITOPE_LEN_LIST, --epitope-len-list EPITOPE_LEN_LIST
                        epitope length for prediction. Default is 9,10.
  -a HLA_ALLELE_LIST, --hla-allele-list HLA_ALLELE_LIST
                        a list of HLA types. Default is
                        HLA-A*01:01,HLA-B*01:01.
  -o OUTDIR, --outdir OUTDIR
                        The output directory.
  --vcf-annotation      Specify local IEDB location if it is installed.
  --assembly ASSEMBLY   Specify the annotation version: GRCh37 or GRCh38.
                        Default is GRCh37.
  --iedb-local IEDB_LOCAL
                        Specify local IEDB location if it is installed.
  --step STEP           Number of entries per time sending to prediction.
                        Default is 50.
  --ic50-cut-off IC50_CUT_OFF
                        Cut-off based on median value of concensus predicted
                        IC50 values. Default is 1000.
  --frame-shift         Parse and predict frame-shift mutation outcomes.
  --protein-ms PROTEIN_MS
                        mzML format is recommended. Currently only support
                        library search. Quantatitive pipeline is under
                        development.
  --protein-reference PROTEIN_REFERENCE
                        fasta format reference protein sequences. Details see
                        MSGFPlus manu.
  --protein-mod PROTEIN_MOD
                        protein modification file needed for labeled MS data
                        library search. Details see MSGFPlus manu.
  --exp-gene-rna EXP_GENE_RNA
                        genes.fpkm_tracking file from Cufflinks
  --exp-isoform-rna EXP_ISOFORM_RNA
                        isoforms.fpkm_tracking file from Cufflinks
  --cov-site-dna COV_SITE_DNA
                        bam-readcount output file for tumor DNA BAM and snvs
  --cov-site-rna COV_SITE_RNA
                        bam-readcount output file for tumor RNA BAM and snvs



##FASTQ mode:

usage: Fastq2Mut.py [-h] [-r REFERENCE] [-i KNOWN_INDELS] [-s KNOWN_SNPS]
                    [-d BINDIR] [-p SAMPLEID]
                    [--snv-calling-method SNV_CALLING_METHOD]
                    readsFilesCase readsFilesCtrl

NeoEpitome-Fastq2Mut (v1.0)

positional arguments:
  readsFilesCase        Tumor sample paired-end fastq files seperated by ",".
  readsFilesCtrl        Tumor sample paired-end fastq files seperated by ",".

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        Reference genome used for short reads mapping and GATK
                        realigments.
  -i KNOWN_INDELS, --known-indels KNOWN_INDELS
                        known indel for GATK use. Sorted VCF format. Detail
                        see GATK.
  -s KNOWN_SNPS, --known-snps KNOWN_SNPS
                        dbSNP for GATK use. Sorted VCF format. Detail see
                        GATK.
  -d BINDIR, --binDir BINDIR
                        Directory for java applications indluding GATK,
                        picard.
  -p SAMPLEID, --sampleID SAMPLEID
                        Sample ID will be used for output folder name and
                        reads group name.
  --snv-calling-method SNV_CALLING_METHOD
                        SNV calling method can be set as VarScan2, MuTect,
                        MuSE, or consensus of the three. Default is VarScan2.

 
 usage: Fastq2Exp.py [-h] [--starGenomeDir STARGENOMEDIR] [--gtf GTF]
                    [-p SAMPLEID]
                    readsFilesCaseRNA

NeoEpitome-Fastq2Mut (v1.0)

positional arguments:
  readsFilesCaseRNA     Tumor sample paired-end fastq files seperated by ",".

optional arguments:
  -h, --help            show this help message and exit
  --starGenomeDir STARGENOMEDIR
                        Reference genome used for short reads mapping and GATK
                        realigments.
  --gtf GTF
  -p SAMPLEID, --sampleID SAMPLEID
                        Sample ID will be used for output folder name and
                        reads group name.


usage: BAM2Cov.py [-h] [-p SAMPLEID] [-r REFERENCE]
                  [--tumor-dna-bam TUMOR_DNA_BAM]
                  [--tumor-rna-bam TUMOR_RNA_BAM]
                  vcf_input

NeoEpitome-Mut2Antigen (v1.0)

positional arguments:
  vcf_input             input annotated somatic mutation VCF file.

optional arguments:
  -h, --help            show this help message and exit
  -p SAMPLEID, --sampleID SAMPLEID
                        Sample ID will be used for output folder name and
                        reads group name.
  -r REFERENCE, --reference REFERENCE
                        Reference genome used for short reads mapping and GATK
                        realigments.
  --tumor-dna-bam TUMOR_DNA_BAM
                        BAM file for tumor DNA seq.
  --tumor-rna-bam TUMOR_RNA_BAM
                        BAM file for tumor RNA seq.


                    
Output format:
