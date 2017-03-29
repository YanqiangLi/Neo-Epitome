Dependency:
1. Python 2.7 (biopython, numpy)
2. R
3. BWA,Samtools,
3. IEDB (optional)

VCF mode:
usage: Mut2Antigen.py [-h] [-j JUNCTION_INPUT] [-e EPITOPE_LEN_LIST]
                      [-a HLA_ALLELE_LIST] [--step STEP]
                      [--ic50-cut-off IC50_CUT_OFF] [-o OUTDIR]
                      [-p PROTEIN_MS] [--protein-reference PROTEIN_REFERENCE]
                      [--protein-mod PROTEIN_MOD]
                      [--exp-gene-rna EXP_GENE_RNA]
                      [--exp-isoform-rna EXP_ISOFORM_RNA]
                      [--cov-site-dna COV_SITE_DNA]
                      [--cov-site-rna COV_SITE_RNA]
                      vcf_input


FASTQ mode:
usage: Fastq2Mut.py [-h] [-r REFERENCE] [-i KNOWN_INDELS] [-s KNOWN_SNPS]
                    [-d BINDIR] [-p SAMPLEID] [--starGenomeDir STARGENOMEDIR] [--gtf GTF]
                    [-p SAMPLEID] [--readsFilesCaseRNA] [-e EPITOPE_LEN_LIST]
                      [-a HLA_ALLELE_LIST] [--step STEP]
                      [--ic50-cut-off IC50_CUT_OFF] [-o OUTDIR]
                      [-p PROTEIN_MS] [--protein-reference PROTEIN_REFERENCE]
                      [--protein-mod PROTEIN_MOD]
                    readsFilesCase readsFilesCtrl
                    
                    
Output format:
