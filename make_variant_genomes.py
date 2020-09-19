"""This short script uses samtools, picard and GATK to produce pseudogenomes"""
import os

samtools = '/opt/Bio/samtools/1.9/bin/samtools faidx'
picard = '/opt/Bio/picard-tools/2.9.0/bin/picard.sh CreateSequenceDictionary'
gatk = '/opt/Bio/GenomeAnalysisTK/3.8-0/bin/GenomeAnalysisTK.sh -T FastaAlternateReferenceMaker'
if not os.path.isfile('Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai'):
    os.system('{} Arabidopsis_thaliana.TAIR10.dna.toplevel.fa'.format(samtools))
if not os.path.isfile('Arabidopsis_thaliana.TAIR10.dna.toplevel.dict'):
    os.system("{} R=Arabidopsis_thaliana.TAIR10.dna.toplevel.fa O= Arabidopsis_thaliana.TAIR10.dna.toplevel.dict".format(picard))

for vcf in os.listdir('/nam-99/ablage/nam/peleke/snp_vcf'):   
    variant = '/nam-99/ablage/nam/peleke/snp_vcf/' + vcf
    ecotype_id = vcf.split('_')[1].split('.')[0]
    file = '/nam-99/ablage/nam/peleke/snp_only_pseudogenomes/' + ecotype_id + '.fa'
    if not os.path.exists(file):
        os.system('{} -R Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -o {} -V {}'.format(gatk, file, variant))
