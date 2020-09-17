"""
This file processes the all_vs_all protein blast output
Download maize proteome from:
 ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/zea_mays/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
Gunzip Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
create database on localhost then blast Arabidopsis_thaliana.TAIR10.pep.all.fa  vs database :
 makeblastdb -in Arabidopsis_thaliana.TAIR10.pep.all.fa -dbtype prot -title A_thaliana_TAIR10 -parse_seqids -hash_index -out A_thaliana_TAIR10
 blastp -db A_thaliana_TAIR10 -query Arabidopsis_thaliana.TAIR10.pep.all.fa -out result_TAIR10 -outfmt 6
"""
import pandas as pd
from Bio import SeqIO
pd.options.display.width = 0

path_to_A_thalianaprot = '/nam-99/ablage/nam/peleke/Arabidopsis_thaliana.TAIR10.pep.all.fa'
path_to_blast_result = '/nam-99/ablage/nam/peleke/result_TAIR10'
write_path = '/nam-99/ablage/nam/peleke/result_TAIR10_filtered.csv'
cols = ['qseqid', 'sseqid', 'pident',
        'length', 'mismatch', 'gapopen', 'qstat', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

blast_result = pd.read_csv(path_to_blast_result, delimiter='\t', names=cols)

gene_list = []
protein_list = []
for rec in SeqIO.parse(path_to_A_thalianaprot, 'fasta'):
    prot_id = rec.id
    gene_id = rec.description.split(' ')[3].split(':')[1]
    protein_list.append(prot_id)
    gene_list.append(gene_id)

info = pd.DataFrame({'gene_list': gene_list, 'protein_list': protein_list}, index=protein_list)

qgene = info.loc[blast_result.qseqid.values.tolist(), ['gene_list']].values
sgene = info.loc[blast_result.sseqid.values.tolist(), ['gene_list']].values
blast_result['qgene'] = qgene
blast_result['sgene'] = sgene

pairs = [i[12] + '_' + i[13]
         if i[12] > i[13]
         else i[13] + '_' + i[12]
         for i in blast_result.values]
blast_result['pairs'] = pairs

blast_result = blast_result[blast_result.evalue < 0.001]
blast_result = blast_result[blast_result.bitscore > 50]
blast_result = blast_result[blast_result.qgene != blast_result.sgene]

blast_result = blast_result.drop_duplicates(subset='pairs', keep='first')
blast_result = blast_result.sort_values(by=['pairs', 'bitscore'], ascending=True)

blast_result.to_csv(write_path, index=False)
