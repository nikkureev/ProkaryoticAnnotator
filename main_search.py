from prom_term import *
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('-in_gff', '--input_gff_file', help='input gff file (annotation)')
parser.add_argument('-in_fasta', '--input_fasta_file', help='input fasta file')
parser.add_argument('-out', '--output_file', help='output directory')
parser.add_argument('-wd', '--working_directory',
                    help='set working directory (without "/" at the end), if None current directory will be used')

args = parser.parse_args()
input_gff = args.input_gff_file
input_fasta = args.input_fasta_file
prom_term_output_file = args.output_file
output_file = args.output_file
working_directory = os.path.abspath(os.curdir)
if args.working_directory:
    working_directory = args.working_directory


def main(input_gff, out_gff, input_fasta, inter_gene_file, prom_term_output_file):
    parsing_gff(input_gff, out_gff)
    take_sequences(out_gff, input_fasta, inter_gene_file)
    names, seqs = fastaParser(inter_gene_file)
    with open(prom_term_output_file, 'w') as f:
        f.write('seq_name ' + 'sequence ' + 'prom_count ' + 'term_count\n')
        for name, seq in zip(names, seqs):
            prom_count = promoter_identifier(seq, name.split()[1].split(':')[1])
            term_count = terminator_identifier(seq, name.split()[1].split(':')[1])
            f.write(str(name.split()[0]) + ' ' + str(prom_count) + ' ' + str(term_count) + '\n')


if __name__ == '__main__':
    main(input_gff=input_gff,
         out_gff=working_directory + 'genes.txt',
         input_fasta=input_fasta,
         inter_gene_file=working_directory + 'inter_genes.txt',
         prom_term_output_file=prom_term_output_file)


# files
# 'C:/AGlab/PROKKA_09172020.gff', 'C:/AGlab/Genes.txt', 'C:/AGlab/NC_000913.fna', 'C:/AGlab/inter_genes.txt',
#          'C:/AGlab/out.txt'
