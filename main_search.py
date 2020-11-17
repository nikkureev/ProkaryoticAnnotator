from prom_term import *
import argparse
import os


# parser = argparse.ArgumentParser()
# parser.add_argument('-in_gff', '--input_gff_file', help='input gff file (annotation)')
# parser.add_argument('-in_fasta', '--input_fasta_file', help='input fasta file')
# parser.add_argument('-out', '--output_file', help='output directory')
# parser.add_argument('-wd', '--working_directory',
#                     help='set working directory (without "/" at the end), if None current directory will be used')
#
# args = parser.parse_args()
# input_fasta_file = args.input
# output_file = args.output_file
# input_error, save_index_file, save_contig_file, save_logging_file = False, False, False, False
# working_directory = os.path.abspath(os.curdir)
# if args.working_directory:
#     working_directory = args.working_directory


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
    main('C:/AGlab/PROKKA_09172020.gff', 'C:/AGlab/Genes.txt', 'C:/AGlab/NC_000913.fna', 'C:/AGlab/inter_genes.txt',
         'C:/AGlab/out.txt')


# files
# 'C:/AGlab/Genes.txt', 'C:/AGlab/NC_000913.fna', 'C:/AGlab/inter_genes.txt'