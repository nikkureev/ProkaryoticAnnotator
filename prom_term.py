from fasta_parse import fastaParser, loop_finder


def promoter_identifier(sequence, strand):
    AT_count = 0
    if strand == '+':
        promoter_region = sequence[-50:-1]
    else:
        promoter_region = sequence[0:50]

    for nucleotide in promoter_region:
        if nucleotide == 'A' or nucleotide == 'T':
            AT_count += 1
    AT_content = AT_count / len(promoter_region)
    if AT_content > 0.5:
        return 1


def terminator_identifier(sequence, strand):
    GC_count = 0
    term_count = 0
    if strand == '+':
        terminator_region = sequence[0:50]
    else:
        terminator_region = sequence[-50:-1]

    for nucleotide in terminator_region:
        if nucleotide == 'G' or nucleotide == 'C':
            GC_count += 1
    AT_content = GC_count / len(terminator_region)
    if AT_content > 0.5:
        term_count -= 1

    term_count -= loop_finder(sequence, 4)
    return term_count


def take_sequences(gene_file, genome_fasta, inter_gene_file):
    start_list, end_list, strand_list = [], [], []
    inter_gene_dict = {}

    with open(gene_file, 'r') as d:
        for lines in d:
            try:
                start_list.append(lines.split()[0])
                end_list.append(lines.split()[1])
                strand_list.append(lines.split()[2])
            except:
                break

    n = 0
    for i in range(len(start_list) - 1):
        name = 'id' + str(n)
        inter_gene_dict[name] = [int(end_list[n]),
                                 int(start_list[n + 1]),
                                 int(start_list[n + 1]) - int(end_list[n]),
                                 strand_list[n + 1]]
        n += 1

    n = 0
    with open(inter_gene_file, 'w') as f:
        source_name, source_sequence = fastaParser(genome_fasta)
        for i, v in inter_gene_dict.items():
            if v[2] >= 100:
                sequence = source_sequence[0][v[0]: v[1]]
                f.write('>inter_gene_id' + str(n) + ' strand:' + v[3] + ' coordinates:' + str(v[0]) + '-' + str(v[1]) + '\n')
                f.write(sequence + '\n')
                n += 1


def parsing_gff(input_file, out_file):
    with open(input_file, 'r') as f:
        with open(out_file, 'w') as o:
            for lines in f:
                if not lines.startswith('#'):
                    o.write(lines.split()[3] + ' ' + lines.split()[4] + ' ' + lines.split()[6] + '\n')



