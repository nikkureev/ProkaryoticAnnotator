def fastaParser(infile):
    seqs = []
    headers = []
    with open(infile, 'r') as f:
        sequence = ""
        header = None
        for line in f:
            if line.startswith('>'):
                headers.append(line[1:-1])
                if header:
                    seqs.append(sequence)
                sequence = ""
                header = line[1:]
            else:
                sequence += line.rstrip()
        seqs.append(sequence)
    return headers, seqs


def rev_comp(sequence):
    rev_seq = ''
    for nucs in reversed(sequence):
        if nucs == 'A':
            rev_seq += 'T'
        elif nucs == 'T':
            rev_seq += 'A'
        elif nucs == 'G':
            rev_seq += 'C'
        elif nucs == 'C':
            rev_seq += 'G'
    return rev_seq


def loop_finder(sequence, k):
    loop_count = 0
    for var, j in zip(list([sequence[i : i + k] for i in range(len(sequence) - k)]), list(range(k, len(sequence)))):
        if rev_comp(var) in sequence[j + 1 : len(sequence)]:
            loop_count += 1
    return loop_count

