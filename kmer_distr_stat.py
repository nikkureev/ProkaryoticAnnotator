from itertools import product
import matplotlib.pyplot as plt


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


core = 5
c = product('ATGC', repeat=core)
kmers = []
for i in c:
    kmers.append(i)

file = 'C:/AGlab/tandem_2000.fasta'

headers, sequences = fastaParser(file)

d = {}
for h, s in zip(headers, sequences):
    seq_list = list(s)
    for kmer in kmers:
        for i in range(len(seq_list) - len(kmer)):
            if tuple(seq_list[i: i+len(kmer)]) == kmer:
                if kmer in d.keys():
                    d[kmer].append((i, i+len(kmer)))
                else:
                    d[kmer] = [(i, i+len(kmer))]


d_k = {}
for i, v in d.items():
    d_k[i] = list([0 for p in range(100)])
    for k in range(100):
        for j in v:
            if j[0] < k < j[1]:
                d_k[i][k] += 1
for i, v in d_k.items():
    plt.plot(list(range(len(v))), v, label=i)
    n = 0
    for k in v:
        if k > 120:
            print(''.join(i), n, '-' + str(100 - n))
        n += 1
plt.show()
    

