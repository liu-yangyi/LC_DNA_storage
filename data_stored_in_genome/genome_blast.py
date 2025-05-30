from Bio import SeqIO

def extract_kmers_from_fasta(fasta_file, k):
    kmers = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if 'N' not in kmer:
                kmers.add(kmer)
    return kmers

original_kmers = extract_kmers_from_fasta('JRYM01.1.fsa_nt', k=11)
mutated_kmers = extract_kmers_from_fasta('mutated_genome.fasta', k=11)

shared = original_kmers & mutated_kmers
print(f"原始: {len(original_kmers)} 个 k-mers")
print(f"突变后: {len(mutated_kmers)} 个 k-mers")
print(f"保留: {len(shared)} 个 k-mers ({len(shared)/len(original_kmers):.2%})")

import json
from Bio import SeqIO

def extract_kmers_from_fasta(fasta_path, k):
    kmers = set()
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if 'N' not in kmer:
                kmers.add(kmer)
    return kmers

# 加载编码生成的 k-mer 序列
encoded_kmers = json.load(open("encoded_output/shakespeare_kmer_sequence.json"))
encoded_kmers_set = set(encoded_kmers)

# 提取突变基因组中的 k-mers
mutated_kmers = extract_kmers_from_fasta("mutated_genome.fasta", k=len(encoded_kmers[0]))

# 比较
shared = encoded_kmers_set & mutated_kmers
print(f"编码使用的k-mers数量: {len(encoded_kmers_set)}")
print(f"突变后仍存在的k-mers: {len(shared)}")
print(f"保留率: {len(shared)/len(encoded_kmers_set):.2%}")

