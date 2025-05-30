from Bio import SeqIO
import gzip

with gzip.open("mutated_genome.fsa.gz", "rt") as in_handle, open("mutated_genome.fasta", "w") as out_handle:
    records = SeqIO.parse(in_handle, "fastq")
    count = SeqIO.write(records, out_handle, "fasta")

print(f"成功将 {count} 条记录写入 mutated_genome.fasta")
