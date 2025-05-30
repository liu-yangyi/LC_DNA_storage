import pickle
from collections import defaultdict
from typing import Dict, List


def build_bidirectional_index(genome_seq: str, k: int):
    kmer_to_positions: Dict[str, List[int]] = defaultdict(list)
    position_to_kmer: Dict[int, str] = {}

    for i in range(len(genome_seq) - k + 1):
        kmer = genome_seq[i:i+k]
        if 'N' in kmer:
            continue
        kmer_to_positions[kmer].append(i)
        position_to_kmer[i] = kmer

    return kmer_to_positions, position_to_kmer


def save_index(kmer_to_positions: Dict[str, List[int]],
               position_to_kmer: Dict[int, str],
               k: int):
    with open(f"kmer_to_positions_k{k}.pkl", "wb") as f1:
        pickle.dump(kmer_to_positions, f1)

    with open(f"position_to_kmer_k{k}.pkl", "wb") as f2:
        pickle.dump(position_to_kmer, f2)

    print(f"✅ 索引已保存：kmer_to_positions_k{k}.pkl, position_to_kmer_k{k}.pkl")
