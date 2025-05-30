from typing import List, Tuple, Dict
from kmer_gpu_counter import count_kmers_gpu


def evaluate_kmer_stats(freq_dict: Dict[int, int]) -> Tuple[int, int, float]:
    """
    计算给定k值下k-mer频率的基本统计信息：
    - 唯一k-mer数量
    - 高频k-mer（出现次数 > 1）数量
    - Shannon熵作为多样性评估指标
    """
    from math import log2

    total = sum(freq_dict.values())
    unique = len(freq_dict)
    high_freq = sum(1 for v in freq_dict.values() if v > 1)

    shannon_entropy = -sum((v / total) * log2(v / total) for v in freq_dict.values() if v > 0)

    return unique, high_freq, shannon_entropy


def find_best_k(genome_seq: str, k_range: Tuple[int, int]) -> Tuple[int, Dict[int, Dict]]:
    """
    在给定k值范围内，自动选择最优k，并返回每个k值对应的统计结果
    """
    k_stat_map = {}

    for k in range(k_range[0], k_range[1] + 1):
        print(f"🔍 分析 k = {k} ...")
        freq_dict = count_kmers_gpu(genome_seq, k)
        unique, high_freq, entropy = evaluate_kmer_stats(freq_dict)
        k_stat_map[k] = {
            "unique_kmers": unique,
            "high_freq_kmers": high_freq,
            "shannon_entropy": entropy
        }

    # 选择高频多、熵高、重复少的综合最优k
    best_k = max(k_stat_map.items(), key=lambda item: (
        item[1]['high_freq_kmers'], item[1]['shannon_entropy'], item[1]['unique_kmers']
    ))[0]

    return best_k, k_stat_map
