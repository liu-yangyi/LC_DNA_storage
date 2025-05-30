import pyopencl as cl
import numpy as np
from collections import defaultdict
from typing import Dict, Tuple


# 碱基字符映射为2-bit编码: A=00, C=01, G=10, T=11
base2bits = {
    'A': 0b00,
    'C': 0b01,
    'G': 0b10,
    'T': 0b11
}


def encode_kmer_window(sequence: str, k: int) -> np.ndarray:
    """
    将k-mer窗口转换为无符号整数数组表示（每个k-mer一个整型）
    """
    kmers = []
    bitmask = (1 << (2 * k)) - 1

    current = 0
    count = 0

    for base in sequence:
        if base not in base2bits:
            count = 0
            current = 0
            continue
        current = ((current << 2) | base2bits[base]) & bitmask
        count += 1
        if count >= k:
            kmers.append(current)

    return np.array(kmers, dtype=np.uint64)


def count_kmers_gpu(sequence: str, k: int) -> Dict[int, int]:
    """
    使用OpenCL并行计算所有k-mer出现频率，返回频率字典（键为编码后的整数）
    """
    encoded_kmers = encode_kmer_window(sequence, k)
    total_kmers = len(encoded_kmers)

    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    mf = cl.mem_flags
    input_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=encoded_kmers)
    output_buf = cl.Buffer(ctx, mf.WRITE_ONLY, size=encoded_kmers.nbytes)

    kernel_code = """
    __kernel void count_kmers(__global const ulong *kmers,
                              __global ulong *output,
                              const ulong count) {
        int gid = get_global_id(0);
        if (gid < count) {
            output[gid] = kmers[gid]; // 拷贝原始kmer值，稍后CPU统计频率
        }
    }
    """

    prg = cl.Program(ctx, kernel_code).build()
    prg.count_kmers(queue, (total_kmers,), None, input_buf, output_buf, np.uint64(total_kmers))

    output = np.empty_like(encoded_kmers)
    cl.enqueue_copy(queue, output, output_buf)

    freq = defaultdict(int)
    for kmer in output:
        freq[kmer] += 1

    return dict(freq)


def decode_kmer(encoded: int, k: int) -> str:
    """将编码后的整数还原为ATCG序列"""
    bits2base = ['A', 'C', 'G', 'T']
    kmer = []
    for _ in range(k):
        kmer.append(bits2base[encoded & 0b11])
        encoded >>= 2
    return ''.join(reversed(kmer))
