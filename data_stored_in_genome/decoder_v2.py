import os
import json
import gzip
import pickle
from bitarray import bitarray
from reedsolo import RSCodec
from collections import defaultdict
from typing import Optional

def build_kmer_index_from_fasta(fasta_gz_path, k=11):
    kmer_to_pos = defaultdict(list)
    print(f"✨ 构建 k-mer 索引 (基因经约 {k})...")
    with gzip.open(fasta_gz_path, 'rt') as f:
        sequence = ''
        for line in f:
            if line.startswith('>'):
                continue
            sequence += line.strip().upper()

    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if 'N' not in kmer:
            kmer_to_pos[kmer].append(i)
    print(f"✅ 完成: {len(kmer_to_pos)} 个唯一 k-mers")
    return kmer_to_pos

def bits_to_bytes(bits):
    b = bytearray()
    b.extend(bits.tobytes())
    return b

def lzw_decompress(compressed_data: bytes, byte_length: int) -> bytes:
    dict_size = 256
    dictionary = {i: bytes([i]) for i in range(dict_size)}
    result = bytearray()

    # 分解编码
    data_len = len(compressed_data)
    codes = []
    for i in range(0, data_len, byte_length):
        index = int.from_bytes(compressed_data[i:i + byte_length], 'big')
        codes.append(index)

    if not codes:
        return b''

    w = dictionary[codes[0]]
    result.extend(w)
    for k in codes[1:]:
        if k in dictionary:
            entry = dictionary[k]
        elif k == dict_size:
            entry = w + w[:1]
        else:
            raise ValueError(f"Invalid LZW index: {k}")
        result.extend(entry)
        dictionary[dict_size] = w + entry[:1]
        dict_size += 1
        w = entry
    return bytes(result)

def decode_file(input_prefix: str, output_path: str, mutated_genome_path: Optional[str] = None):
    print("🔄 加载元数据...")
    with open(f"{input_prefix}_metadata.json", "r") as f:
        metadata = json.load(f)

    byte_length = metadata.get("byte_length")
    if byte_length is None:
        print("⚠️ metadata 中缺失 byte_length, 无法解码")
        return

    # 加载 k-mer 映射表
    with open(f"{input_prefix}_kmer_reverse_map.pkl", "rb") as f:
        kmer_reverse_map = pickle.load(f)

    reverse_map = kmer_reverse_map  # kmer => id
    id_map = {v: k for k, v in reverse_map.items()}  # id => kmer

    top_n_kmers = len(id_map)
    bits_per_entry = metadata.get("bits_per_entry") or (top_n_kmers - 1).bit_length()
    print(f"✅ 推断 bits_per_entry = {bits_per_entry}")

    # 加载编码 kmer 序列
    with open(f"{input_prefix}_kmer_sequence.json", "r") as f:
        kmer_sequence = json.load(f)

    # 构建 k-mer => 编码值 map
    kmer_to_id = {kmer: idx for idx, kmer in id_map.items()}

    print("🐍 解码 k-mer 序列为低位段 bitstream...")
    bits = bitarray()
    missing_kmers = 0
    for kmer in kmer_sequence:
        if kmer in kmer_to_id:
            idx = kmer_to_id[kmer]
            bits.extend(f"{{idx:0{{bits_per_entry}}b}}".format(idx=idx, bits_per_entry=bits_per_entry))
        else:
            missing_kmers += 1
            bits.extend('0' * bits_per_entry)  # 用 0 填充

    print(f"❗ 未找到的 k-mers: {missing_kmers} / {len(kmer_sequence)}")
    print("💡 转换 bitarray 为 bytes...")
    corrected = bits_to_bytes(bits)

    print("🔧 使用 Reed-Solomon 进行错误修复...")
    rs = RSCodec(64)
    corrected_bytes = bytes(corrected)
    decoded = rs.decode(corrected_bytes)[0]  # 返回元素 part

    print("📦 扫描 LZW 解压...")
    decompressed = lzw_decompress(decoded, byte_length)

    with open(output_path, 'wb') as f:
        f.write(decompressed)
    print(f"✅ 解码完成，输出保存到: {output_path}")

    if mutated_genome_path:
        print(f"📘 如需比较连代后的基因经差异请使用: {mutated_genome_path}")
