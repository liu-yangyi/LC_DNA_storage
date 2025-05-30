import pickle
import json
import os
from bitarray import bitarray
from reedsolo import RSCodec
from typing import List, Dict
import math

def lzw_compress(data: bytes, max_dict_size=2**20) -> List[int]:
    dict_size = 256
    dictionary = {bytes([i]): i for i in range(dict_size)}
    w = b""
    result = []

    for c in data:
        wc = w + bytes([c])
        if wc in dictionary:
            w = wc
        else:
            result.append(dictionary[w])
            if dict_size < max_dict_size:
                dictionary[wc] = dict_size
                dict_size += 1
            w = bytes([c])
    if w:
        result.append(dictionary[w])
    return result

def load_kmer_index(kmer_to_pos_path):
    with open(kmer_to_pos_path, 'rb') as f:
        kmer_to_pos = pickle.load(f)
    return kmer_to_pos

def build_full_kmer_dict(kmer_to_pos):
    all_kmers = list(kmer_to_pos.keys())
    kmer_id_map = {i: kmer for i, kmer in enumerate(all_kmers)}
    kmer_reverse_map = {v: k for k, v in kmer_id_map.items()}
    return kmer_id_map, kmer_reverse_map

def compress_and_encode(file_path, ecc_symbols=64):
    with open(file_path, 'rb') as f:
        raw_data = f.read()
    compressed = lzw_compress(raw_data)

    max_index = max(compressed) if compressed else 0
    byte_length = max((max_index.bit_length() + 7) // 8, 1)
    byte_length = min(byte_length, 4)

    byte_stream = bytearray()
    for num in compressed:
        byte_stream += num.to_bytes(byte_length, byteorder='big')

    rs = RSCodec(ecc_symbols)
    rs_encoded = rs.encode(byte_stream)

    bits = bitarray()
    bits.frombytes(rs_encoded)
    return bits, byte_length

def bits_to_kmer_sequence(bits: bitarray, kmer_id_map: Dict[int, str], bit_per_kmer: int) -> List[str]:
    kmer_sequence = []
    total_bits = len(bits)
    for i in range(0, total_bits, bit_per_kmer):
        chunk = bits[i:i + bit_per_kmer]
        if len(chunk) < bit_per_kmer:
            chunk += bitarray('0' * (bit_per_kmer - len(chunk)))
        idx = int(chunk.to01(), 2)
        kmer = kmer_id_map.get(idx)
        if not kmer:
            raise ValueError(f"无对应 k-mer ID：{idx}")
        kmer_sequence.append(kmer)
    return kmer_sequence

def encode_file_v2(input_file, kmer_to_pos_path, top_n_kmers=131072, output_prefix="encoded_output", ecc_symbols=64):
    print("🔄 加载 k-mer 索引...")
    kmer_to_pos = load_kmer_index(kmer_to_pos_path)

    print(f"📊 构建 top-{top_n_kmers} 高频 k-mer ID 映射...")
    kmer_id_map, kmer_reverse_map = build_full_kmer_dict(kmer_to_pos)

    print("🧬 LZW 压缩 + RS 编码...")
    bitarray_obj, byte_length = compress_and_encode(input_file, ecc_symbols=ecc_symbols)

    bit_per_kmer = math.ceil(math.log2(top_n_kmers))
    print(f"🔗 将 bitstream 映射为 k-mer 序列 (每组 {bit_per_kmer} bits)...")
    kmer_sequence = bits_to_kmer_sequence(bitarray_obj, kmer_id_map, bit_per_kmer)

    # === 保存输出文件 ===
    json.dump(kmer_sequence, open(f"{output_prefix}_kmer_sequence.json", "w"))

    with open(f"{output_prefix}_kmer_id_map.pkl", "wb") as f:
        pickle.dump(kmer_id_map, f)

    with open(f"{output_prefix}_kmer_reverse_map.pkl", "wb") as f:
        pickle.dump(kmer_reverse_map, f)

    metadata = {
        "byte_length": byte_length,
        "bit_per_kmer": bit_per_kmer,
        "top_n_kmers": top_n_kmers,
        "ecc_symbols": ecc_symbols
    }
    json.dump(metadata, open(f"{output_prefix}_metadata.json", "w"))

    print("✅ 编码完成！生成文件：")
    print(f"  ➤ k-mer 序列: {output_prefix}_kmer_sequence.json")
    print(f"  ➤ k-mer ID 映射: {output_prefix}_kmer_id_map.pkl")
    print(f"  ➤ 元数据: {output_prefix}_metadata.json")
