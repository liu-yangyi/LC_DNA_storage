import gzip
from typing import Union
from Bio import SeqIO
import os


def load_genome(file_path: str) -> str:
    """
    加载.fasta/.fsa/.gz基因组文件，返回纯DNA序列字符串（只含A/T/C/G/N）

    参数:
        file_path: str - 基因组文件路径

    返回:
        dna_sequence: str - 纯碱基序列，所有行拼接，去掉描述头
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"文件未找到: {file_path}")

    is_gz = file_path.endswith(".gz")

    # 打开文件
    open_func = gzip.open if is_gz else open

    with open_func(file_path, "rt") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        if not records:
            raise ValueError("未能解析出FASTA格式记录，请检查输入文件格式")

        # 多个contig时拼接所有序列
        sequence = "".join(str(record.seq).upper() for record in records)

    # 保留合法碱基字符
    allowed_bases = {"A", "T", "C", "G", "N"}
    filtered_seq = "".join(base for base in sequence if base in allowed_bases)

    return filtered_seq
