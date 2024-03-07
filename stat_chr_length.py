from Bio import SeqIO

import sys
fa_file = sys.argv[1] ## genome file

records = SeqIO.parse(fa_file, "fasta")

# 统计每条染色体的长度
chromosome_lengths = {}
for record in records:
    chromosome_name = record.id
    chromosome_length = len(record)
    chromosome_lengths[chromosome_name] = chromosome_length

# 输出结果
for chromosome_name, chromosome_length in chromosome_lengths.items():
    print(f"{chromosome_name}\t1\t{chromosome_length}")
