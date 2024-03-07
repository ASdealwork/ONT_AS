import sys


def count_aaaccct_in_fasta_file(file_path):
    with open(file_path, 'r') as f:
        fasta_sequences = f.read().split('>')[1:]
        for sequence in fasta_sequences:
            sequence_lines = sequence.split('\n')
            header = sequence_lines[0]
            sequence = ''.join(sequence_lines[1:])
            CCCTAAA_count = sequence.count('CCCTAAA')
            TTTAGGG_count = sequence.count('TTTAGGG')
            CCCTAAA_TTTAGGG_count = int(CCCTAAA_count)+int(TTTAGGG_count)
            print(f"{header}\tCCCTAAA\t{CCCTAAA_count}\tTTTAGGG\t{TTTAGGG_count}\tCCCTAAA_TTTAGGG\t{CCCTAAA_TTTAGGG_count}")


file_path = sys.argv[1] ## fasta file that needs to be processed
# 调用函数统计每条fasta序列中AAACCCT的数量
count_aaaccct_in_fasta_file(file_path)
